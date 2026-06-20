// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use flate2::read::MultiGzDecoder;
use memmap2::Mmap;
use std::{
    fs::File,
    io::{BufRead, BufReader, Cursor, Read, Seek, Write},
    path::Path,
};
use twobit::TwoBitFile;

use crate::{
    cli::is_stdio_path,
    error::{GenomeMaskError, Result},
};

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const TWOBIT_MAGIC: [u8; 4] = [0x43, 0x27, 0x41, 0x1a];
const TWOBIT_REV_MAGIC: [u8; 4] = [0x1a, 0x41, 0x27, 0x43];

/// Statistics from processing a genome file.
#[derive(Debug, Default, Clone, Copy, Eq, PartialEq)]
pub struct ProcessStats {
    /// Number of FASTA records processed
    pub records: u64,
    /// Number of masking events applied
    pub events: u64,
}

impl ProcessStats {
    fn add_record(&mut self, events: u64) {
        self.records += 1;
        self.events += events;
    }
}

/// Trait for transforming genome records during processing.
///
/// Implement this trait to create custom masking strategies.
pub trait RecordTransformer {
    /// Transforms a single FASTA record's sequence.
    ///
    /// # Arguments
    /// * `header` - The FASTA header (without >)
    /// * `sequence` - The sequence data (mutated in place)
    /// * `record_index` - Index of the current record
    ///
    /// # Returns
    /// * `Ok(u64)` - Number of events (bases masked/replaced)
    fn transform_record(
        &mut self,
        header: &[u8],
        sequence: &mut Vec<u8>,
        record_index: u64,
    ) -> Result<u64>;

    /// Called after all records have been processed.
    ///
    /// Use to validate that all expected data was processed.
    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

/// Input format detected from file magic bytes.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum InputFormat {
    /// Plain FASTA format
    Fasta,
    /// Gzip-compressed FASTA
    GzipFasta,
    /// 2bit binary format
    TwoBit,
}

/// Processes an input genome file and writes transformed FASTA records.
///
/// Supports FASTA, gzipped FASTA, and 2bit formats.
/// Automatically detects format from file magic bytes.
///
/// # Arguments
/// * `input` - Input file path or "-" for stdin
/// * `writer` - Output writer for FASTA records
/// * `transformer` - Record transformer for masking
///
/// # Returns
/// * `Ok(ProcessStats)` - Statistics on records and events
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_input_to_fasta(Path::new("genome.fa"), &mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
pub fn process_input_to_fasta<T: RecordTransformer, W: Write>(
    input: &Path,
    writer: &mut W,
    transformer: &mut T,
    preserve_mask: bool,
) -> Result<ProcessStats> {
    let stats = if is_stdio_path(input) {
        process_stdin_to_fasta(writer, transformer, preserve_mask)?
    } else {
        process_file_to_fasta(input, writer, transformer, preserve_mask)?
    };

    transformer.finish()?;
    Ok(stats)
}

/// Returns true when a filesystem input path is a 2bit file.
pub fn is_twobit_input_path(input: &Path) -> Result<bool> {
    if is_stdio_path(input) {
        return Ok(false);
    }

    let probe = read_probe(input)?;
    Ok(matches!(sniff_input(&probe), Some(InputFormat::TwoBit)))
}

/// Processes stdin to FASTA format.
///
/// Reads stdin and processes each record.
///
/// # Arguments
/// * `writer` - Output writer for FASTA records
/// * `transformer` - Record transformer for masking
///
/// # Returns
/// * `Ok(ProcessStats)` - Statistics on records and events
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_stdin_to_fasta(&mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
fn process_stdin_to_fasta<T: RecordTransformer, W: Write>(
    writer: &mut W,
    transformer: &mut T,
    preserve_mask: bool,
) -> Result<ProcessStats> {
    let stdin = std::io::stdin();
    let mut reader = BufReader::new(stdin.lock());
    let format = {
        let prefix = reader
            .fill_buf()
            .map_err(|err| GenomeMaskError::io("cannot inspect stdin", err))?;
        if prefix.is_empty() {
            return Err(GenomeMaskError::EmptyInput);
        }
        sniff_input(prefix).ok_or_else(|| {
            GenomeMaskError::UnsupportedInput(
                "unsupported stdin format; expected FASTA, gzipped FASTA, or 2bit".to_string(),
            )
        })?
    };

    match format {
        InputFormat::Fasta => process_fasta_reader(reader, writer, transformer, preserve_mask),
        InputFormat::GzipFasta => {
            let decoder = MultiGzDecoder::new(reader);
            process_fasta_reader(BufReader::new(decoder), writer, transformer, preserve_mask)
        }
        InputFormat::TwoBit => {
            let mut data = Vec::new();
            reader
                .read_to_end(&mut data)
                .map_err(|err| GenomeMaskError::io("cannot read 2bit stdin", err))?;
            let genome = TwoBitFile::from_buf(data)
                .map_err(|err| {
                    GenomeMaskError::InvalidTwoBit(format!("cannot read 2bit stdin: {err}"))
                })?
                .enable_softmask(true);
            process_twobit_reader(genome, writer, transformer, preserve_mask)
        }
    }
}

/// Process a file to FASTA format.
///
/// Reads the input file and processes each record.
/// Uses a probe to determine the input format.
///
/// # Arguments
/// * `input` - Path to the input file
/// * `writer` - Output writer for FASTA records
/// * `transformer` - Record transformer for masking
///
/// # Returns
/// * `Ok(ProcessStats)` - Statistics on records and events
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_file_to_fasta(Path::new("genome.fa"), &mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
fn process_file_to_fasta<T: RecordTransformer, W: Write>(
    input: &Path,
    writer: &mut W,
    transformer: &mut T,
    preserve_mask: bool,
) -> Result<ProcessStats> {
    let probe = read_probe(input)?;
    if probe.is_empty() {
        return Err(GenomeMaskError::EmptyInput);
    }

    match sniff_input(&probe) {
        Some(InputFormat::GzipFasta) => {
            let file = File::open(input).map_err(|err| {
                GenomeMaskError::io(format!("cannot open {}", input.display()), err)
            })?;
            let decoder = MultiGzDecoder::new(file);
            process_fasta_reader(BufReader::new(decoder), writer, transformer, preserve_mask)
        }
        Some(InputFormat::TwoBit) => {
            let genome = TwoBitFile::open(input)
                .map_err(|err| {
                    GenomeMaskError::InvalidTwoBit(format!(
                        "cannot open 2bit file {}: {err}",
                        input.display()
                    ))
                })?
                .enable_softmask(true);
            process_twobit_reader(genome, writer, transformer, preserve_mask)
        }
        Some(InputFormat::Fasta) | None => {
            let file = File::open(input).map_err(|err| {
                GenomeMaskError::io(format!("cannot open {}", input.display()), err)
            })?;
            let mmap = unsafe { Mmap::map(&file) }.map_err(|err| {
                GenomeMaskError::io(format!("cannot mmap {}", input.display()), err)
            })?;
            match sniff_input(&mmap[..]) {
                Some(InputFormat::Fasta) => process_fasta_reader(
                    BufReader::new(Cursor::new(&mmap[..])),
                    writer,
                    transformer,
                    preserve_mask,
                ),
                Some(InputFormat::GzipFasta) => Err(GenomeMaskError::UnsupportedInput(format!(
                    "unexpected gzipped data in plain FASTA path {}",
                    input.display()
                ))),
                Some(InputFormat::TwoBit) => Err(GenomeMaskError::UnsupportedInput(format!(
                    "unexpected 2bit data in plain FASTA path {}",
                    input.display()
                ))),
                None => Err(GenomeMaskError::UnsupportedInput(format!(
                    "unsupported input format for {}",
                    input.display()
                ))),
            }
        }
    }
}

/// Processes a FASTA reader to FASTA format.
///
/// Reads the input reader and processes each record.
///
/// # Arguments
/// * `reader` - Buffered reader for the input file
/// * `writer` - Output writer for FASTA records
/// * `transformer` - Record transformer for masking
///
/// # Returns
/// * `Ok(ProcessStats)` - Statistics on records and events
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_fasta_reader(BufReader::new(File::open("genome.fa")?), &mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
fn process_fasta_reader<R: BufRead, T: RecordTransformer, W: Write>(
    mut reader: R,
    writer: &mut W,
    transformer: &mut T,
    preserve_mask: bool,
) -> Result<ProcessStats> {
    let mut line = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut sequence = Vec::new();
    let mut stats = ProcessStats::default();

    loop {
        line.clear();
        let bytes_read = reader
            .read_until(b'\n', &mut line)
            .map_err(|err| GenomeMaskError::io("cannot read FASTA input", err))?;

        if bytes_read == 0 {
            break;
        }

        trim_line_endings(&mut line);
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if let Some(previous_header) = header.replace(line[1..].to_vec()) {
                emit_record(
                    previous_header,
                    &mut sequence,
                    writer,
                    transformer,
                    &mut stats,
                    preserve_mask,
                )?;
            }
        } else {
            let current_header = header.as_ref().ok_or_else(|| {
                GenomeMaskError::InvalidFasta(
                    "sequence data encountered before the first FASTA header".to_string(),
                )
            })?;

            if line.iter().any(|byte| byte.is_ascii_whitespace()) {
                return Err(GenomeMaskError::InvalidFasta(format!(
                    "whitespace found inside sequence for FASTA record '{}'",
                    header_display(current_header)
                )));
            }

            sequence.extend_from_slice(&line);
        }
    }

    match header {
        Some(last_header) => {
            emit_record(
                last_header,
                &mut sequence,
                writer,
                transformer,
                &mut stats,
                preserve_mask,
            )?;
            Ok(stats)
        }
        None => Err(GenomeMaskError::EmptyInput),
    }
}

/// Processes a 2bit reader to FASTA format.
///
/// Reads the input reader and processes each record.
///
/// # Arguments
/// * `genome` - 2bit reader for the input file
/// * `writer` - Output writer for FASTA records
/// * `transformer` - Record transformer for masking
///
/// # Returns
/// * `Ok(ProcessStats)` - Statistics on records and events
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_twobit_reader(TwoBitFile::open("genome.2bit")?, &mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
fn process_twobit_reader<R: Read + Seek, T: RecordTransformer, W: Write>(
    mut genome: TwoBitFile<R>,
    writer: &mut W,
    transformer: &mut T,
    preserve_mask: bool,
) -> Result<ProcessStats> {
    let chrom_names = genome.chrom_names();
    if chrom_names.is_empty() {
        return Err(GenomeMaskError::EmptyInput);
    }

    let mut stats = ProcessStats::default();
    for chrom_name in chrom_names {
        let header = chrom_name.as_bytes().to_vec();
        let mut sequence = genome
            .read_sequence(&chrom_name, ..)
            .map_err(|err| {
                GenomeMaskError::InvalidTwoBit(format!(
                    "cannot read 2bit sequence '{}': {err}",
                    chrom_name
                ))
            })?
            .into_bytes();

        let events = transformer.transform_record(&header, &mut sequence, stats.records)?;
        write_fasta_record(writer, &header, &sequence, preserve_mask)?;
        stats.add_record(events);
    }

    Ok(stats)
}

/// Emits a single FASTA record to an output writer.
///
/// # Arguments
/// * `header` - FASTA header (without >)
/// * `sequence` - Sequence data
/// * `writer` - Output writer
/// * `transformer` - Record transformer for masking
/// * `stats` - Process statistics
///
/// # Example
///
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let stats = process_fasta_reader(BufReader::new(File::open("genome.fa")?), &mut writer, &mut transformer)?;
/// println!("Processed {} records with {} events", stats.records, stats.events);
/// ```
fn emit_record<T: RecordTransformer, W: Write>(
    header: Vec<u8>,
    sequence: &mut Vec<u8>,
    writer: &mut W,
    transformer: &mut T,
    stats: &mut ProcessStats,
    preserve_mask: bool,
) -> Result<()> {
    let events = transformer.transform_record(&header, sequence, stats.records)?;
    write_fasta_record(writer, &header, sequence, preserve_mask)?;
    sequence.clear();
    stats.add_record(events);
    Ok(())
}

/// Writes a single FASTA record to an output writer.
///
/// # Arguments
/// * `writer` - Output writer
/// * `header` - FASTA header (without >)
/// * `sequence` - Sequence data
pub fn write_fasta_record<W: Write>(
    writer: &mut W,
    header: &[u8],
    sequence: &[u8],
    preserve_mask: bool,
) -> Result<()> {
    writer
        .write_all(b">")
        .and_then(|_| writer.write_all(header))
        .and_then(|_| writer.write_all(b"\n"))
        .and_then(|_| {
            if preserve_mask {
                writer.write_all(sequence)
            } else {
                writer.write_all(sequence.to_ascii_uppercase().as_slice())
            }
        })
        .and_then(|_| writer.write_all(b"\n"))
        .map_err(|err| GenomeMaskError::io("cannot write FASTA output", err))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_writer_normalizes_lowercase_bases_by_default() {
        let mut output = Vec::new();

        write_fasta_record(&mut output, b"chr1", b"ACgtNn", false).expect("write FASTA");

        assert_eq!(output, b">chr1\nACGTNN\n");
    }

    #[test]
    fn fasta_writer_preserves_lowercase_bases_when_requested() {
        let mut output = Vec::new();

        write_fasta_record(&mut output, b"chr1", b"ACgtNn", true).expect("write FASTA");

        assert_eq!(output, b">chr1\nACgtNn\n");
    }
}

/// Reads a probe file to determine input format.
///
/// # Arguments
/// * `input` - Path to the input file
///
/// # Returns
/// * `Ok(Vec<u8>)` - Probe data
fn read_probe(input: &Path) -> Result<Vec<u8>> {
    let mut file = File::open(input)
        .map_err(|err| GenomeMaskError::io(format!("cannot open {}", input.display()), err))?;
    let mut probe = vec![0u8; 8192];
    let bytes_read = file
        .read(&mut probe)
        .map_err(|err| GenomeMaskError::io(format!("cannot read {}", input.display()), err))?;
    probe.truncate(bytes_read);
    Ok(probe)
}

/// Detects the input format from file magic bytes.
///
/// Checks for gzip magic (0x1f 0x8b), 2bit magic, or FASTA header.
fn sniff_input(bytes: &[u8]) -> Option<InputFormat> {
    if bytes.starts_with(&GZIP_MAGIC) {
        return Some(InputFormat::GzipFasta);
    }

    if bytes.len() >= TWOBIT_MAGIC.len()
        && (bytes[..TWOBIT_MAGIC.len()] == TWOBIT_MAGIC
            || bytes[..TWOBIT_MAGIC.len()] == TWOBIT_REV_MAGIC)
    {
        return Some(InputFormat::TwoBit);
    }

    if bytes
        .iter()
        .copied()
        .find(|byte| !byte.is_ascii_whitespace())
        .is_some_and(|byte| byte == b'>')
    {
        return Some(InputFormat::Fasta);
    }

    None
}

/// Trims line endings from a FASTA header.
fn trim_line_endings(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

/// Returns a human-readable description of a FASTA header.
fn header_display(header: &[u8]) -> String {
    String::from_utf8_lossy(header).into_owned()
}
