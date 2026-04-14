// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use memmap2::Mmap;
use std::{
    collections::HashMap, error::Error as StdError, fmt, fs::File, io::Write, ops::Range,
    path::Path, sync::Arc,
};
use twobit::{
    convert::{to_2bit, Nucleotides, SequenceLength, SequenceRead},
    nucleotide::Nucleotide,
};

use crate::error::{GenomeMaskError, Result};

/// Metadata for a FASTA record in a normalized 2bit conversion.
#[derive(Debug)]
struct RecordMeta {
    /// Start position of sequence data in the file
    sequence_start: usize,
    /// Length of the sequence
    sequence_len: usize,
    /// Hard-masked blocks (N bases)
    hard_blocks: Vec<Range<usize>>,
    /// Soft-masked blocks (lowercase bases)
    soft_blocks: Vec<Range<usize>>,
}

/// Metadata for a transformed sequence stored in a raw temporary sequence store.
#[derive(Debug, Clone)]
pub struct StoredSequenceRecord {
    pub header: String,
    pub sequence_start: usize,
    pub sequence_len: usize,
    pub hard_blocks: Vec<Range<usize>>,
    pub soft_blocks: Vec<Range<usize>>,
}

/// Writes 2bit data from a raw temporary sequence store and precomputed metadata.
pub fn write_two_bit_from_normalized_fasta<W: Write>(path: &Path, writer: &mut W) -> Result<()> {
    let reader = NormalizedFastaReader::open(path)?;
    to_2bit(writer, &reader)
        .map_err(|err| GenomeMaskError::InvalidTwoBit(format!("cannot write 2bit output: {err}")))
}

/// Analyzes a transformed sequence for 2bit writing.
pub fn analyze_sequence_for_twobit(header: &[u8], sequence: &[u8]) -> Result<StoredSequenceRecord> {
    validate_header(header)?;
    let header = std::str::from_utf8(header).map_err(|_| GenomeMaskError::InvalidTwoBitHeader {
        header: String::from_utf8_lossy(header).into_owned(),
        reason: "header must be valid UTF-8 ASCII".to_string(),
    })?;
    let (hard_blocks, soft_blocks) = scan_blocks(header, sequence)?;
    Ok(StoredSequenceRecord {
        header: header.to_string(),
        sequence_start: 0,
        sequence_len: sequence.len(),
        hard_blocks,
        soft_blocks,
    })
}
/// Writes 2bit data from a raw temporary sequence store and precomputed metadata.
pub fn write_two_bit_from_sequence_store<W: Write>(
    path: &Path,
    records: Vec<StoredSequenceRecord>,
    writer: &mut W,
) -> Result<()> {
    let reader = StoredSequenceReader::open(path, records)?;
    to_2bit(writer, &reader)
        .map_err(|err| GenomeMaskError::InvalidTwoBit(format!("cannot write 2bit output: {err}")))
}
/// Reader for normalized FASTA files that implements the 2bit conversion trait.
///
/// Normalized means: only headers and sequences (no blank lines),
/// with soft-masked bases in lowercase.
struct NormalizedFastaReader {
    /// Memory-mapped file data
    data: Arc<Mmap>,
    /// Lengths of each sequence record
    sequence_lengths: Vec<SequenceLength>,
    /// Metadata for each record
    records: Vec<RecordMeta>,
    /// Index from header to record position
    index: HashMap<String, usize>,
    /// Empty range collection for records without hard blocks
    empty_blocks: Vec<Range<usize>>,
}
impl NormalizedFastaReader {
    /// Opens a normalized FASTA file and indexes its records.
    ///
    /// # Arguments
    /// * `path` - Path to the FASTA file
    ///
    /// # Returns
    /// * `Ok(Self)` - Initialized reader with indexed records
    /// * `Err(GenomeMaskError)` - If file is invalid or empty
    fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|err| GenomeMaskError::io(format!("cannot open {}", path.display()), err))?;
        let data =
            Arc::new(unsafe { Mmap::map(&file) }.map_err(|err| {
                GenomeMaskError::io(format!("cannot mmap {}", path.display()), err)
            })?);
        if data.is_empty() {
            return Err(GenomeMaskError::EmptyInput);
        }
        let mut sequence_lengths = Vec::new();
        let mut records = Vec::new();
        let mut index = HashMap::new();
        let mut cursor = 0usize;
        while cursor < data.len() {
            let header_line = read_line(&data, &mut cursor);
            if header_line.is_empty() {
                continue;
            }
            if header_line[0] != b'>' {
                return Err(GenomeMaskError::InvalidFasta(
                    "temporary FASTA does not start with a header".to_string(),
                ));
            }
            let header_bytes = &header_line[1..];
            validate_header(header_bytes)?;
            let header = std::str::from_utf8(header_bytes).map_err(|_| {
                GenomeMaskError::InvalidTwoBitHeader {
                    header: String::from_utf8_lossy(header_bytes).into_owned(),
                    reason: "header must be valid UTF-8 ASCII".to_string(),
                }
            })?;
            if index.contains_key(header) {
                return Err(GenomeMaskError::InvalidTwoBitHeader {
                    header: header.to_string(),
                    reason: "duplicate header".to_string(),
                });
            }
            let sequence_start = cursor;
            let sequence_line = read_line(&data, &mut cursor);
            let (hard_blocks, soft_blocks) = scan_blocks(header, sequence_line)?;
            index.insert(header.to_string(), records.len());
            sequence_lengths.push(SequenceLength::new(&header, sequence_line.len()));
            records.push(RecordMeta {
                sequence_start,
                sequence_len: sequence_line.len(),
                hard_blocks,
                soft_blocks,
            });
        }
        if records.is_empty() {
            return Err(GenomeMaskError::EmptyInput);
        }
        Ok(Self {
            data,
            sequence_lengths,
            records,
            index,
            empty_blocks: Vec::new(),
        })
    }

    /// Returns a record by header.
    fn record(&self, header: &str) -> std::result::Result<&RecordMeta, Box<dyn StdError>> {
        self.index
            .get(header)
            .and_then(|index| self.records.get(*index))
            .ok_or_else(|| {
                Box::new(SequenceReadError(format!("missing record '{header}'")))
                    as Box<dyn StdError>
            })
    }
}

impl<'a> SequenceRead<'a> for NormalizedFastaReader {
    /// Returns the sequence lengths for each record.
    fn sequence_lengths(&'a self) -> std::result::Result<&'a [SequenceLength], Box<dyn StdError>> {
        Ok(&self.sequence_lengths)
    }

    /// Returns a nucleotides reader for a given sequence.
    fn nucleotides(
        &self,
        chr: &str,
    ) -> std::result::Result<Box<dyn Nucleotides>, Box<dyn StdError>> {
        let record = self.record(chr)?;
        Ok(Box::new(NormalizedFastaNucleotides {
            data: Arc::clone(&self.data),
            cursor: record.sequence_start,
            end: record.sequence_start + record.sequence_len,
        }))
    }

    /// Returns the soft-masked blocks for a given sequence.
    fn soft_masked_blocks(
        &'a self,
        chr: &str,
    ) -> std::result::Result<&'a [Range<usize>], Box<dyn StdError>> {
        let record = self.record(chr)?;
        Ok(&record.soft_blocks)
    }

    /// Returns the hard-masked blocks for a given sequence.
    fn hard_masked_blocks(
        &'a self,
        chr: &str,
    ) -> std::result::Result<&'a [Range<usize>], Box<dyn StdError>> {
        let record = self.record(chr)?;
        if record.hard_blocks.is_empty() {
            Ok(&self.empty_blocks)
        } else {
            Ok(&record.hard_blocks)
        }
    }
}

/// Reader for raw temporary sequence store files.
struct StoredSequenceReader {
    data: Arc<Mmap>,
    sequence_lengths: Vec<SequenceLength>,
    records: Vec<StoredSequenceRecord>,
    index: HashMap<String, usize>,
    empty_blocks: Vec<Range<usize>>,
}

impl StoredSequenceReader {
    /// Opens a temporary sequence store file and indexes its records.
    ///
    /// # Arguments
    /// * `path` - Path to the temporary sequence store file
    /// * `records` - Records to index
    ///
    /// # Returns
    /// * `Ok(Self)` - Initialized reader with indexed records
    /// * `Err(GenomeMaskError)` - If file is invalid or empty
    ///
    /// # Example
    /// ```rust,ignore
    /// let file = File::create("output.2bit")?;
    /// write_two_bit_from_normalized_fasta(Path::new("input.fa"), &mut file)?;
    /// ```
    fn open(path: &Path, records: Vec<StoredSequenceRecord>) -> Result<Self> {
        let file = File::open(path)
            .map_err(|err| GenomeMaskError::io(format!("cannot open {}", path.display()), err))?;
        let data =
            Arc::new(unsafe { Mmap::map(&file) }.map_err(|err| {
                GenomeMaskError::io(format!("cannot mmap {}", path.display()), err)
            })?);
        if records.is_empty() {
            return Err(GenomeMaskError::EmptyInput);
        }
        let mut sequence_lengths = Vec::with_capacity(records.len());
        let mut index = HashMap::with_capacity(records.len());
        for (record_index, record) in records.iter().enumerate() {
            if index.contains_key(&record.header) {
                return Err(GenomeMaskError::InvalidTwoBitHeader {
                    header: record.header.clone(),
                    reason: "duplicate header".to_string(),
                });
            }
            index.insert(record.header.clone(), record_index);
            sequence_lengths.push(SequenceLength::new(&record.header, record.sequence_len));
        }
        Ok(Self {
            data,
            sequence_lengths,
            records,
            index,
            empty_blocks: Vec::new(),
        })
    }

    /// Returns a record by header.
    fn record(
        &self,
        header: &str,
    ) -> std::result::Result<&StoredSequenceRecord, Box<dyn StdError>> {
        self.index
            .get(header)
            .and_then(|index| self.records.get(*index))
            .ok_or_else(|| {
                Box::new(SequenceReadError(format!("missing record '{header}'")))
                    as Box<dyn StdError>
            })
    }
}
impl<'a> SequenceRead<'a> for StoredSequenceReader {
    /// Returns the sequence lengths for each record.
    fn sequence_lengths(&'a self) -> std::result::Result<&'a [SequenceLength], Box<dyn StdError>> {
        Ok(&self.sequence_lengths)
    }

    /// Returns a nucleotides reader for a given sequence.
    fn nucleotides(
        &self,
        chr: &str,
    ) -> std::result::Result<Box<dyn Nucleotides>, Box<dyn StdError>> {
        let record = self.record(chr)?;
        Ok(Box::new(NormalizedFastaNucleotides {
            data: Arc::clone(&self.data),
            cursor: record.sequence_start,
            end: record.sequence_start + record.sequence_len,
        }))
    }

    /// Returns the soft-masked blocks for a given sequence.
    fn soft_masked_blocks(
        &'a self,
        chr: &str,
    ) -> std::result::Result<&'a [Range<usize>], Box<dyn StdError>> {
        let record = self.record(chr)?;
        Ok(&record.soft_blocks)
    }

    /// Returns the hard-masked blocks for a given sequence.
    fn hard_masked_blocks(
        &'a self,
        chr: &str,
    ) -> std::result::Result<&'a [Range<usize>], Box<dyn StdError>> {
        let record = self.record(chr)?;
        if record.hard_blocks.is_empty() {
            Ok(&self.empty_blocks)
        } else {
            Ok(&record.hard_blocks)
        }
    }
}

/// Nucleotides reader for a normalized FASTA sequence.
struct NormalizedFastaNucleotides {
    data: Arc<Mmap>,
    cursor: usize,
    end: usize,
}

impl Nucleotides for NormalizedFastaNucleotides {
    /// Reads a chunk of nucleotides from the sequence.
    ///
    /// # Arguments
    /// * `buf` - The buffer to write nucleotides to
    ///
    /// # Returns
    /// * `Ok(Option<usize>)` - Number of nucleotides read (or None if at end)
    /// * `Err(Box<dyn StdError>)` - Error for invalid nucleotides
    ///
    /// # Example
    /// ```rust,ignore
    /// let mut buf = [Nucleotide::A; 10];
    /// let count = nucleotides.read_chunk(&mut buf)?;
    /// ```
    fn read_chunk(
        &mut self,
        buf: &mut [Nucleotide],
    ) -> std::result::Result<Option<usize>, Box<dyn StdError>> {
        if self.cursor >= self.end {
            return Ok(None);
        }
        let count = (self.end - self.cursor).min(buf.len());
        for (offset, slot) in buf.iter_mut().take(count).enumerate() {
            *slot = as_nucleotide(self.data[self.cursor + offset])?;
        }
        self.cursor += count;
        Ok(Some(count))
    }
}

/// Scans a sequence for 2bit block boundaries.
///
/// # Arguments
/// * `header` - The FASTA header (without >)
/// * `sequence` - The sequence data (mutated in place)
///
/// # Returns
/// * `Ok((Vec<Range<usize>>, Vec<Range<usize>>))` - Hard and soft block ranges
///
/// # Example
/// ```rust,ignore
/// let (hard_blocks, soft_blocks) = scan_blocks(b"chr1", b"ACGTNNNNNN")?;
/// ```
#[allow(clippy::type_complexity)]
fn scan_blocks(header: &str, sequence: &[u8]) -> Result<(Vec<Range<usize>>, Vec<Range<usize>>)> {
    let mut hard_blocks = Vec::new();
    let mut soft_blocks = Vec::new();
    let mut hard_start = None;
    let mut soft_start = None;
    for (offset, byte) in sequence.iter().copied().enumerate() {
        match byte {
            b'A' | b'C' | b'G' | b'T' | b'N' => {}
            b'a' | b'c' | b'g' | b't' | b'n' => {}
            _ => {
                return Err(GenomeMaskError::InvalidTwoBitBase {
                    header: header.to_string(),
                    offset,
                    byte,
                })
            }
        }
        let hard = matches!(byte, b'N' | b'n');
        let soft = byte.is_ascii_lowercase();
        update_block(offset, hard, &mut hard_start, &mut hard_blocks);
        update_block(offset, soft, &mut soft_start, &mut soft_blocks);
    }
    close_block(sequence.len(), &mut hard_start, &mut hard_blocks);
    close_block(sequence.len(), &mut soft_start, &mut soft_blocks);
    Ok((hard_blocks, soft_blocks))
}
/// Updates block tracking based on current position.
///
/// # Arguments
/// * `offset` - Current position in sequence
/// * `enabled` - Whether this block type is active at this position
/// * `current_start` - Current block start position (if any)
/// * `blocks` - Vector of completed blocks
fn update_block(
    offset: usize,
    enabled: bool,
    current_start: &mut Option<usize>,
    blocks: &mut Vec<Range<usize>>,
) {
    match (*current_start, enabled) {
        (None, true) => *current_start = Some(offset),
        (Some(start), false) => {
            blocks.push(start..offset);
            *current_start = None;
        }
        _ => {}
    }
}
/// Closes any open block at the end of a sequence.
///
/// # Arguments
/// * `end` - End position of sequence
/// * `current_start` - Current block start position (if any)
/// * `blocks` - Vector of completed blocks
fn close_block(end: usize, current_start: &mut Option<usize>, blocks: &mut Vec<Range<usize>>) {
    if let Some(start) = current_start.take() {
        blocks.push(start..end);
    }
}

/// Reads a line from data starting at cursor position.
///
/// Handles both Unix (\n) and Windows (\r\n) line endings.
/// Updates cursor to position after the line.
///
/// # Arguments
/// * `data` - The input data bytes
/// * `cursor` - Mutable cursor position (updated in place)
///
/// # Returns
/// * `&[u8]` - The line bytes (excluding line ending)
fn read_line<'a>(data: &'a [u8], cursor: &mut usize) -> &'a [u8] {
    if *cursor >= data.len() {
        return &[];
    }

    let start = *cursor;
    let mut end = start;
    while end < data.len() && data[end] != b'\n' {
        end += 1;
    }

    *cursor = if end < data.len() { end + 1 } else { end };

    if end > start && data[end - 1] == b'\r' {
        &data[start..end - 1]
    } else {
        &data[start..end]
    }
}

/// Validates a FASTA header for 2bit compatibility.
///
/// Checks: not empty, ASCII-only, length <= 255 bytes.
///
/// # Arguments
/// * `header` - Raw header bytes (without >)
fn validate_header(header: &[u8]) -> Result<()> {
    let display = String::from_utf8_lossy(header).into_owned();
    if header.is_empty() {
        return Err(GenomeMaskError::InvalidTwoBitHeader {
            header: display,
            reason: "header must not be empty".to_string(),
        });
    }
    if !header.is_ascii() {
        return Err(GenomeMaskError::InvalidTwoBitHeader {
            header: display,
            reason: "header must be ASCII".to_string(),
        });
    }
    if header.len() > u8::MAX as usize {
        return Err(GenomeMaskError::InvalidTwoBitHeader {
            header: display,
            reason: "header is longer than 255 bytes".to_string(),
        });
    }
    Ok(())
}

/// Converts a byte to a 2bit Nucleotide.
///
/// Supports A, C, G, T, N (both cases). Other bytes return an error.
///
/// # Arguments
/// * `byte` - The byte to convert
///
/// # Returns
/// * `Ok(Nucleotide)` - Valid nucleotide
/// * `Err(Box<dyn StdError>` - Error for invalid bytes
fn as_nucleotide(byte: u8) -> std::result::Result<Nucleotide, Box<dyn StdError>> {
    match byte {
        b'A' | b'a' => Ok(Nucleotide::A),
        b'C' | b'c' => Ok(Nucleotide::C),
        b'G' | b'g' => Ok(Nucleotide::G),
        b'T' | b't' => Ok(Nucleotide::T),
        b'N' | b'n' => Ok(Nucleotide::N),
        _ => Err(Box::new(SequenceReadError(format!(
            "unsupported nucleotide {}",
            byte
        )))),
    }
}

#[derive(Debug)]
struct SequenceReadError(String);

impl fmt::Display for SequenceReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl StdError for SequenceReadError {}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs, process,
        time::{SystemTime, UNIX_EPOCH},
    };
    use twobit::TwoBitFile;

    #[test]
    fn normalized_fasta_writer_preserves_lowercase_n() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir =
            std::env::temp_dir().join(format!("genomemask-twobit-{}-{unique}", process::id()));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let fasta = temp_dir.join("tmp.fa");
        let twobit = temp_dir.join("tmp.2bit");
        fs::write(&fasta, b">chr1\nACNnTG\n").expect("write FASTA");

        let mut writer = File::create(&twobit).expect("create 2bit");
        write_two_bit_from_normalized_fasta(&fasta, &mut writer).expect("write 2bit");

        let mut genome = TwoBitFile::open(&twobit)
            .expect("open 2bit")
            .enable_softmask(true);
        assert_eq!(
            genome.read_sequence("chr1", ..).expect("read chr1"),
            "ACNnTG"
        );

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }
}
