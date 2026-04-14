// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use flate2::read::MultiGzDecoder;
use genepred::{
    bed::BedFormat, Bed12, Bed3, Bed4, Bed5, Bed6, Bed8, Bed9, GenePred, Gff, Gtf, Reader,
    ReaderOptions,
};
use rayon::prelude::*;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::{
    cli::{MaskFeature, MaskSelection, ReplacementSpec},
    commands::ns::replacement_base_at,
    error::{GenomeMaskError, Result},
    io::RecordTransformer,
};

/// Minimum sequence length to enable parallel processing.
const PARALLEL_MIN_LEN: usize = 1 << 20;
/// Chunk size for parallel processing.
const PARALLEL_CHUNK_LEN: usize = 1 << 20;

/// Represents a genomic interval from a regions file.
///
/// # Fields
/// * `chrom` - Chromosome name
/// * `start` - Start position (0-based)
/// * `end` - End position (exclusive)
/// * `line_number` - Line number in file for error reporting
#[derive(Debug, Clone, Eq, PartialEq)]
struct Interval {
    chrom: Vec<u8>,
    start: usize,
    end: usize,
    line_number: usize,
}

/// Index of intervals grouped by chromosome.
#[derive(Debug, Default)]
struct IntervalIndex {
    by_chrom: HashMap<Vec<u8>, Vec<Interval>>,
}

/// Detected regions file format.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum RegionsFormat {
    /// BED format with variant and extra fields count
    Bed(BedKind, usize),
    /// GTF format
    Gtf,
    /// GFF format
    Gff,
}

/// BED format variant based on field count.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum BedKind {
    /// 3-field BED (chrom, start, end)
    Bed3,
    /// 4-field BED (+ name)
    Bed4,
    /// 5-field BED (+ score)
    Bed5,
    /// 6-field BED (+ strand)
    Bed6,
    /// 8-field BED (thick start/end)
    Bed8,
    /// 9-field BED (item RGB)
    Bed9,
    /// 12-field BED (blocks)
    Bed12,
}

/// Transformer that masks genomic intervals.
///
/// # Example
/// ```rust,ignore
/// let mut transformer = MaskTransformer::from_regions(
///     Path::new("exons.bed"),
///     MaskSelection::Feature(MaskFeature::Exon),
///     ReplacementSpec::Fixed(b'N')
/// )?;
/// let events = transformer.transform_record(b"chr1", &mut sequence, 0)?;
/// transformer.finish()?;
/// ```
pub struct MaskTransformer {
    interval_index: IntervalIndex,
    replacement_spec: ReplacementSpec,
}

impl MaskTransformer {
    /// Creates a new MaskTransformer from a regions file.
    ///
    /// # Arguments
    /// * `path` - Path to regions file (BED, GTF, GFF)
    /// * `selection` - Selection mode (specific intervals or feature-derived)
    /// * `replacement_spec` - How to replace bases in intervals
    ///
    /// # Returns
    /// * `Ok(Self)` - Initialized transformer
    /// * `Err(GenomeMaskError)` - If regions file is invalid
    pub fn from_regions(
        path: &Path,
        selection: MaskSelection,
        replacement_spec: ReplacementSpec,
    ) -> Result<Self> {
        Ok(Self {
            interval_index: IntervalIndex::from_regions(path, selection)?,
            replacement_spec,
        })
    }

    #[cfg(test)]
    fn from_intervals(
        intervals: &[(&str, usize, usize)],
        replacement_spec: ReplacementSpec,
    ) -> Self {
        Self {
            interval_index: IntervalIndex::from_intervals(intervals),
            replacement_spec,
        }
    }
}

impl RecordTransformer for MaskTransformer {
    fn transform_record(
        &mut self,
        header: &[u8],
        sequence: &mut Vec<u8>,
        record_index: u64,
    ) -> Result<u64> {
        let intervals = self.interval_index.take_for_header(header);
        mask_intervals(
            sequence,
            header,
            &intervals,
            self.replacement_spec,
            record_index,
        )
    }

    fn finish(&mut self) -> Result<()> {
        self.interval_index.ensure_consumed()
    }
}

impl IntervalIndex {
    /// Loads intervals from a regions file based on selection mode.
    ///
    /// # Arguments
    /// * `path` - Path to regions file
    /// * `selection` - How to interpret the regions
    ///
    /// # Returns
    /// * `Ok(Self)` - Indexed intervals
    /// * `Err(GenomeMaskError)` - If format doesn't match selection
    fn from_regions(path: &Path, selection: MaskSelection) -> Result<Self> {
        let format = detect_regions_format(path)?;
        let intervals = match (selection, format) {
            (MaskSelection::Specific, RegionsFormat::Bed(kind, extra_fields)) => {
                load_specific_intervals(path, kind, extra_fields)?
            }
            (MaskSelection::Specific, RegionsFormat::Gtf | RegionsFormat::Gff) => {
                return Err(GenomeMaskError::InvalidArgument(
                    "--specific is only supported for BED inputs".to_string(),
                ))
            }
            (MaskSelection::Feature(_), RegionsFormat::Bed(kind, _)) if kind != BedKind::Bed12 => {
                return Err(GenomeMaskError::InvalidArgument(
                    "transcript-derived feature masking from BED requires BED12; use --specific for direct BED intervals"
                        .to_string(),
                ))
            }
            (MaskSelection::Feature(feature), RegionsFormat::Bed(BedKind::Bed12, extra_fields)) => {
                load_feature_intervals::<Bed12>(path, feature, extra_fields)?
            }
            (MaskSelection::Feature(feature), RegionsFormat::Gtf) => {
                load_feature_intervals::<Gtf>(path, feature, 0)?
            }
            (MaskSelection::Feature(feature), RegionsFormat::Gff) => {
                load_feature_intervals::<Gff>(path, feature, 0)?
            }
            (MaskSelection::Feature(_), RegionsFormat::Bed(_, _)) => unreachable!("handled above"),
        };

        if intervals.is_empty() {
            let label = match selection {
                MaskSelection::Specific => "direct BED intervals".to_string(),
                MaskSelection::Feature(feature) => {
                    format!("{} intervals", mask_feature_name(feature))
                }
            };
            return Err(GenomeMaskError::InvalidRegions(format!(
                "no {label} were derived from {}",
                path.display()
            )));
        }

        Ok(Self::from_loaded_intervals(intervals))
    }

    /// Creates an index from pre-loaded intervals.
    ///
    /// Sorts intervals by position and merges overlapping intervals.
    fn from_loaded_intervals(intervals: Vec<Interval>) -> Self {
        let mut by_chrom: HashMap<Vec<u8>, Vec<Interval>> = HashMap::new();

        for interval in intervals {
            by_chrom
                .entry(interval.chrom.clone())
                .or_default()
                .push(interval);
        }

        for intervals in by_chrom.values_mut() {
            intervals.sort_by_key(|interval| (interval.start, interval.end));
            *intervals = merge_intervals(intervals);
        }

        Self { by_chrom }
    }

    /// Takes and removes all intervals for a given header.
    fn take_for_header(&mut self, header: &[u8]) -> Vec<Interval> {
        self.by_chrom
            .remove(sequence_key(header))
            .unwrap_or_default()
    }

    /// Ensures all intervals were consumed (matched to genome records).
    fn ensure_consumed(&self) -> Result<()> {
        if self.by_chrom.is_empty() {
            return Ok(());
        }

        let mut unresolved: Vec<String> = self
            .by_chrom
            .iter()
            .map(|(chrom, intervals)| {
                let first = &intervals[0];
                format!(
                    "{}:{}-{}",
                    String::from_utf8_lossy(chrom),
                    first.start,
                    first.end
                )
            })
            .collect();
        unresolved.sort();

        Err(GenomeMaskError::InvalidRegions(format!(
            "the following mask intervals were not found in the genome headers: {}",
            unresolved.join(", ")
        )))
    }

    #[cfg(test)]
    fn from_intervals(intervals: &[(&str, usize, usize)]) -> Self {
        let loaded = intervals
            .iter()
            .enumerate()
            .map(|(index, (chrom, start, end))| Interval {
                chrom: chrom.as_bytes().to_vec(),
                start: *start,
                end: *end,
                line_number: index + 1,
            })
            .collect();
        Self::from_loaded_intervals(loaded)
    }
}

/// Loads specific (direct) intervals from a BED file.
///
/// Dispatches to the appropriate reader based on BED variant.
fn load_specific_intervals(
    path: &Path,
    kind: BedKind,
    extra_fields: usize,
) -> Result<Vec<Interval>> {
    match kind {
        BedKind::Bed3 => load_specific_intervals_from_reader::<Bed3>(path, extra_fields),
        BedKind::Bed4 => load_specific_intervals_from_reader::<Bed4>(path, extra_fields),
        BedKind::Bed5 => load_specific_intervals_from_reader::<Bed5>(path, extra_fields),
        BedKind::Bed6 => load_specific_intervals_from_reader::<Bed6>(path, extra_fields),
        BedKind::Bed8 => load_specific_intervals_from_reader::<Bed8>(path, extra_fields),
        BedKind::Bed9 => load_specific_intervals_from_reader::<Bed9>(path, extra_fields),
        BedKind::Bed12 => load_specific_intervals_from_reader::<Bed12>(path, extra_fields),
    }
}

/// Loads specific intervals using a generic BED reader.
fn load_specific_intervals_from_reader<R>(path: &Path, extra_fields: usize) -> Result<Vec<Interval>>
where
    R: BedFormat + Into<GenePred>,
{
    let mut reader = open_regions_reader::<R>(path, extra_fields)?;

    let mut intervals = Vec::new();
    for (index, record) in reader.records().enumerate() {
        let record = record.map_err(|err| {
            GenomeMaskError::InvalidRegions(format!(
                "cannot parse regions file {} at logical record {}: {err}",
                path.display(),
                index + 1
            ))
        })?;
        intervals.push(interval_from_bounds(
            record.chrom(),
            record.start(),
            record.end(),
            index + 1,
        )?);
    }

    Ok(intervals)
}

/// Loads feature-derived intervals from a regions file.
///
/// Extracts CDS, exon, intron, or UTR intervals from transcript records.
fn load_feature_intervals<R>(
    path: &Path,
    feature: MaskFeature,
    extra_fields: usize,
) -> Result<Vec<Interval>>
where
    R: BedFormat + Into<GenePred>,
{
    let mut reader = open_regions_reader::<R>(path, extra_fields)?;

    let mut intervals = Vec::new();
    for (index, record) in reader.records().enumerate() {
        let record = record.map_err(|err| {
            GenomeMaskError::InvalidRegions(format!(
                "cannot parse regions file {} at logical record {}: {err}",
                path.display(),
                index + 1
            ))
        })?;

        for (start, end) in feature_intervals(&record, feature) {
            intervals.push(interval_from_bounds(record.chrom(), start, end, index + 1)?);
        }
    }

    Ok(intervals)
}

/// Extracts intervals for a specific feature from a GenePred record.
fn feature_intervals(record: &GenePred, feature: MaskFeature) -> Vec<(u64, u64)> {
    match feature {
        MaskFeature::Cds => record.coding_exons(),
        MaskFeature::Exon => record.exons(),
        MaskFeature::Intron => record.introns(),
        MaskFeature::Utr => record.utr_exons(),
    }
}

/// Creates an Interval from genomic coordinates.
///
/// Validates that start < end and coordinates fit in usize.
fn interval_from_bounds(
    chrom: &[u8],
    start: u64,
    end: u64,
    line_number: usize,
) -> Result<Interval> {
    let start = usize::try_from(start).map_err(|_| {
        GenomeMaskError::InvalidRegions(format!(
            "interval {}:{}-{} is too large to fit in memory indexing",
            String::from_utf8_lossy(chrom),
            start,
            end
        ))
    })?;
    let end = usize::try_from(end).map_err(|_| {
        GenomeMaskError::InvalidRegions(format!(
            "interval {}:{}-{} is too large to fit in memory indexing",
            String::from_utf8_lossy(chrom),
            start,
            end
        ))
    })?;

    if end < start {
        return Err(GenomeMaskError::InvalidRegions(format!(
            "interval {}:{}-{} has end < start",
            String::from_utf8_lossy(chrom),
            start,
            end
        )));
    }

    if start == end {
        return Err(GenomeMaskError::InvalidRegions(format!(
            "interval {}:{}-{} is empty",
            String::from_utf8_lossy(chrom),
            start,
            end
        )));
    }

    Ok(Interval {
        chrom: chrom.to_vec(),
        start,
        end,
        line_number,
    })
}

/// Merges overlapping and adjacent intervals.
///
/// # Arguments
/// * `intervals` - Sorted list of intervals
///
/// # Returns
/// * `Vec<Interval>` - Merged non-overlapping intervals
fn merge_intervals(intervals: &[Interval]) -> Vec<Interval> {
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut merged = Vec::with_capacity(intervals.len());
    let mut current = intervals[0].clone();

    for interval in &intervals[1..] {
        if interval.start <= current.end {
            current.end = current.end.max(interval.end);
        } else {
            merged.push(current);
            current = interval.clone();
        }
    }

    merged.push(current);
    merged
}

/// Masks all intervals in a sequence.
///
/// Validates intervals are within sequence bounds.
fn mask_intervals(
    sequence: &mut [u8],
    header: &[u8],
    intervals: &[Interval],
    replacement_spec: ReplacementSpec,
    record_index: u64,
) -> Result<u64> {
    let mut masked_bases = 0u64;

    for interval in intervals {
        if interval.end > sequence.len() {
            return Err(GenomeMaskError::InvalidRegions(format!(
                "regions line {} for {}:{}-{} is out of bounds for record '{}' (length {})",
                interval.line_number,
                String::from_utf8_lossy(&interval.chrom),
                interval.start,
                interval.end,
                String::from_utf8_lossy(header),
                sequence.len()
            )));
        }

        masked_bases += mask_interval_range(
            &mut sequence[interval.start..interval.end],
            replacement_spec,
            record_index,
            interval.start as u64,
        );
    }

    Ok(masked_bases)
}

/// Masks a range of sequence with parallel processing support.
///
/// Uses parallel chunks for sequences >= 1MB when multiple threads available.
fn mask_interval_range(
    sequence: &mut [u8],
    replacement_spec: ReplacementSpec,
    record_index: u64,
    start_offset: u64,
) -> u64 {
    if sequence.len() >= PARALLEL_MIN_LEN && rayon::current_num_threads() > 1 {
        sequence
            .par_chunks_mut(PARALLEL_CHUNK_LEN)
            .enumerate()
            .map(|(chunk_index, chunk)| {
                let chunk_offset = start_offset + (chunk_index * PARALLEL_CHUNK_LEN) as u64;
                mask_interval_chunk(chunk, replacement_spec, record_index, chunk_offset)
            })
            .sum()
    } else {
        mask_interval_chunk(sequence, replacement_spec, record_index, start_offset)
    }
}

/// Masks a chunk of sequence data.
fn mask_interval_chunk(
    sequence: &mut [u8],
    replacement_spec: ReplacementSpec,
    record_index: u64,
    start_offset: u64,
) -> u64 {
    for (offset, base) in sequence.iter_mut().enumerate() {
        let lowercase = base.is_ascii_lowercase();
        *base = replacement_base_at(
            replacement_spec,
            record_index,
            start_offset + offset as u64,
            lowercase,
        );
    }

    sequence.len() as u64
}

/// Opens a regions file reader with optional additional fields.
///
/// # Arguments
/// * `path` - Path to the regions file
/// * `extra_fields` - Number of additional fields beyond standard BED format
fn open_regions_reader<R>(path: &Path, extra_fields: usize) -> Result<Reader<R>>
where
    R: BedFormat + Into<GenePred>,
{
    if extra_fields == 0 {
        Reader::<R>::from_path(path).map_err(|err| {
            GenomeMaskError::InvalidRegions(format!(
                "cannot read regions file {}: {err}",
                path.display()
            ))
        })
    } else {
        let options = ReaderOptions::new().additional_fields(extra_fields);
        Reader::<R>::from_path_with_custom_fields(path, options).map_err(|err| {
            GenomeMaskError::InvalidRegions(format!(
                "cannot read regions file {}: {err}",
                path.display()
            ))
        })
    }
}

/// Detects the regions file format from file extension.
///
/// Supports: .bed, .bed.gz, .gtf, .gtf.gz, .gff, .gff.gz, .gff3, .gff3.gz
fn detect_regions_format(path: &Path) -> Result<RegionsFormat> {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| {
            GenomeMaskError::InvalidArgument(format!(
                "cannot determine regions format from '{}'",
                path.display()
            ))
        })?
        .to_ascii_lowercase();

    if name.ends_with(".gtf") || name.ends_with(".gtf.gz") {
        return Ok(RegionsFormat::Gtf);
    }

    if name.ends_with(".gff")
        || name.ends_with(".gff.gz")
        || name.ends_with(".gff3")
        || name.ends_with(".gff3.gz")
    {
        return Ok(RegionsFormat::Gff);
    }

    if name.ends_with(".bed") || name.ends_with(".bed.gz") {
        let (kind, extra_fields) = detect_bed_kind(path)?;
        return Ok(RegionsFormat::Bed(kind, extra_fields));
    }

    Err(GenomeMaskError::InvalidArgument(format!(
        "unsupported regions format for {}; expected BED, GTF, GFF, GTF.GZ, or GFF.GZ",
        path.display()
    )))
}

/// Detects BED format variant by counting fields in first non-empty line.
///
/// Skips comments (starting with #, track, browser).
/// Returns the BedKind and number of extra fields beyond standard.
fn detect_bed_kind(path: &Path) -> Result<(BedKind, usize)> {
    let mut reader = open_text_reader(path)?;
    let mut line = String::new();

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line).map_err(|err| {
            GenomeMaskError::io(format!("cannot inspect {}", path.display()), err)
        })?;
        if bytes_read == 0 {
            return Err(GenomeMaskError::InvalidRegions(format!(
                "regions file {} is empty",
                path.display()
            )));
        }

        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track")
            || trimmed.starts_with("browser")
        {
            continue;
        }

        let fields = trimmed
            .split('\t')
            .count()
            .max(trimmed.split_ascii_whitespace().count());
        return match fields {
            3 => Ok((BedKind::Bed3, 0)),
            4 => Ok((BedKind::Bed4, 0)),
            5 => Ok((BedKind::Bed5, 0)),
            6 => Ok((BedKind::Bed6, 0)),
            7 => Ok((BedKind::Bed6, 1)),
            8 => Ok((BedKind::Bed8, 0)),
            9 => Ok((BedKind::Bed9, 0)),
            10 => Ok((BedKind::Bed9, 1)),
            11 => Ok((BedKind::Bed9, 2)),
            n if n >= 12 => Ok((BedKind::Bed12, n - 12)),
            _ => Err(GenomeMaskError::InvalidRegions(format!(
                "unsupported BED field count {fields} in {}",
                path.display()
            ))),
        };
    }
}

/// Opens a text file reader, handling gzip compression automatically.
///
/// # Arguments
/// * `path` - Path to the text file
///
/// # Returns
/// * `Ok(Box<dyn BufRead>)` - Buffered reader for the file
fn open_text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .map_err(|err| GenomeMaskError::io(format!("cannot open {}", path.display()), err))?;
    let gzipped = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"));

    if gzipped {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Extracts the first whitespace-separated token from a FASTA header.
///
/// This is used as the key for matching genome records to intervals.
fn sequence_key(header: &[u8]) -> &[u8] {
    header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or(header)
}

/// Returns the name of a mask feature for error messages.
fn mask_feature_name(feature: MaskFeature) -> &'static str {
    match feature {
        MaskFeature::Cds => "cds",
        MaskFeature::Exon => "exon",
        MaskFeature::Intron => "intron",
        MaskFeature::Utr => "utr",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_interval_masking_preserves_case() {
        let mut sequence = b"ACGTacgt".to_vec();
        let mut transformer =
            MaskTransformer::from_intervals(&[("chr1", 2, 6)], ReplacementSpec::Fixed(b'G'));

        let masked = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mask intervals");

        assert_eq!(masked, 4);
        assert_eq!(sequence, b"ACGGgggt");
    }

    #[test]
    fn overlapping_intervals_are_merged() {
        let mut sequence = b"NNNNNN".to_vec();
        let mut transformer = MaskTransformer::from_intervals(
            &[("chr1", 1, 4), ("chr1", 3, 6)],
            ReplacementSpec::Fixed(b'A'),
        );

        let masked = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mask intervals");

        assert_eq!(masked, 5);
        assert_eq!(sequence, b"NAAAAA");
    }
}
