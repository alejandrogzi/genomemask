// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use genepred::{Bed3, Reader};
use log::{info, warn};
use std::{collections::HashMap, path::Path};

use crate::{
    cli::ReplacementSpec,
    error::{GenomeMaskError, Result},
    io::RecordTransformer,
};

/// Represents a selenocysteine site from a BED3 file.
///
/// # Fields
/// * `chrom` - Chromosome name
/// * `start` - Start position (0-based)
/// * `end` - End position (exclusive)
/// * `line_number` - Line number in BED file for error reporting
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SelenoSite {
    pub chrom: Vec<u8>,
    pub start: usize,
    pub end: usize,
    pub line_number: usize,
}

/// Index of selenocysteine sites grouped by chromosome.
#[derive(Debug, Default)]
struct SiteIndex {
    by_chrom: HashMap<Vec<u8>, Vec<SelenoSite>>,
}

impl SiteIndex {
    /// Creates a site index from a BED3 file.
    ///
    /// # Arguments
    /// * `path` - Path to BED3 file with TGA codon coordinates
    ///
    /// # Returns
    /// * `Ok(Self)` - Indexed sites grouped by chromosome
    /// * `Err(GenomeMaskError)` - If BED3 is invalid or malformed
    fn from_bed3(path: &Path) -> Result<Self> {
        let mut reader = Reader::<Bed3>::from_mmap(path).map_err(|err| {
            GenomeMaskError::InvalidBed(format!("cannot read BED3 {}: {err}", path.display()))
        })?;

        let mut by_chrom: HashMap<Vec<u8>, Vec<SelenoSite>> = HashMap::new();

        for (index, record) in reader.records().enumerate() {
            let line_number = index + 1;
            let record = record.map_err(|err| {
                GenomeMaskError::InvalidBed(format!(
                    "cannot parse BED3 {} at logical record {}: {err}",
                    path.display(),
                    line_number
                ))
            })?;

            let start = usize::try_from(record.start()).map_err(|_| {
                GenomeMaskError::InvalidBed(format!(
                    "BED3 record {}:{}-{} is too large to fit in memory indexing",
                    String::from_utf8_lossy(record.chrom()),
                    record.start(),
                    record.end()
                ))
            })?;
            let end = usize::try_from(record.end()).map_err(|_| {
                GenomeMaskError::InvalidBed(format!(
                    "BED3 record {}:{}-{} is too large to fit in memory indexing",
                    String::from_utf8_lossy(record.chrom()),
                    record.start(),
                    record.end()
                ))
            })?;

            if end < start {
                return Err(GenomeMaskError::InvalidBed(format!(
                    "BED3 record {}:{}-{} has end < start",
                    String::from_utf8_lossy(record.chrom()),
                    record.start(),
                    record.end()
                )));
            }

            if end - start != 3 {
                warn!(
                    "skipping BED3 record {}:{}-{}: interval is not length 3",
                    String::from_utf8_lossy(record.chrom()),
                    record.start(),
                    record.end()
                );
                continue;
            }

            by_chrom
                .entry(record.chrom().to_vec())
                .or_default()
                .push(SelenoSite {
                    chrom: record.chrom().to_vec(),
                    start,
                    end,
                    line_number,
                });
        }

        for sites in by_chrom.values_mut() {
            sites.sort_by_key(|site| (site.start, site.end));

            let mut deduped = Vec::with_capacity(sites.len());
            for site in std::mem::take(sites) {
                if deduped.last().is_some_and(|last: &SelenoSite| {
                    last.start == site.start && last.end == site.end
                }) {
                    warn!(
                        "skipping duplicate BED3 entry for {}:{}-{} (line {})",
                        String::from_utf8_lossy(&site.chrom),
                        site.start,
                        site.end,
                        site.line_number
                    );
                } else {
                    deduped.push(site);
                }
            }
            *sites = deduped;
        }

        Ok(Self { by_chrom })
    }

    /// Takes and removes all sites for a given header.
    fn take_for_header(&mut self, header: &[u8]) -> Vec<SelenoSite> {
        self.by_chrom
            .remove(sequence_key(header))
            .unwrap_or_default()
    }

    /// Warns about any sites that were never matched to a genome record.
    ///
    /// Sites whose chromosome was absent from the genome (e.g. unplaced scaffolds)
    /// are reported as warnings and silently skipped rather than aborting.
    fn ensure_consumed(&self) {
        if self.by_chrom.is_empty() {
            return;
        }

        let mut chroms: Vec<&Vec<u8>> = self.by_chrom.keys().collect();
        chroms.sort();
        for chrom in chroms {
            let sites = &self.by_chrom[chrom];
            warn!(
                "chromosome '{}' not found in genome: skipping {} BED3 site(s)",
                String::from_utf8_lossy(chrom),
                sites.len()
            );
        }
    }

    #[cfg(test)]
    fn from_sites(sites: &[(&str, usize, usize)]) -> Self {
        let mut by_chrom: HashMap<Vec<u8>, Vec<SelenoSite>> = HashMap::new();

        for (index, (chrom, start, end)) in sites.iter().enumerate() {
            by_chrom
                .entry(chrom.as_bytes().to_vec())
                .or_default()
                .push(SelenoSite {
                    chrom: chrom.as_bytes().to_vec(),
                    start: *start,
                    end: *end,
                    line_number: index + 1,
                });
        }

        for sites in by_chrom.values_mut() {
            sites.sort_by_key(|site| (site.start, site.end));
        }

        Self { by_chrom }
    }
}

/// Transformer that masks selenocysteine TGA codons.
///
/// # Example
/// ```rust,ignore
/// let mut transformer = SelenoTransformer::from_bed3(
///     Path::new("seleno.bed"),
///     ReplacementSpec::Fixed(b'A')
/// )?;
/// let events = transformer.transform_record(b"chr1", &mut sequence, 0)?;
/// transformer.finish()?; // Verify all sites consumed
/// ```
pub struct SelenoTransformer {
    site_index: SiteIndex,
    replacement_spec: ReplacementSpec,
}

impl SelenoTransformer {
    /// Creates a new SelenoTransformer from a BED3 file.
    ///
    /// # Arguments
    /// * `path` - Path to BED3 file containing TGA codon coordinates
    /// * `replacement_spec` - How to replace the TGA codons
    ///
    /// # Returns
    /// * `Ok(Self)` - Initialized transformer
    /// * `Err(GenomeMaskError)` - If BED3 file is invalid
    pub fn from_bed3(path: &Path, replacement_spec: ReplacementSpec) -> Result<Self> {
        Ok(Self {
            site_index: SiteIndex::from_bed3(path)?,
            replacement_spec,
        })
    }

    #[cfg(test)]
    fn from_sites(sites: &[(&str, usize, usize)], replacement_spec: ReplacementSpec) -> Self {
        Self {
            site_index: SiteIndex::from_sites(sites),
            replacement_spec,
        }
    }
}

impl RecordTransformer for SelenoTransformer {
    fn transform_record(
        &mut self,
        header: &[u8],
        sequence: &mut Vec<u8>,
        record_index: u64,
    ) -> Result<u64> {
        let sites = self.site_index.take_for_header(header);
        mask_selenocysteines(
            sequence,
            header,
            &sites,
            self.replacement_spec,
            record_index,
        )
    }

    fn finish(&mut self) -> Result<()> {
        self.site_index.ensure_consumed();
        Ok(())
    }
}

/// Represents a resolved selenocysteine site with position adjustment.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
struct ResolvedSite {
    /// Resolved start position
    start: usize,
    /// Resolved end position (exclusive)
    end: usize,
    /// Shift applied from original position (-1, 0, or +1)
    shift: isize,
    /// Original start from BED file
    original_start: usize,
    /// Original end from BED file
    original_end: usize,
}

/// Masks selenocysteine TGA codons in a sequence.
///
/// Validates that sites don't overlap after resolution.
fn mask_selenocysteines(
    sequence: &mut [u8],
    header: &[u8],
    sites: &[SelenoSite],
    replacement_spec: ReplacementSpec,
    record_index: u64,
) -> Result<u64> {
    if sites.is_empty() {
        return Ok(0);
    }

    let mut resolved_sites = Vec::with_capacity(sites.len());
    for site in sites {
        resolved_sites.push(resolve_site(sequence, header, site)?);
    }

    resolved_sites.sort_by_key(|site| (site.start, site.end));
    for pair in resolved_sites.windows(2) {
        if pair[1].start < pair[0].end {
            return Err(GenomeMaskError::InvalidSelenocysteine(format!(
                "resolved selenocysteine sites overlap on '{}' at {}-{} and {}-{}",
                String::from_utf8_lossy(header),
                pair[0].start,
                pair[0].end,
                pair[1].start,
                pair[1].end
            )));
        }
    }

    for (site_index, site) in resolved_sites.iter().enumerate() {
        if site.shift != 0 {
            info!(
                "adjusted selenocysteine site on '{}' from {}-{} to {}-{}",
                String::from_utf8_lossy(header),
                site.original_start,
                site.original_end,
                site.start,
                site.end
            );
        }

        let codon = replacement_codon(
            replacement_spec,
            record_index,
            site.start as u64,
            site_index as u64,
        );
        let original = &sequence[site.start..site.end];
        let replacement = apply_case(codon, original);
        sequence[site.start..site.end].copy_from_slice(&replacement);
    }

    Ok(resolved_sites.len() as u64)
}

/// Resolves a BED3 site to an actual TGA, TAG, or TAA codon in the sequence.
///
/// Tries exact position first, then -1 and +1 shifts.
/// Returns error if site is ambiguous (multiple TGA matches).
fn resolve_site(sequence: &[u8], header: &[u8], site: &SelenoSite) -> Result<ResolvedSite> {
    if site.end > sequence.len() {
        return Err(GenomeMaskError::InvalidSelenocysteine(format!(
            "BED3 line {} for {}:{}-{} is out of bounds for record '{}' (length {})",
            site.line_number,
            String::from_utf8_lossy(&site.chrom),
            site.start,
            site.end,
            String::from_utf8_lossy(header),
            sequence.len()
        )));
    }

    if is_stop_codon_either_strand(&sequence[site.start..site.end]) {
        return Ok(ResolvedSite {
            start: site.start,
            end: site.end,
            shift: 0,
            original_start: site.start,
            original_end: site.end,
        });
    }

    let mut shifted_matches = Vec::new();
    for shift in [-1isize, 1isize] {
        if let Some(start) = shift_start(site.start, shift) {
            let end = start.saturating_add(3);
            if end <= sequence.len() && is_stop_codon_either_strand(&sequence[start..end]) {
                shifted_matches.push((start, end, shift));
            }
        }
    }

    match shifted_matches.as_slice() {
        [(start, end, shift)] => Ok(ResolvedSite {
            start: *start,
            end: *end,
            shift: *shift,
            original_start: site.start,
            original_end: site.end,
        }),
        [] => Err(GenomeMaskError::InvalidSelenocysteine(format!(
            "BED3 line {} for {}:{}-{} did not resolve to a stop codon on either strand for record '{}' (exact={}, -1={}, +1={})",
            site.line_number,
            String::from_utf8_lossy(&site.chrom),
            site.start,
            site.end,
            String::from_utf8_lossy(header),
            codon_display(sequence.get(site.start..site.end)),
            shifted_display(sequence, site.start, -1),
            shifted_display(sequence, site.start, 1),
        ))),
        _ => Err(GenomeMaskError::InvalidSelenocysteine(format!(
            "BED3 line {} for {}:{}-{} is ambiguous on record '{}' because both +/-1 resolve to a stop codon on either strand",
            site.line_number,
            String::from_utf8_lossy(&site.chrom),
            site.start,
            site.end,
            String::from_utf8_lossy(header),
        ))),
    }
}

/// Generates a replacement codon based on the replacement spec.
fn replacement_codon(
    replacement_spec: ReplacementSpec,
    record_index: u64,
    start: u64,
    site_index: u64,
) -> [u8; 3] {
    match replacement_spec {
        ReplacementSpec::Fixed(base) => [base, base, base],
        ReplacementSpec::Stochastic { seed } => {
            stochastic_codon(seed, record_index, start, site_index)
        }
    }
}

/// Generates a stochastic codon that is not a stop codon.
///
/// Uses SplitMix64 for deterministic randomness based on seed,
/// record index, start position, and site index.
fn stochastic_codon(seed: u64, record_index: u64, start: u64, site_index: u64) -> [u8; 3] {
    let mut remaining = (splitmix64(
        seed ^ record_index.wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ start.wrapping_mul(0xBF58_476D_1CE4_E5B9)
            ^ site_index.wrapping_mul(0x94D0_49BB_1331_11EB),
    ) % 61) as usize;

    const BASES: [u8; 4] = [b'A', b'T', b'C', b'G'];

    for first in BASES {
        for second in BASES {
            for third in BASES {
                let codon = [first, second, third];
                if is_stop_codon(&codon) {
                    continue;
                }
                if remaining == 0 {
                    return codon;
                }
                remaining -= 1;
            }
        }
    }

    unreachable!("61 non-stop codons must always exist")
}

/// Applies the case style from original bases to replacement codon.
///
/// If original bases are lowercase, the replacement is also lowercase.
fn apply_case(mut replacement: [u8; 3], original: &[u8]) -> [u8; 3] {
    for (slot, original_base) in replacement.iter_mut().zip(original.iter().copied()) {
        if original_base.is_ascii_lowercase() {
            *slot = slot.to_ascii_lowercase();
        }
    }
    replacement
}

/// Calculates a shifted start position with overflow checking.
///
/// # Arguments
/// * `start` - Original start position
/// * `shift` - Shift amount (-1 or +1)
///
/// # Returns
/// * `Some(usize)` - Shifted position if valid
/// * `None` - If shift would underflow/overflow
fn shift_start(start: usize, shift: isize) -> Option<usize> {
    if shift.is_negative() {
        start.checked_sub(shift.unsigned_abs())
    } else {
        start.checked_add(shift as usize)
    }
}

/// Formats a shifted codon for error messages.
fn shifted_display(sequence: &[u8], start: usize, shift: isize) -> String {
    shift_start(start, shift)
        .and_then(|candidate_start| {
            let candidate_end = candidate_start.checked_add(3)?;
            sequence.get(candidate_start..candidate_end)
        })
        .map_or_else(
            || "out-of-bounds".to_string(),
            |codon| codon_display(Some(codon)),
        )
}

/// Formats a codon for error messages.
fn codon_display(codon: Option<&[u8]>) -> String {
    codon
        .map(|bytes| String::from_utf8_lossy(bytes).into_owned())
        .unwrap_or_else(|| "out-of-bounds".to_string())
}

/// Checks if a codon is a stop codon (TAA, TAG, or TGA), case-insensitively.
fn is_stop_codon(codon: &[u8]) -> bool {
    codon.len() == 3
        && matches!(
            [
                codon[0].to_ascii_uppercase(),
                codon[1].to_ascii_uppercase(),
                codon[2].to_ascii_uppercase(),
            ],
            [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
        )
}

/// Returns the uppercase DNA complement of a single base.
fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        other => other,
    }
}

/// Returns the reverse complement of a 3-base codon (always uppercase).
fn rev_comp_codon(codon: [u8; 3]) -> [u8; 3] {
    [
        complement(codon[2]),
        complement(codon[1]),
        complement(codon[0]),
    ]
}

/// Checks if a codon is a stop codon on either the forward or reverse-complement strand.
fn is_stop_codon_either_strand(codon: &[u8]) -> bool {
    if codon.len() != 3 {
        return false;
    }
    is_stop_codon(codon) || is_stop_codon(&rev_comp_codon([codon[0], codon[1], codon[2]]))
}

/// SplitMix64 PRNG for deterministic stochastic selection.
///
/// A fast, good-quality PRNG used for reproducible random base selection.
fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

/// Extracts the first whitespace-separated token from a FASTA header.
///
/// This is used as the key for matching genome records to selenocysteine sites.
fn sequence_key(header: &[u8]) -> &[u8] {
    header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or(header)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequence_key_uses_first_token() {
        assert_eq!(sequence_key(b"chr1 some description"), b"chr1");
        assert_eq!(sequence_key(b"chr2"), b"chr2");
    }

    #[test]
    fn fixed_replacement_masks_tga_triplet() {
        let mut sequence = b"CCCTGAaaa".to_vec();
        let mut transformer =
            SelenoTransformer::from_sites(&[("chr1", 3, 6)], ReplacementSpec::Fixed(b'G'));

        let count = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mask sites");

        assert_eq!(count, 1);
        assert_eq!(sequence, b"CCCGGGaaa");
    }

    #[test]
    fn shifted_coordinate_is_rescued() {
        let mut sequence = b"ACTTGACCC".to_vec();
        let mut transformer =
            SelenoTransformer::from_sites(&[("chr1", 4, 7)], ReplacementSpec::Fixed(b'A'));

        let count = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mask sites");

        assert_eq!(count, 1);
        assert_eq!(sequence, b"ACTAAACCC");
    }

    #[test]
    fn stochastic_replacement_never_emits_stop_codon() {
        let mut sequence = b"TGA".to_vec();
        let mut transformer = SelenoTransformer::from_sites(
            &[("chr1", 0, 3)],
            ReplacementSpec::Stochastic { seed: 7 },
        );

        transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mask sites");

        assert!(!matches!(sequence.as_slice(), b"TAA" | b"TAG" | b"TGA"));
    }

    // --- non-length-3: warn and skip ---

    #[test]
    fn non_length_3_record_is_skipped_not_an_error() {
        use std::io::Write as _;
        let path = std::env::temp_dir().join("genomemask_test_non3.bed");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "chr1\t100\t104").unwrap(); // length 4 — should be skipped
            writeln!(f, "chr1\t200\t203").unwrap(); // length 3 — kept
        }
        let index = SiteIndex::from_bed3(&path).expect("non-length-3 should not be an error");
        let sites = index.by_chrom.get(b"chr1".as_slice()).map(|v| v.len());
        assert_eq!(sites, Some(1), "only the length-3 record should be kept");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn non_length_3_only_record_yields_empty_index() {
        use std::io::Write as _;
        let path = std::env::temp_dir().join("genomemask_test_non3_only.bed");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "chr1\t100\t106").unwrap(); // length 6 — skipped
        }
        let index = SiteIndex::from_bed3(&path).expect("non-length-3 should not be an error");
        assert!(index.by_chrom.is_empty(), "no valid records should remain");
        let _ = std::fs::remove_file(&path);
    }

    // --- duplicates: warn and deduplicate ---

    #[test]
    fn duplicate_entries_are_deduplicated_not_an_error() {
        use std::io::Write as _;
        let path = std::env::temp_dir().join("genomemask_test_dup.bed");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "chr1\t100\t103").unwrap();
            writeln!(f, "chr1\t100\t103").unwrap(); // duplicate
            writeln!(f, "chr1\t200\t203").unwrap();
        }
        let index = SiteIndex::from_bed3(&path).expect("duplicates should not be an error");
        let count = index.by_chrom.get(b"chr1".as_slice()).map(|v| v.len());
        assert_eq!(
            count,
            Some(2),
            "duplicate should be dropped, leaving 2 unique sites"
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn three_copies_reduce_to_one() {
        use std::io::Write as _;
        let path = std::env::temp_dir().join("genomemask_test_dup3.bed");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            for _ in 0..3 {
                writeln!(f, "chr1\t50\t53").unwrap();
            }
        }
        let index = SiteIndex::from_bed3(&path).expect("triplicate should not be an error");
        let count = index.by_chrom.get(b"chr1".as_slice()).map(|v| v.len());
        assert_eq!(count, Some(1), "only one copy should survive deduplication");
        let _ = std::fs::remove_file(&path);
    }

    // --- two-strand stop codon validation ---

    #[test]
    fn is_stop_codon_either_strand_forward() {
        assert!(is_stop_codon_either_strand(b"TGA"));
        assert!(is_stop_codon_either_strand(b"TAA"));
        assert!(is_stop_codon_either_strand(b"TAG"));
        // soft-masked (lowercase) — is_stop_codon normalizes to uppercase
        assert!(is_stop_codon_either_strand(b"tga"));
        assert!(is_stop_codon_either_strand(b"taa"));
        assert!(is_stop_codon_either_strand(b"tag"));
        // mixed case (e.g. TgA)
        assert!(is_stop_codon_either_strand(b"TgA"));
        assert!(is_stop_codon_either_strand(b"tAg"));
    }

    #[test]
    fn unresolved_chromosome_does_not_error() {
        use std::io::Write as _;
        let path = std::env::temp_dir().join("genomemask_test_unresolved.bed");
        {
            let mut f = std::fs::File::create(&path).unwrap();
            writeln!(f, "chrUn_scaffold42\t0\t3").unwrap();
        }
        // Loading succeeds
        let mut transformer =
            SelenoTransformer::from_bed3(&path, ReplacementSpec::Fixed(b'N')).unwrap();
        // Processing a different chromosome produces no events
        let mut seq = b"ATGTGA".to_vec();
        let count = transformer
            .transform_record(b"chr1", &mut seq, 0)
            .expect("unrelated record should process cleanly");
        assert_eq!(count, 0);
        // finish() should not error even though chrUn_scaffold42 was never seen
        transformer
            .finish()
            .expect("unresolved chromosome should only warn, not error");
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn is_stop_codon_either_strand_reverse_complement() {
        // rc(TGA) = TCA, rc(TAA) = TTA, rc(TAG) = CTA
        assert!(is_stop_codon_either_strand(b"TCA"));
        assert!(is_stop_codon_either_strand(b"TTA"));
        assert!(is_stop_codon_either_strand(b"CTA"));
    }

    #[test]
    fn is_stop_codon_either_strand_non_stop() {
        assert!(!is_stop_codon_either_strand(b"ATG"));
        assert!(!is_stop_codon_either_strand(b"GGG"));
        assert!(!is_stop_codon_either_strand(b"CCC"));
    }

    #[test]
    fn reverse_complement_stop_codon_is_masked() {
        // TCA is rc(TGA), so a BED record pointing to TCA should be accepted and masked
        let mut sequence = b"CCCTCA".to_vec();
        let mut transformer =
            SelenoTransformer::from_sites(&[("chr1", 3, 6)], ReplacementSpec::Fixed(b'N'));

        let count = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("reverse-complement stop codon should be accepted");

        assert_eq!(count, 1);
        assert_eq!(sequence, b"CCCNNN");
    }

    #[test]
    fn forward_and_reverse_stop_codons_coexist() {
        // TAA at 0..3 (forward), TCA at 3..6 (rc of TGA)
        let mut sequence = b"TAATCA".to_vec();
        let mut transformer = SelenoTransformer::from_sites(
            &[("chr1", 0, 3), ("chr1", 3, 6)],
            ReplacementSpec::Fixed(b'G'),
        );

        let count = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("mixed forward/reverse stop codons should be accepted");

        assert_eq!(count, 2);
        assert_eq!(sequence, b"GGGGGG");
    }
}
