// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use rayon::prelude::*;

use crate::{cli::ReplacementSpec, io::RecordTransformer, Result};

/// Minimum sequence length to enable parallel processing.
const PARALLEL_MIN_LEN: usize = 1 << 20;
/// Chunk size for parallel processing.
const PARALLEL_CHUNK_LEN: usize = 1 << 20;

/// Transformer that replaces N bases in genome sequences.
///
/// # Example
/// ```rust,ignore
/// let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));
/// let events = transformer.transform_record(b"chr1", &mut sequence, 0)?;
/// ```
pub struct NsTransformer {
    /// Specification for nucleotide replacement
    replacement_spec: ReplacementSpec,
}

impl NsTransformer {
    /// Creates a new NTransformer.
    ///
    /// # Arguments
    /// * `replacement_spec` - How to replace N bases
    pub fn new(replacement_spec: ReplacementSpec) -> Self {
        Self { replacement_spec }
    }
}

impl RecordTransformer for NsTransformer {
    fn transform_record(
        &mut self,
        _header: &[u8],
        sequence: &mut Vec<u8>,
        record_index: u64,
    ) -> Result<u64> {
        Ok(transform_range_in_place(
            sequence,
            self.replacement_spec,
            record_index,
            0,
        ))
    }
}

/// Transforms a sequence by replacing N bases in-place.
///
/// Uses parallel processing for sequences >= 1MB when multiple threads available.
pub(crate) fn transform_range_in_place(
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
                let offset = start_offset + (chunk_index * PARALLEL_CHUNK_LEN) as u64;
                transform_chunk(chunk, replacement_spec, record_index, offset)
            })
            .sum()
    } else {
        transform_chunk(sequence, replacement_spec, record_index, start_offset)
    }
}

/// Transforms a chunk of sequence data by replacing N bases.
fn transform_chunk(
    sequence: &mut [u8],
    replacement_spec: ReplacementSpec,
    record_index: u64,
    start_offset: u64,
) -> u64 {
    let mut replacements = 0u64;

    for (offset, base) in sequence.iter_mut().enumerate() {
        let absolute_offset = start_offset + offset as u64;
        let replacement = match *base {
            b'N' => Some(replacement_base_at(
                replacement_spec,
                record_index,
                absolute_offset,
                false,
            )),
            b'n' => Some(replacement_base_at(
                replacement_spec,
                record_index,
                absolute_offset,
                true,
            )),
            _ => None,
        };

        if let Some(replacement) = replacement {
            *base = replacement;
            replacements += 1;
        }
    }

    replacements
}

/// Calculates the replacement base for a given position.
///
/// # Arguments
/// * `replacement_spec` - Fixed or stochastic replacement
/// * `record_index` - Current record index for seeding
/// * `offset` - Position within the sequence
/// * `lowercase` - Whether original base was lowercase (preserve case)
///
/// # Returns
/// * `u8` - Replacement nucleotide
pub(crate) fn replacement_base_at(
    replacement_spec: ReplacementSpec,
    record_index: u64,
    offset: u64,
    lowercase: bool,
) -> u8 {
    let base = match replacement_spec {
        ReplacementSpec::Fixed(base) => base,
        ReplacementSpec::Stochastic { seed } => stochastic_base(seed, record_index, offset),
    };

    if lowercase {
        base.to_ascii_lowercase()
    } else {
        base
    }
}

/// Generates a deterministic random nucleotide using SplitMix64.
///
/// # Arguments
/// * `seed` - User-provided seed
/// * `record_index` - Record index for additional mixing
/// * `offset` - Position within sequence
///
/// # Returns
/// * `u8` - Random nucleotide (A, T, C, or G)
fn stochastic_base(seed: u64, record_index: u64, offset: u64) -> u8 {
    let mixed = splitmix64(
        seed ^ record_index.wrapping_mul(0x9E37_79B9_7F4A_7C15)
            ^ offset.wrapping_mul(0xBF58_476D_1CE4_E5B9),
    );

    match mixed & 0b11 {
        0 => b'A',
        1 => b'T',
        2 => b'C',
        _ => b'G',
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixed_replacement_preserves_case() {
        let mut sequence = b"ACGTNnRY".to_vec();
        let mut transformer = NsTransformer::new(ReplacementSpec::Fixed(b'G'));

        let replacements = transformer
            .transform_record(b"chr1", &mut sequence, 0)
            .expect("transform record");

        assert_eq!(replacements, 2);
        assert_eq!(sequence, b"ACGTGgRY");
    }

    #[test]
    fn stochastic_replacement_is_deterministic() {
        let mut first = b"NNNNnnnn".to_vec();
        let mut second = b"NNNNnnnn".to_vec();
        let mut first_transformer = NsTransformer::new(ReplacementSpec::Stochastic { seed: 7 });
        let mut second_transformer = NsTransformer::new(ReplacementSpec::Stochastic { seed: 7 });

        let first_count = first_transformer
            .transform_record(b"chr1", &mut first, 3)
            .expect("first transform");
        let second_count = second_transformer
            .transform_record(b"chr1", &mut second, 3)
            .expect("second transform");

        assert_eq!(first_count, second_count);
        assert_eq!(first, second);
    }
}
