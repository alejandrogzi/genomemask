// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use log::info;
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

use crate::{
    cli::{Command, CommonConfig, MaskArgs, NsArgs, OutputFormat, ReplacementSpec, SelenoArgs},
    commands::{mask::MaskTransformer, ns::NsTransformer, seleno::SelenoTransformer},
    io::{is_twobit_input_path, process_input_to_fasta, ProcessStats, RecordTransformer},
    output::{create_temp_file, describe_output_target, OutputWriter},
    twobit_writer::{
        analyze_sequence_for_twobit, write_two_bit_from_normalized_fasta,
        write_two_bit_from_sequence_store, StoredSequenceRecord,
    },
    Result,
};
use twobit::TwoBitFile;

/// Runs the appropriate command based on CLI input.
///
/// # Arguments
/// * `command` - The parsed command (Ns, Seleno, or Mask)
///
/// # Returns
/// * `Ok(())` - Command completed successfully
/// * `Err(GenomeMaskError)` - If processing failed
///
/// # Example
/// ```rust,ignore
/// let args = Args::parse();
/// run(args.command)?;
/// ```
pub fn run(command: Command) -> Result<()> {
    match command {
        Command::Ns(args) => run_ns(args),
        Command::Seleno(args) => run_seleno(args),
        Command::Mask(args) => run_mask(args),
    }
}

/// Runs the N base replacement command.
fn run_ns(args: NsArgs) -> Result<()> {
    let config = args.into_config()?;
    log_replacement(config.common.replacement_spec);
    let mut transformer = NsTransformer::new(config.common.replacement_spec);

    run_command(
        &config.common,
        &mut transformer,
        "replaced",
        "bases",
        "genomemask-ns",
    )
}

/// Runs the selenocysteinee masking command.
fn run_seleno(args: SelenoArgs) -> Result<()> {
    let config = args.into_config()?;
    log_replacement(config.common.replacement_spec);
    info!(
        "loading selenocysteinee BED3 from {}",
        config.selenocysteine.display()
    );

    let mut transformer =
        SelenoTransformer::from_bed3(&config.selenocysteine, config.common.replacement_spec)?;

    run_command(
        &config.common,
        &mut transformer,
        "masked",
        "codons",
        "genomemask-seleno",
    )
}

/// Runs the mask command with regions.
fn run_mask(args: MaskArgs) -> Result<()> {
    let config = args.into_config()?;
    log_replacement(config.common.replacement_spec);
    info!("loading mask regions from {}", config.regions.display());

    let mut transformer = MaskTransformer::from_regions(
        &config.regions,
        config.selection,
        config.common.replacement_spec,
    )?;

    run_command(
        &config.common,
        &mut transformer,
        "masked",
        "bases",
        "genomemask-mask",
    )
}

/// Routes to the appropriate output format handler.
fn run_command<T: RecordTransformer>(
    common: &CommonConfig,
    transformer: &mut T,
    verb: &str,
    unit: &str,
    temp_prefix: &str,
) -> Result<()> {
    match common.output_format {
        OutputFormat::Fasta | OutputFormat::Stdout => run_to_fasta(common, transformer, verb, unit),
        OutputFormat::TwoBit => run_to_twobit(common, transformer, verb, unit, temp_prefix),
    }
}

/// Processes input and writes to FASTA format.
fn run_to_fasta<T: RecordTransformer>(
    common: &CommonConfig,
    transformer: &mut T,
    verb: &str,
    unit: &str,
) -> Result<()> {
    info!(
        "writing FASTA to {}",
        describe_output_target(&common.output_target)
    );

    let mut writer = OutputWriter::new(&common.output_target, common.gzip)?;
    let stats = process_input_to_fasta(&common.sequence, &mut writer, transformer)?;
    writer.finish()?;
    log_summary(stats, verb, unit);
    Ok(())
}

/// Processes input and writes to 2bit format via a temporary FASTA file.
fn run_to_twobit<T: RecordTransformer>(
    common: &CommonConfig,
    transformer: &mut T,
    verb: &str,
    unit: &str,
    temp_prefix: &str,
) -> Result<()> {
    info!(
        "writing 2bit to {}",
        describe_output_target(&common.output_target)
    );

    if is_twobit_input_path(&common.sequence)? {
        return run_to_twobit_from_twobit_input(common, transformer, verb, unit, temp_prefix);
    }

    let (temp_path, temp_file) = create_temp_file(temp_prefix, ".tmp.fa")?;
    let result = run_to_twobit_with_temp(
        common,
        transformer,
        temp_path.as_path(),
        temp_file,
        verb,
        unit,
    );
    let cleanup_result = fs::remove_file(&temp_path).map_err(|err| {
        crate::error::GenomeMaskError::io(format!("cannot remove {}", temp_path.display()), err)
    });

    match (result, cleanup_result) {
        (Ok(()), Ok(())) => Ok(()),
        (Err(err), _) => Err(err),
        (Ok(()), Err(err)) => Err(err),
    }
}

/// Writes to 2bit format via a temporary FASTA file.
fn run_to_twobit_with_temp<T: RecordTransformer>(
    common: &CommonConfig,
    transformer: &mut T,
    temp_path: &Path,
    temp_file: File,
    verb: &str,
    unit: &str,
) -> Result<()> {
    info!(
        "streaming transformed FASTA through {}",
        temp_path.display()
    );

    let mut temp_writer = BufWriter::new(temp_file);
    let stats = process_input_to_fasta(&common.sequence, &mut temp_writer, transformer)?;
    temp_writer.flush().map_err(|err| {
        crate::error::GenomeMaskError::io(format!("cannot flush {}", temp_path.display()), err)
    })?;

    info!("encoding temporary FASTA as 2bit");
    let mut writer = OutputWriter::new(&common.output_target, false)?;
    write_two_bit_from_normalized_fasta(temp_path, &mut writer)?;
    writer.finish()?;
    log_summary(stats, verb, unit);
    Ok(())
}

/// Processes input and writes to 2bit format via a temporary FASTA file.
fn run_to_twobit_from_twobit_input<T: RecordTransformer>(
    common: &CommonConfig,
    transformer: &mut T,
    verb: &str,
    unit: &str,
    temp_prefix: &str,
) -> Result<()> {
    info!("reading input 2bit directly without temporary FASTA");

    let mut genome = TwoBitFile::open(&common.sequence)
        .map_err(|err| {
            crate::error::GenomeMaskError::InvalidTwoBit(format!(
                "cannot open 2bit file {}: {err}",
                common.sequence.display()
            ))
        })?
        .enable_softmask(true);

    let chrom_names = genome.chrom_names();
    if chrom_names.is_empty() {
        return Err(crate::error::GenomeMaskError::EmptyInput);
    }

    let (temp_path, temp_file) = create_temp_file(temp_prefix, ".tmp.seq")?;
    let result = run_to_twobit_from_twobit_input_with_store(
        &mut genome,
        chrom_names,
        common,
        transformer,
        temp_path.as_path(),
        temp_file,
        verb,
        unit,
    );
    let cleanup_result = fs::remove_file(&temp_path).map_err(|err| {
        crate::error::GenomeMaskError::io(format!("cannot remove {}", temp_path.display()), err)
    });

    match (result, cleanup_result) {
        (Ok(()), Ok(())) => Ok(()),
        (Err(err), _) => Err(err),
        (Ok(()), Err(err)) => Err(err),
    }
}

/// Processes input and writes to 2bit format via a temporary FASTA file.
#[allow(clippy::too_many_arguments)]
fn run_to_twobit_from_twobit_input_with_store<T: RecordTransformer>(
    genome: &mut TwoBitFile<impl std::io::Read + std::io::Seek>,
    chrom_names: Vec<String>,
    common: &CommonConfig,
    transformer: &mut T,
    temp_path: &Path,
    temp_file: File,
    verb: &str,
    unit: &str,
) -> Result<()> {
    let mut temp_writer = BufWriter::new(temp_file);
    let mut stats = ProcessStats::default();
    let mut records = Vec::<StoredSequenceRecord>::with_capacity(chrom_names.len());
    let mut sequence_start = 0usize;

    for chrom_name in chrom_names {
        let header = chrom_name.as_bytes().to_vec();
        let mut sequence = genome
            .read_sequence(&chrom_name, ..)
            .map_err(|err| {
                crate::error::GenomeMaskError::InvalidTwoBit(format!(
                    "cannot read 2bit sequence '{}': {err}",
                    chrom_name
                ))
            })?
            .into_bytes();

        let events = transformer.transform_record(&header, &mut sequence, stats.records)?;
        let mut record = analyze_sequence_for_twobit(&header, &sequence)?;
        record.sequence_start = sequence_start;

        temp_writer.write_all(&sequence).map_err(|err| {
            crate::error::GenomeMaskError::io(format!("cannot write {}", temp_path.display()), err)
        })?;

        sequence_start += sequence.len();
        stats.records += 1;
        stats.events += events;
        records.push(record);
    }

    transformer.finish()?;
    temp_writer.flush().map_err(|err| {
        crate::error::GenomeMaskError::io(format!("cannot flush {}", temp_path.display()), err)
    })?;

    info!("encoding transformed sequence store as 2bit");
    let mut writer = OutputWriter::new(&common.output_target, false)?;
    write_two_bit_from_sequence_store(temp_path, records, &mut writer)?;
    writer.finish()?;
    log_summary(stats, verb, unit);
    Ok(())
}

/// Logs the replacement specification at startup.
fn log_replacement(replacement_spec: ReplacementSpec) {
    match replacement_spec {
        ReplacementSpec::Stochastic { seed } => info!("using stochastic replacement seed {seed}"),
        ReplacementSpec::Fixed(base) => info!("using fixed replacement base {}", base as char),
    }
}

/// Logs a summary of processing statistics.
fn log_summary(stats: ProcessStats, verb: &str, unit: &str) {
    info!(
        "processed {} records and {} {} {}",
        stats.records, verb, stats.events, unit
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Args;
    use clap::Parser;
    use std::{
        fs, process,
        time::{SystemTime, UNIX_EPOCH},
    };
    use twobit::TwoBitFile;

    #[test]
    fn app_writes_twobit_output_from_ns_input() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir =
            std::env::temp_dir().join(format!("genomemask-ns-test-{}-{unique}", process::id()));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input = temp_dir.join("input.fa");
        let output = temp_dir.join("input.masked.2bit");
        fs::write(&input, b">chr1\nACNnTG\n").expect("write FASTA");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "ns",
            "--sequence",
            input.to_str().expect("input path"),
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "G",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "ACGgTG");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }

    #[test]
    fn app_writes_twobit_output_from_seleno_input() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir =
            std::env::temp_dir().join(format!("genomemask-seleno-test-{}-{unique}", process::id()));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input = temp_dir.join("input.fa");
        let bed = temp_dir.join("sites.bed");
        let output = temp_dir.join("input.masked.2bit");
        fs::write(&input, b">chr1\nACTTGACCC\n").expect("write FASTA");
        fs::write(&bed, b"chr1\t4\t7\n").expect("write BED");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "seleno",
            "--sequence",
            input.to_str().expect("input path"),
            "--selenocysteine",
            bed.to_str().expect("bed path"),
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "A",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "ACTAAACCC");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }

    #[test]
    fn app_writes_twobit_output_from_twobit_input_without_fasta_roundtrip() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir = std::env::temp_dir().join(format!(
            "genomemask-direct-twobit-test-{}-{unique}",
            process::id()
        ));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input_fasta = temp_dir.join("input.fa");
        let input_twobit = temp_dir.join("input.2bit");
        let output_twobit = temp_dir.join("input.masked.2bit");
        fs::write(&input_fasta, b">chr1\nACNnTG\n").expect("write FASTA");
        let mut seed_writer = File::create(&input_twobit).expect("create input 2bit");
        write_two_bit_from_normalized_fasta(&input_fasta, &mut seed_writer)
            .expect("seed input 2bit");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "ns",
            "--sequence",
            input_twobit.to_str().expect("input path"),
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "G",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output_twobit)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "ACGgTG");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }

    #[test]
    fn app_writes_twobit_output_from_mask_specific_input() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir =
            std::env::temp_dir().join(format!("genomemask-mask-test-{}-{unique}", process::id()));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input = temp_dir.join("input.fa");
        let regions = temp_dir.join("regions.bed");
        let output = temp_dir.join("input.masked.2bit");
        fs::write(&input, b">chr1\nAACCGGTT\n").expect("write FASTA");
        fs::write(&regions, b"chr1\t2\t6\n").expect("write BED");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "mask",
            "--sequence",
            input.to_str().expect("input path"),
            "--regions",
            regions.to_str().expect("regions path"),
            "--specific",
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "T",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "AATTTTTT");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }

    #[test]
    fn app_writes_twobit_output_from_mask_introns() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir = std::env::temp_dir().join(format!(
            "genomemask-mask-intron-test-{}-{unique}",
            process::id()
        ));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input = temp_dir.join("input.fa");
        let regions = temp_dir.join("regions.bed");
        let output = temp_dir.join("input.masked.2bit");
        fs::write(&input, b">chr1\nAAACCCGGGTTT\n").expect("write FASTA");
        fs::write(
            &regions,
            b"chr1\t0\t9\ttx1\t0\t+\t0\t9\t0,0,0\t2\t3,3,\t0,6,\n",
        )
        .expect("write BED12");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "mask",
            "--sequence",
            input.to_str().expect("input path"),
            "--regions",
            regions.to_str().expect("regions path"),
            "--feature",
            "intron",
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "A",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "AAAAAAGGGTTT");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }

    #[test]
    fn app_writes_twobit_output_from_mask_gtf_introns() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let temp_dir = std::env::temp_dir().join(format!(
            "genomemask-mask-gtf-test-{}-{unique}",
            process::id()
        ));
        fs::create_dir_all(&temp_dir).expect("create temp dir");

        let input = temp_dir.join("input.fa");
        let regions = temp_dir.join("regions.gtf");
        let output = temp_dir.join("input.masked.2bit");
        fs::write(&input, b">chr1\nAAACCCGGGTTT\n").expect("write FASTA");
        fs::write(
            &regions,
            concat!(
                "chr1\tsrc\texon\t1\t3\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
                "chr1\tsrc\texon\t7\t9\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\";\n",
            ),
        )
        .expect("write GTF");

        let args = Args::try_parse_from([
            "genomemask",
            "--threads",
            "1",
            "--level",
            "info",
            "mask",
            "--sequence",
            input.to_str().expect("input path"),
            "--regions",
            regions.to_str().expect("regions path"),
            "--feature",
            "intron",
            "--outdir",
            temp_dir.to_str().expect("temp dir"),
            "--output-format",
            "2bit",
            "--nucleotide",
            "T",
        ])
        .expect("parse args");

        run(args.command).expect("run app");

        let mut genome = TwoBitFile::open(&output)
            .expect("open 2bit")
            .enable_softmask(true);
        let sequence = genome.read_sequence("chr1", ..).expect("read chr1");
        assert_eq!(sequence, "AAATTTGGGTTT");

        fs::remove_dir_all(&temp_dir).expect("remove temp dir");
    }
}
