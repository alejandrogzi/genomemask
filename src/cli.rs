// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::{ArgAction, Args as ClapArgs, Parser, Subcommand, ValueEnum};
use log::Level;
use std::path::{Path, PathBuf};

use crate::error::{GenomeMaskError, Result};

/// Output format for masked genome sequences.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// FASTA format output
    Fasta,
    /// 2bit binary format output
    #[value(name = "2bit")]
    TwoBit,
    /// Write to stdout
    Stdout,
}

/// Log level for diagnostic output.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl From<LogLevel> for Level {
    fn from(level: LogLevel) -> Self {
        match level {
            LogLevel::Error => Level::Error,
            LogLevel::Warn => Level::Warn,
            LogLevel::Info => Level::Info,
            LogLevel::Debug => Level::Debug,
            LogLevel::Trace => Level::Trace,
        }
    }
}

/// Specification for nucleotide replacement during masking.
///
/// # Variants
/// - `Fixed`: Replace all bases with a single nucleotide
/// - `Stochastic`: Replace bases with random nucleotides using a seed
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ReplacementSpec {
    /// Fixed nucleotide replacement
    Fixed(u8),
    /// Stochastic replacement with seed for determinism
    Stochastic { seed: u64 },
}

/// Output destination for masked sequences.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum OutputTarget {
    /// Write to stdout
    Stdout,
    /// Write to a file path
    Path(PathBuf),
}

/// Feature type for transcript-derived masking.
#[derive(Debug, Clone, Copy, Eq, PartialEq, ValueEnum)]
pub enum MaskFeature {
    /// Coding sequences (CDS)
    Cds,
    /// Exon regions
    Exon,
    /// Intron regions
    Intron,
    /// Untranslated regions (UTR)
    Utr,
}

/// Selection mode for mask intervals.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum MaskSelection {
    /// Use BED intervals directly as mask regions
    Specific,
    /// Derive intervals from transcript features
    Feature(MaskFeature),
}

/// Common configuration shared across all commands.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct CommonConfig {
    pub sequence: PathBuf,
    pub output_format: OutputFormat,
    pub output_target: OutputTarget,
    pub replacement_spec: ReplacementSpec,
    pub gzip: bool,
}

/// Configuration for the `ns` command.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct NsConfig {
    pub common: CommonConfig,
}

/// Configuration for the `seleno` command.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SelenoConfig {
    pub common: CommonConfig,
    /// Path to BED3 file with TGA codon coordinates
    pub selenocysteine: PathBuf,
}

/// Configuration for the `mask` command.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct MaskConfig {
    pub common: CommonConfig,
    /// Path to regions file (BED, GTF, GFF)
    pub regions: PathBuf,
    /// Selection mode for mask intervals
    pub selection: MaskSelection,
}

#[derive(Debug, Parser)]
#[command(
    name = "genomemask",
    about = "Genome masking toolkit",
    version = env!("CARGO_PKG_VERSION"),
    author = env!("CARGO_PKG_AUTHORS"),
)]
pub struct Args {
    #[arg(
        short = 't',
        long = "threads",
        help = "Number of Rayon worker threads",
        value_name = "THREADS",
        default_value_t = num_cpus::get().max(1),
        global = true
    )]
    pub threads: usize,

    #[arg(
        short = 'l',
        long = "level",
        help = "Log level",
        value_name = "LEVEL",
        value_enum,
        default_value_t = LogLevel::Info,
        global = true
    )]
    pub level: LogLevel,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    #[command(about = "Replace genomic N bases in FASTA and 2bit inputs")]
    Ns(NsArgs),
    #[command(about = "Mask selenocysteinee TGA codons in FASTA and 2bit genomes")]
    Seleno(SelenoArgs),
    #[command(about = "Mask bases in direct or transcript-derived genomic regions")]
    Mask(MaskArgs),
}

#[derive(Debug, ClapArgs)]
pub struct NsArgs {
    #[command(flatten)]
    pub common: CommonCommandArgs,
}

impl NsArgs {
    /// Converts CLI arguments to runtime configuration.
    pub fn into_config(self) -> Result<NsConfig> {
        Ok(NsConfig {
            common: self.common.into_config()?,
        })
    }
}

#[derive(Debug, ClapArgs)]
pub struct SelenoArgs {
    #[command(flatten)]
    pub common: CommonCommandArgs,

    #[arg(
        long = "selenocysteine",
        help = "BED3 file containing TGA codon coordinates",
        value_name = "BED3",
        required = true
    )]
    pub selenocysteine: PathBuf,
}

impl SelenoArgs {
    /// Converts CLI arguments to runtime configuration.
    ///
    /// # Arguments
    /// * `self` - The parsed CLI arguments
    ///
    /// # Returns
    /// * `Result<SelenoConfig>` - Configuration or error
    ///
    /// # Example
    /// ```rust,ignore
    /// let args = SelenoArgs::parse_from(["seleno", "--selenocysteine", "sites.bed", ...]);
    /// let config = args.into_config()?;
    /// ```
    pub fn into_config(self) -> Result<SelenoConfig> {
        Ok(SelenoConfig {
            common: self.common.into_config()?,
            selenocysteine: self.selenocysteine,
        })
    }
}

#[derive(Debug, ClapArgs)]
pub struct MaskArgs {
    #[command(flatten)]
    pub common: CommonCommandArgs,

    #[arg(
        long = "regions",
        help = "Path to the regions file (BED, GTF, GFF, .gz variants for text formats)",
        value_name = "REGIONS",
        required = true
    )]
    pub regions: PathBuf,

    #[arg(
        long = "feature",
        help = "Transcript-derived feature to extract from the regions file",
        value_name = "FEATURE",
        value_enum
    )]
    pub feature: Option<MaskFeature>,

    #[arg(
        long = "specific",
        help = "Treat BED records as direct mask intervals instead of deriving transcript features",
        action = ArgAction::SetTrue
    )]
    pub specific: bool,
}

impl MaskArgs {
    /// Converts CLI arguments to runtime configuration.
    ///
    /// Validates that `--specific` and `--feature` are not both specified.
    pub fn into_config(self) -> Result<MaskConfig> {
        let selection = match (self.specific, self.feature) {
            (true, Some(_)) => {
                return Err(GenomeMaskError::InvalidArgument(
                    "--feature cannot be used with --specific".to_string(),
                ))
            }
            (true, None) => MaskSelection::Specific,
            (false, Some(feature)) => MaskSelection::Feature(feature),
            (false, None) => {
                return Err(GenomeMaskError::InvalidArgument(
                    "--feature is required unless --specific is enabled".to_string(),
                ))
            }
        };

        Ok(MaskConfig {
            common: self.common.into_config()?,
            regions: self.regions,
            selection,
        })
    }
}

#[derive(Debug, Clone, ClapArgs)]
pub struct CommonCommandArgs {
    #[arg(
        short = 's',
        long = "sequence",
        help = "Path to the input genome (.fa, .fa.gz, .fna, .fasta, .2bit) or '-' for stdin",
        value_name = "SEQUENCE",
        default_value = "-"
    )]
    pub sequence: PathBuf,

    #[arg(
        short = 'o',
        long = "outdir",
        help = "Directory for file outputs. Ignored when --output-format=stdout.",
        value_name = "OUTDIR",
        default_value = "."
    )]
    pub outdir: PathBuf,

    #[arg(
        short = 'f',
        long = "output-format",
        help = "Output format",
        value_name = "FORMAT",
        value_enum,
        required = true
    )]
    pub output_format: OutputFormat,

    #[arg(
        short = 'n',
        long = "nucleotide",
        help = "Replacement nucleotide when stochastic mode is off",
        value_name = "NUCLEOTIDE",
        value_parser = parse_nucleotide
    )]
    pub nucleotide: Option<u8>,

    #[arg(
        short = 'S',
        long = "stochastic",
        help = "Enable deterministic stochastic replacement",
        action = ArgAction::SetTrue
    )]
    pub stochastic: bool,

    #[arg(
        long = "seed",
        help = "Seed used for deterministic stochastic replacement",
        value_name = "SEED",
        default_value_t = 0
    )]
    pub seed: u64,

    #[arg(
        short = 'z',
        long = "gzip",
        help = "Compress FASTA output with gzip",
        action = ArgAction::SetTrue
    )]
    pub gzip: bool,
}

impl CommonCommandArgs {
    /// Converts common CLI arguments to runtime configuration.
    ///
    /// Validates that `--nucleotide` and `--stochastic` are not both specified.
    /// Also validates output format compatibility with gzip option.
    fn into_config(self) -> Result<CommonConfig> {
        let replacement_spec = match (self.stochastic, self.nucleotide) {
            (true, Some(_)) => {
                return Err(GenomeMaskError::InvalidArgument(
                    "--nucleotide cannot be used with --stochastic".to_string(),
                ))
            }
            (true, None) => ReplacementSpec::Stochastic { seed: self.seed },
            (false, Some(base)) => ReplacementSpec::Fixed(base),
            (false, None) => {
                return Err(GenomeMaskError::InvalidArgument(
                    "--nucleotide is required unless --stochastic is enabled".to_string(),
                ))
            }
        };

        if self.gzip && self.output_format != OutputFormat::Fasta {
            return Err(GenomeMaskError::InvalidArgument(
                "--gzip is only valid with --output-format=fasta".to_string(),
            ));
        }

        let output_target = match self.output_format {
            OutputFormat::Stdout => OutputTarget::Stdout,
            OutputFormat::Fasta | OutputFormat::TwoBit => OutputTarget::Path(derive_output_path(
                &self.sequence,
                &self.outdir,
                self.output_format,
                self.gzip,
            )?),
        };

        Ok(CommonConfig {
            sequence: self.sequence,
            output_format: self.output_format,
            output_target,
            replacement_spec,
            gzip: self.gzip,
        })
    }
}

/// Parses a single nucleotide character (A, T, C, G).
///
/// # Arguments
/// * `raw` - The raw string input
///
/// # Returns
/// * `Ok(u8)` - Uppercase nucleotide byte
/// * `Err(String)` - Error if not a valid nucleotide
fn parse_nucleotide(raw: &str) -> std::result::Result<u8, String> {
    if raw.len() != 1 {
        return Err("nucleotide must be a single base in {A,T,C,G}".to_string());
    }

    match raw.as_bytes()[0].to_ascii_uppercase() {
        b'A' | b'T' | b'C' | b'G' => Ok(raw.as_bytes()[0].to_ascii_uppercase()),
        _ => Err("nucleotide must be one of {A,T,C,G}".to_string()),
    }
}

/// Checks if a path represents stdin/stdout ("-").
///
/// # Arguments
/// * `path` - The path to check
///
/// # Returns
/// * `true` if the path is "-"
pub fn is_stdio_path(path: &Path) -> bool {
    path == Path::new("-")
}

/// Derives an output path from an input path based on the output format.
///
/// # Arguments
/// * `input` - Input file path
/// * `format` - Output format (fasta, 2bit)
/// * `gzip` - Whether to apply gzip compression
///
/// # Returns
/// * `Result<PathBuf>` - Derived output path
///
/// # Example
/// ```rust,ignore
/// let path = derive_output_path(Path::new("genome.fa"), OutputFormat::Fasta, false)?;
/// // Returns: "genome.masked.fa"
/// ```
fn derive_output_path(
    input: &Path,
    outdir: &Path,
    format: OutputFormat,
    gzip: bool,
) -> Result<PathBuf> {
    let stem = if is_stdio_path(input) {
        "stdin"
    } else {
        let file_name = input
            .file_name()
            .and_then(|value| value.to_str())
            .ok_or_else(|| {
                GenomeMaskError::InvalidArgument(format!(
                    "cannot derive an output path from '{}'",
                    input.display()
                ))
            })?;
        strip_known_suffix(file_name)
    };

    let suffix = match format {
        OutputFormat::Fasta if gzip => ".masked.fa.gz",
        OutputFormat::Fasta => ".masked.fa",
        OutputFormat::TwoBit => ".masked.2bit",
        OutputFormat::Stdout => unreachable!("stdout output does not derive a path"),
    };

    Ok(outdir.join(format!("{stem}{suffix}")))
}

/// Strips known file extensions from a filename.
///
/// Removes common genome file suffixes like .fa, .fasta, .fna, .gz, .2bit.
///
/// # Arguments
/// * `file_name` - The filename to strip
///
/// # Returns
/// * `&str` - Filename with known suffix removed
///
/// # Example
/// ```rust,ignore
/// let base = strip_known_suffix("genome.fa.gz"); // Returns "genome"
/// ```
fn strip_known_suffix(file_name: &str) -> &str {
    for suffix in [
        ".fa.gz",
        ".fasta.gz",
        ".fna.gz",
        ".fa",
        ".fasta",
        ".fna",
        ".2bit",
        ".gz",
    ] {
        if let Some(prefix) = file_name.strip_suffix(suffix) {
            return prefix;
        }
    }

    file_name
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[test]
    fn ns_requires_nucleotide_without_stochastic() {
        let result = Args::try_parse_from(["genomemask", "ns", "--output-format", "fasta"]);

        assert!(result.is_ok());
        let args = result.expect("parse args");
        let error = match args.command {
            Command::Ns(command) => command.into_config().expect_err("config error"),
            _ => panic!("expected ns command"),
        };

        assert_eq!(
            error.to_string(),
            "--nucleotide is required unless --stochastic is enabled"
        );
    }

    #[test]
    fn seleno_rejects_nucleotide_with_stochastic() {
        let result = Args::try_parse_from([
            "genomemask",
            "seleno",
            "--output-format",
            "fasta",
            "--selenocysteine",
            "sites.bed",
            "--nucleotide",
            "A",
            "--stochastic",
        ]);

        assert!(result.is_ok());
        let args = result.expect("parse args");
        let error = match args.command {
            Command::Seleno(command) => command.into_config().expect_err("config error"),
            _ => panic!("expected seleno command"),
        };

        assert_eq!(
            error.to_string(),
            "--nucleotide cannot be used with --stochastic"
        );
    }

    #[test]
    fn mask_requires_feature_without_specific() {
        let result = Args::try_parse_from([
            "genomemask",
            "mask",
            "--output-format",
            "fasta",
            "--regions",
            "regions.bed",
            "--nucleotide",
            "A",
        ]);

        assert!(result.is_ok());
        let args = result.expect("parse args");
        let error = match args.command {
            Command::Mask(command) => command.into_config().expect_err("config error"),
            _ => panic!("expected mask command"),
        };

        assert_eq!(
            error.to_string(),
            "--feature is required unless --specific is enabled"
        );
    }

    #[test]
    fn mask_rejects_feature_with_specific() {
        let result = Args::try_parse_from([
            "genomemask",
            "mask",
            "--output-format",
            "fasta",
            "--regions",
            "regions.bed",
            "--nucleotide",
            "A",
            "--feature",
            "exon",
            "--specific",
        ]);

        assert!(result.is_ok());
        let args = result.expect("parse args");
        let error = match args.command {
            Command::Mask(command) => command.into_config().expect_err("config error"),
            _ => panic!("expected mask command"),
        };

        assert_eq!(
            error.to_string(),
            "--feature cannot be used with --specific"
        );
    }

    #[test]
    fn derives_output_path_from_input_basename_and_outdir() {
        let args = Args::try_parse_from([
            "genomemask",
            "ns",
            "--sequence",
            "/tmp/here",
            "--outdir",
            "/tmp/out",
            "--output-format",
            "2bit",
            "--nucleotide",
            "A",
        ])
        .expect("parse args");

        let config = match args.command {
            Command::Ns(command) => command.into_config().expect("config"),
            _ => panic!("expected ns command"),
        };

        assert_eq!(
            config.common.output_target,
            OutputTarget::Path(PathBuf::from("/tmp/out/here.masked.2bit"))
        );
    }

    #[test]
    fn derives_stdin_output_name_when_reading_from_stdin() {
        let args = Args::try_parse_from([
            "genomemask",
            "ns",
            "--outdir",
            "/tmp/out",
            "--output-format",
            "fasta",
            "--nucleotide",
            "C",
        ])
        .expect("parse args");

        let config = match args.command {
            Command::Ns(command) => command.into_config().expect("config"),
            _ => panic!("expected ns command"),
        };

        assert_eq!(
            config.common.output_target,
            OutputTarget::Path(PathBuf::from("/tmp/out/stdin.masked.fa"))
        );
    }
}
