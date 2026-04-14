// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use std::{error::Error as StdError, fmt, io};

/// Result type alias for GenomeMask operations.
pub type Result<T> = std::result::Result<T, GenomeMaskError>;

/// Error type for GenomeMask operations.
///
/// # Variants
/// - `EmptyInput`: Input data is empty
/// - `InvalidArgument`: Invalid command-line argument
/// - `InvalidBed`: Invalid BED file format
/// - `InvalidFasta`: Invalid FASTA file format
/// - `InvalidRegions`: Invalid regions file format
/// - `InvalidSelenocysteine`: Invalid selenocysteine site specification
/// - `InvalidTwoBit`: Invalid 2bit file format
/// - `InvalidTwoBitHeader`: Invalid 2bit header
/// - `InvalidTwoBitBase`: Invalid nucleotide in 2bit conversion
/// - `Io`: I/O error with context
/// - `Logger`: Logger initialization error
/// - `UnsupportedInput`: Unsupported input format
#[derive(Debug)]
pub enum GenomeMaskError {
    EmptyInput,
    InvalidArgument(String),
    InvalidBed(String),
    InvalidFasta(String),
    InvalidRegions(String),
    InvalidSelenocysteine(String),
    InvalidTwoBit(String),
    InvalidTwoBitHeader {
        header: String,
        reason: String,
    },
    InvalidTwoBitBase {
        header: String,
        offset: usize,
        byte: u8,
    },
    Io {
        context: String,
        source: io::Error,
    },
    Logger(String),
    UnsupportedInput(String),
}

impl GenomeMaskError {
    /// Creates a new I/O error with context.
    ///
    /// # Arguments
    /// * `context` - Description of the operation that failed
    /// * `source` - The underlying I/O error
    ///
    /// # Example
    /// ```rust,ignore
    /// let err = GenomeMaskError::io("cannot open file", io::Error::new(...));
    /// ```
    pub fn io(context: impl Into<String>, source: io::Error) -> Self {
        Self::Io {
            context: context.into(),
            source,
        }
    }
}

impl fmt::Display for GenomeMaskError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "input is empty"),
            Self::InvalidArgument(message) => write!(f, "{message}"),
            Self::InvalidBed(message) => write!(f, "{message}"),
            Self::InvalidFasta(message) => write!(f, "{message}"),
            Self::InvalidRegions(message) => write!(f, "{message}"),
            Self::InvalidSelenocysteine(message) => write!(f, "{message}"),
            Self::InvalidTwoBit(message) => write!(f, "{message}"),
            Self::InvalidTwoBitHeader { header, reason } => {
                write!(f, "cannot write 2bit record '{header}': {reason}")
            }
            Self::InvalidTwoBitBase {
                header,
                offset,
                byte,
            } => write!(
                f,
                "cannot write 2bit record '{header}': unsupported nucleotide {} at sequence offset {}",
                display_byte(*byte),
                offset
            ),
            Self::Io { context, source } => write!(f, "{context}: {source}"),
            Self::Logger(message) => write!(f, "{message}"),
            Self::UnsupportedInput(message) => write!(f, "{message}"),
        }
    }
}

impl StdError for GenomeMaskError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match self {
            Self::Io { source, .. } => Some(source),
            _ => None,
        }
    }
}

fn display_byte(byte: u8) -> String {
    if byte.is_ascii_graphic() || byte == b' ' {
        format!("'{}' (0x{byte:02X})", byte as char)
    } else {
        format!("0x{byte:02X}")
    }
}
