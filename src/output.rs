// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use flate2::{write::GzEncoder, Compression};
use std::{
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    path::PathBuf,
    process,
    time::{SystemTime, UNIX_EPOCH},
};

use crate::{
    cli::OutputTarget,
    error::{GenomeMaskError, Result},
};

/// Output writer that supports plain or gzip-compressed writes.
pub enum OutputWriter {
    /// Plain text output
    Plain(BufWriter<Box<dyn Write>>),
    /// Gzip-compressed output
    Gzip(GzEncoder<BufWriter<Box<dyn Write>>>),
}

impl OutputWriter {
    /// Creates a new output writer.
    ///
    /// # Arguments
    /// * `target` - Output destination (stdout or file path)
    /// * `gzip` - Whether to enable gzip compression
    ///
    /// # Returns
    /// * `Ok(Self)` - New writer instance
    pub fn new(target: &OutputTarget, gzip: bool) -> Result<Self> {
        let sink = open_output_sink(target)?;
        let writer = BufWriter::new(sink);

        if gzip {
            Ok(Self::Gzip(GzEncoder::new(writer, Compression::default())))
        } else {
            Ok(Self::Plain(writer))
        }
    }

    /// Finishes writing and flushes any buffered output.
    ///
    /// For gzip output, this finalizes the gzip stream.
    pub fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => writer
                .flush()
                .map_err(|err| GenomeMaskError::io("cannot flush output", err)),
            Self::Gzip(mut writer) => writer
                .try_finish()
                .map_err(|err| GenomeMaskError::io("cannot finish gzip stream", err)),
        }
    }
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

/// Returns a human-readable description of the output target.
///
/// # Arguments
/// * `target` - Output target to describe
///
/// # Returns
/// * `String` - Description ("stdout" or file path)
pub fn describe_output_target(target: &OutputTarget) -> String {
    match target {
        OutputTarget::Stdout => "stdout".to_string(),
        OutputTarget::Path(path) => path.display().to_string(),
    }
}

/// Creates a temporary file.
///
/// # Arguments
/// * `prefix` - Prefix for the temporary filename
/// * `suffix` - Suffix for the temporary filename
///
/// # Returns
/// * `Ok((PathBuf, File))` - Path and file handle for the temp file
///
/// # Example
/// ```rust,ignore
/// let (path, file) = create_temp_file("genomemask", ".tmp.fa")?;
/// // Creates: /tmp/genomemask-{pid}-{timestamp}-{attempt}.tmp.fa
/// ```
pub fn create_temp_file(prefix: &str, suffix: &str) -> Result<(PathBuf, File)> {
    let temp_dir = std::env::temp_dir();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();

    for attempt in 0..64u32 {
        let path = temp_dir.join(format!(
            "{prefix}-{}-{timestamp}-{attempt}{suffix}",
            process::id()
        ));
        match OpenOptions::new().create_new(true).write(true).open(&path) {
            Ok(file) => return Ok((path, file)),
            Err(err) if err.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(err) => {
                return Err(GenomeMaskError::io(
                    format!("cannot create temporary file {}", path.display()),
                    err,
                ))
            }
        }
    }

    Err(GenomeMaskError::InvalidArgument(
        "cannot allocate a temporary FASTA path".to_string(),
    ))
}

/// Opens an output sink for writing.
///
/// # Arguments
/// * `target` - Output target (stdout or file path)
///
/// # Returns
/// * `Ok(Box<dyn Write>)` - Output sink for writing
///
/// # Example
/// ```rust,ignore
/// let file = File::create("output.fa")?;
/// let sink = open_output_sink(OutputTarget::Path(Path::new("output.fa")))?;
/// ```
fn open_output_sink(target: &OutputTarget) -> Result<Box<dyn Write>> {
    match target {
        OutputTarget::Stdout => Ok(Box::new(std::io::stdout())),
        OutputTarget::Path(path) => {
            if let Some(parent) = path.parent() {
                fs::create_dir_all(parent).map_err(|err| {
                    GenomeMaskError::io(
                        format!("cannot create output directory {}", parent.display()),
                        err,
                    )
                })?;
            }
            let file = File::create(path).map_err(|err| {
                GenomeMaskError::io(format!("cannot create {}", path.display()), err)
            })?;
            Ok(Box::new(file))
        }
    }
}
