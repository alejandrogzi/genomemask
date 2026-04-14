// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

pub mod app;
pub mod cli;
pub mod commands;
pub mod error;
pub mod io;
pub mod output;
pub mod twobit_writer;

pub use crate::error::{GenomeMaskError, Result};
