// Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

use clap::Parser;
use log::{error, info};
use simple_logger::init_with_level;

use genomemask::{app, cli::Args, error::GenomeMaskError};

fn main() {
    let start = std::time::Instant::now();
    let args = Args::parse();

    if let Err(err) = init_with_level(args.level.into()) {
        eprintln!(
            "{}",
            GenomeMaskError::Logger(format!("failed to initialize logger: {err}"))
        );
        std::process::exit(1);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build()
        .unwrap_or_else(|err| {
            error!("failed to initialize Rayon thread pool: {err}");
            std::process::exit(1);
        });

    if let Err(err) = app::run(args.command) {
        error!("{err}");
        std::process::exit(1);
    }

    info!("finished in {:?}", start.elapsed());
}
