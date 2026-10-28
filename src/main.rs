//! simpletrim: a small, portable NGS read trimmer.

mod adapter;
mod analyze;
mod cli;
mod detect;
mod fastq;
mod filter;
mod io;
mod known_adapters;
mod log;
mod matcher;
mod overlap;
mod platform;
mod polyg;
mod quality;
mod short;
mod stats;

use std::process::ExitCode;

use clap::Parser;
use clap::error::ErrorKind;

use cli::{Cli, Command};
use log::Logger;

/// Errors are reported to the user as messages.
pub type Result<T> = std::result::Result<T, String>;

fn main() -> ExitCode {
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(e) => {
            let usage_error =
                !matches!(e.kind(), ErrorKind::DisplayHelp | ErrorKind::DisplayVersion);
            if usage_error && std::env::args().any(|a| a == "--json") {
                let message = e.kind().to_string();
                let detail = e.to_string();
                eprintln!(
                    "{}",
                    serde_json::json!({ "event": "error", "kind": "usage", "message": message, "detail": detail.trim() })
                );
                return ExitCode::from(2);
            }
            e.exit();
        }
    };

    let log = Logger::new(cli.json, cli.quiet);
    let threads = cli.threads as usize;
    let result = match &cli.command {
        Command::Short(args) => short::run(args, threads, &log),
        Command::Analyze(args) => analyze::run(args, threads, &log),
        Command::Long(_) => {
            Err("the 'long' subcommand (long-read trimming) is not implemented yet".into())
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            log.error(&message);
            ExitCode::FAILURE
        }
    }
}
