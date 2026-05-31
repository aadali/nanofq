mod all_test;
mod amplicons;
mod bam;
mod fastq2;
mod filter;
mod input_type;
mod primer_barcode;
mod stats;
mod subseq;
mod summary;
mod utils;

use crate::amplicons::{amplicons_cmd, run_amplicons};
use crate::filter::{filter_cmd, run_filter};
use crate::stats::{run_stats, stats_cmd};
use crate::subseq::{run_subseq, subseq_cmd};
use chrono::Local;
use clap::Command;
use colored::Colorize;
use env_logger::fmt::Formatter;
use log::Level;
use log::{LevelFilter, Record};
use std::io::Write;
use std::time::Instant;

fn main() {
    colog::default_builder()
        .filter_level(LevelFilter::Info)
        .format(|buf: &mut Formatter, record: &Record| {
            let now = Local::now().format("%Y-%m-%d %H:%M:%S%.3f");
            writeln!(
                buf,
                "[{}] [{}] [{}] {}",
                now,
                match record.level() {
                    Level::Error => {
                        "ERR".red()
                    }
                    Level::Warn => {
                        "WAR".yellow()
                    }
                    Level::Info => {
                        "INF".green()
                    }
                    Level::Debug => {
                        "DBG".green()
                    }
                    Level::Trace => {
                        "TRC".magenta()
                    }
                },
                record.target(),
                record.args()
            )
        })
        .init();

    let start = Instant::now();
    let cmd = Command::new("nanofq")
        .version("0.4.0")
        .about(
            "A simple program for nanopore long reads to stats, get consensus from amplicons, filter, subseq...",
        )
        .arg_required_else_help(true)
        .subcommand(stats_cmd())
        .subcommand(amplicons_cmd())
        .subcommand(filter_cmd())
        .subcommand(subseq_cmd());
    let matches = cmd.get_matches();

    if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        run_stats(stats_cmd);
    } else if let Some(filter_cmd) = matches.subcommand_matches("filter") {
        run_filter(filter_cmd)
    } else if let Some(amplicons_cmd) = matches.subcommand_matches("amplicon") {
        run_amplicons(amplicons_cmd)
    } else if let Some(subseq_cmd) = matches.subcommand_matches("subseq") {
        run_subseq(subseq_cmd)
    } else {
        unreachable!()
    };
    eprintln!("{:?}", start.elapsed());
}
