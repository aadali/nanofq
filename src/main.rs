mod amplicons;
mod bam;
mod fastq;
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
use clap::Command;
use std::time::Instant;


fn main() {

    let start = Instant::now();
    let cmd = Command::new("nanofq")
        .version("0.4.1")
        .about(
            "A simple program for nanopore long reads to stats, generate draft consensus from amplicons, filter, subseq...",
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
