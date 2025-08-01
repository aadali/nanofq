#![allow(unused_assignments)]
// #![allow(unused_mut)]
#![allow(dead_code)]
// #![allow(unused_imports)]
// #![allow(unused_variables)]
mod alignment;
mod arguments;
mod fastq;
mod run;
mod summary;
mod trim;
mod utils;
mod amplicon;

use crate::run::run_entry::{run_filter, run_stats, run_trim};
use std::time::Instant;

fn main() -> Result<(), anyhow::Error> {
    let start = Instant::now();
    let matches = arguments::parse_arguments();
    if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        run_stats(stats_cmd)?;
    } else if let Some(filter_cmd) = matches.subcommand_matches("filter") {
        run_filter(filter_cmd)?;
    } else if let Some(trim_cmd) = matches.subcommand_matches("trim") {
        run_trim(trim_cmd)?;
    }
    let dur = start.elapsed();
    println!("Elapsed time: {:6?}", dur);
    Ok(())
}
