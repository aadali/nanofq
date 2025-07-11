// #![allow(unused_assignments)]
// #![allow(unused_mut)]
// #![allow(dead_code)]
// #![allow(unused_imports)]
// #![allow(unused_variables)]
mod alignment;
mod arguments;
mod fastq;
mod run;
mod summary;
mod trim;
mod utils;

use fastq::{NanoRead};
use rayon::prelude::*;
use std::io::prelude::*;
use crate::run::{run_filter, run_stats};
use std::time::Instant;
use bio::alignment::pairwise::Aligner;
use trim::test_trim;
use alignment::test_alignment;

fn main() -> Result<(), anyhow::Error> {
    // test_trim();
    test_alignment();
    // alignment::pairwise_align::test_local_align();
    // Ok(())
    
    /*
    let start = Instant::now();
    let matches = arguments::parse_arguments();
    if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        run_stats(stats_cmd)?;
    } else if let Some(filter_cmd) = matches.subcommand_matches("filter") {
        run_filter(filter_cmd)?;
    } else {
        println!("error");
    }
    let dur = start.elapsed();
    println!("Elapsed time: {:6?}", dur);
     */
    Ok(())
    
}