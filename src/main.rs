// #![allow(unused_assignments)]
// #![allow(unused_mut)]
// #![allow(dead_code)]
// #![allow(unused_imports)]
// #![allow(unused_variables)]
mod alignment;
mod amplicon;
mod arguments;
mod fastq;
mod run;
mod summary;
mod trim;
mod utils;
mod bam;
mod input_type;
mod subseq;
mod fastq2;
mod stats;
mod filter;

use std::time::Instant;
use crate::subseq::run_subseq;
use crate::filter::run_filter;
use crate::stats::run_stats;

fn main()  {
    let start = Instant::now();
    let matches = arguments::parse_arguments();
     if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        run_stats(stats_cmd);
    } else if let Some(filter_cmd) = matches.subcommand_matches("filter") {
        run_filter(filter_cmd)
    } else if let Some(trim_cmd) = matches.subcommand_matches("trim") {
        // run_trim(trim_cmd)
    } else if let Some(amplicon_cmd) = matches.subcommand_matches("amplicon") {
        // run_amplicon(amplicon_cmd)
    }  else if let Some(subseq_cmd) = matches.subcommand_matches("subseq") {
        run_subseq(subseq_cmd)
    } else {
        unreachable!()
    };
    eprintln!("{:?}", start.elapsed());
}
