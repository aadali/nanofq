#![allow(unused_assignments)]
#![allow(unused_mut)]
#![allow(dead_code)]
#![allow(unused_imports)]
#![allow(unused_variables)]
mod alignment;
mod arguments;
mod fastq;
mod run;
mod summary;
mod trim;

use fastq::{ ReadStats};
use rayon::prelude::*;
use std::io::prelude::*;
use crate::run::{run_filter, run_stats};
use std::time::Instant;

fn main() -> Result<(), anyhow::Error> {
    // alignment::pairwise_align::test_local_align();
    // Ok(())
    let start = Instant::now();
    let matches = arguments::parse_arguments();
    if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        // println!("{:?}", stats_cmd);
        // let formats = stats_cmd.get_many::<Vec<String>>("format");
        // let formats = stats_cmd.get_one::<String>("format");
        // let formats = stats_cmd.get_many::<String>("format")
        //     .unwrap_or_default()
        //     .collect::<Vec<&String>>();
        // let plot = stats_cmd.get_flag("plot");
        // let topn= stats_cmd.get_one::<u16>("topn");
        // let quality = stats_cmd.get_one::<Vec<f64>>("quality");
        // let lengths = stats_cmd.get_one::<Vec<usize>>("length");
        // let input = stats_cmd.get_one::<String>("input");
        run_stats(stats_cmd)?;
    } else if let Some(filter_cmd) = matches.subcommand_matches("filter") {
        run_filter(filter_cmd)?;
    } else {
        println!("error");
    }
    let dur = start.elapsed();
    println!("Elapsed time: {:6?}", dur);
    Ok(())
    
}

#[test]
// #[ignore]
fn get_seq_info() {
    // use crate::trim::adapter::get_seq_info;
    // let x = get_seq_info();
    // println!("{:?}", x);
    let x = 18446744073709551615usize;
    let x = "safaksfjdakljfslajfasdlj".as_bytes();
    let sub_x = &x[10..];
    // let y = 18446744073709552000usize;
    println!("{}", std::str::from_utf8(sub_x).unwrap());
}