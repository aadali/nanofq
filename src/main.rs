#![allow(unused_assignments)]
#![allow(unused_mut)]
#![allow(dead_code)]
// #![allow(unused_imports)]
#![allow(unused_variables)]
mod arguments;
mod fastq;
mod stats;

use fastq::{EachStats, FastqReader, ReadStats};
use rayon::prelude::*;
use seq_io::fastq::{Reader, RecordSet};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};
use std::thread::JoinHandle;
use std::{io, thread};

use std::io::BufReader;
use std::time::Instant;

fn stats_one_fastq_multi_threads() {
    // multi threads time with debug: 27.91277025s
    // multi threads time with release: 10.424908583s
    rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global()
        .unwrap();
    let start = Instant::now();
    let fq = "/Users/aadali/test_data/big.fastq";
    // let fq = "/home/a/big/ycq/projects/RustProjects/nanofq/test_data/big.fastq";
    // let chunk = 10000usize;
    let mut all_stats = vec![];
    // let mut fq_reader = FastqReader(seq_io::fastq::Reader::from_path(fq).unwrap());
    let mut fq_reader = FastqReader(seq_io::fastq::Reader::with_capacity(
        std::io::BufReader::new(std::fs::File::open(fq).unwrap()),
        1024 * 1024
    ));
    let mut round = 1;
    let mut record_set = seq_io::fastq::RecordSet::default();
    loop {
        if let None = fq_reader.read_record_set(&mut record_set) {
            println!("read finished");
            break;
        }
        println!("round: {}", round);
        round += 1;

        let record_vec: Vec<seq_io::fastq::RefRecord> = record_set.into_iter().collect();

        all_stats.extend(
            record_vec
                .into_par_iter()
                .map(|x| x.stats(false))
                .collect::<Vec<EachStats>>(),
        );
    }
    
    let output_file = Path::new("/Users/aadali/stats_one_big_fastq_multi_threads_with_release.tsv");
    let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
    for x in &all_stats {
        writer
            .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
            .unwrap()
    }
    writer.flush().unwrap();
    let dur = start.elapsed();
    println!("Elapsed time: {:6?}", dur);
}
fn main() {
    // println!("Hello, world!");
    // let x: Result<(), String> = Err("hello".to_string());
    // x.expect("today is a good day");
    // stats_one_fastq_multi_threads()
    let lengths = "990,2341,23415,432532";
    let x = lengths.split(",")
        .into_iter()
        .map(|x| x.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();
    println!("{:?}", x);
    let x = "234";
    let x = x.parse::<usize>().unwrap();
    println!("{x}");
    // let  y = x.parse::<f64>();
    let matches = arguments::parse_arguments();
    if let Some(stats_cmd) = matches.subcommand_matches("stats") {
        println!("{:?}", stats_cmd);
        // let formats = stats_cmd.get_many::<Vec<String>>("format");
        let formats = stats_cmd.get_one::<String>("format");
        let formats = stats_cmd.get_many::<String>("format")
            .unwrap_or_default()
            .collect::<Vec<&String>>();
        let plot = stats_cmd.get_flag("plot");
        let topn= stats_cmd.get_one::<u16>("topn");
        let quality = stats_cmd.get_one::<Vec<f64>>("quality");
        let lengths = stats_cmd.get_one::<Vec<usize>>("length");
        println!("format: {:?}", formats);
        println!("plot: {:?}", plot);
        println!("topn: {:?}", topn);
        println!("quality: {:?}", quality);
        println!("length: {:?}", lengths);
        
    }
    // let stats_cmd = matches.subcommand_matches("stats")
}
fn read_big_fastq() -> (Receiver<RecordSet>, JoinHandle<io::Result<()>>) {
    let fq = "/Users/aadali/test_data/big.fastq";
    let chunk = 50000usize;
    let (sender, receiver): (Sender<RecordSet>, Receiver<RecordSet>) = mpsc::channel();
    let mut fq_reader = FastqReader(Reader::from_path_with_capacity(fq, 1024 * 1024 * 4).unwrap());
    // let mut fq_reader = FqReader::new(File::open(fq).unwrap());
    let handle = thread::spawn(move || {
        let mut record_set = RecordSet::default();
        fq_reader.read_record_set_exact(&mut record_set, Some(chunk));
        if sender.send(record_set).is_err() {
            return Ok(());
        }
        Ok(())
    });
    (receiver, handle)
}

fn stats_receiver(receiver: Receiver<RecordSet>) -> Vec<EachStats> {
    let mut all_stats: Vec<EachStats> = vec![];
    for record_set in receiver {
        let mut ref_record_vec = vec![];
        for record in record_set.into_iter() {
            ref_record_vec.push(record);
        }
        let chunk_stats: Vec<EachStats> = ref_record_vec
            .into_par_iter()
            .map(|x| x.stats(false))
            .collect();
        chunk_stats.into_iter().for_each(|x| all_stats.push(x));
    }
    all_stats
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fastq::{EachStats, FastqReader};
    use crate::stats::get_summary;
    use seq_io::fastq;
    use std::io;
    use std::sync::mpsc::{Receiver, Sender};
    use std::thread::JoinHandle;
    use std::time::{Duration, Instant};

    #[test]
    #[ignore]
    fn test_stats() {
        let fastq = "/Users/aadali/Downloads/barcode0002.fastq";
        let mut fastq_reader = FastqReader(seq_io::fastq::Reader::from_path(fastq).unwrap());
        let mut x = fastq_reader.stats(false);
        let summary = get_summary(
            &mut x,
            &mut vec![20000usize, 15000, 10000, 9000, 5000],
            &mut vec![10.0, 9.0, 12.0, 15.0, 20.0],
            10,
        );
        println!("{summary}");
    }

    #[test]
    #[ignore]
    fn stats_one_fastq_by_mpsc_multi_threads() {
        // multi threads time with debug: 24.88566475s
        // multi threads time with release: 6.834277375s
        let start = Instant::now();
        let (receiver, handle) = read_big_fastq();
        let all_stats = stats_receiver(receiver);
        let _ = handle.join().unwrap();

        // let output_file =
        //     Path::new("/Users/aadali/stats_one_big_fastq_mpsc_multi_threads_with_debug.tsv");
        let output_file =
            Path::new("/Users/aadali/stats_one_big_fastq_mpsc_multi_threads_with_release.tsv");
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
        for x in &all_stats {
            writer
                .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
                .unwrap()
        }
        writer.flush().unwrap();
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }

    #[test]
    // #[ignore]
    fn stats_one_fastq_multi_threads() {
        // multi threads time with debug: 27.91277025s
        // multi threads time with release: 10.424908583s
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build_global()
            .unwrap();
        let start = Instant::now();
        let fq = "/Users/aadali/test_data/big.fastq";
        let chunk = 20000usize;
        let mut all_stats = vec![];
        let mut fq_reader = FastqReader(seq_io::fastq::Reader::from_path(fq).unwrap());
        // let mut records = vec![];
        // while let Some(record) = fq_reader.next(){
        //     records.push(record.unwrap())
        // }
        loop {
            let mut record_set = seq_io::fastq::RecordSet::default();
            if let None = fq_reader.read_record_set_exact(&mut record_set, Some(chunk)) {
                println!("read finished");
                break;
            }
            let record_vec: Vec<fastq::RefRecord> = record_set.into_iter().collect();
            all_stats.extend(
                record_vec
                    .into_par_iter()
                    .map(|x| x.stats(false))
                    .collect::<Vec<EachStats>>(),
            );
        }
        //
        // 'outer: loop {
        //     let mut this_vec_fastq: Vec<FqRecord> = vec![];
        //     'inner: loop {
        //         for record in &mut fq_reader {
        //             this_vec_fastq.push(record.unwrap());
        //             if this_vec_fastq.len() >= chunk {
        //                 let chunk_stats: Vec<EachStats> = this_vec_fastq
        //                     .into_par_iter()
        //                     .map(|f| f.stats(false))
        //                     .collect();
        //                 chunk_stats.into_iter().for_each(|x| all_stats.push(x));
        //                 break 'inner;
        //             }
        //         }
        //         let chunk_stats: Vec<EachStats> = this_vec_fastq
        //             .into_par_iter()
        //             .map(|f| f.stats(false))
        //             .collect();
        //         chunk_stats.into_iter().for_each(|x| all_stats.push(x));
        //         break 'outer;
        //     }
        // }

        let output_file =
            Path::new("/Users/aadali/stats_one_big_fastq_multi_threads_with_release.tsv");
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
        for x in &all_stats {
            writer
                .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
                .unwrap()
        }
        writer.flush().unwrap();
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }

    #[test]
    // #[ignore]
    fn stats_one_fastq_one_thread() {
        // one thread time with debug: 112.585139916s
        // one thread time with release: 24.186320875s
        let start = Instant::now();
        let fq = "/Users/aadali/test_data/big.fastq";
        let mut fq_reader = FastqReader(seq_io::fastq::Reader::from_path(fq).unwrap());
        let all_stats = fq_reader.stats(false);
        // let output_file = Path::new("/Users/aadali/stats_one_big_fastq_one_thread_with_debug.tsv");
        let output_file =
            Path::new("/Users/aadali/stats_one_big_fastq_one_thread_with_release.tsv");
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
        for x in &all_stats {
            writer
                .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
                .unwrap()
        }
        writer.flush().unwrap();
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }

    #[test]
    // #[ignore]
    fn stats_multi_fastqs_multi_threads() {
        // multi threads time with debug: 21.9353635s
        // multi threads time with release: 6.078449167s
        rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build_global()
            .unwrap();
        let fastqs_path = Path::new("/Users/aadali/test_data/fastqs");
        let mut fastqs = vec![];
        for fs in fastqs_path.read_dir().unwrap() {
            if let Ok(fq) = fs {
                fastqs.push(fq.path())
            }
        }
        let start = Instant::now();

        let all_stats: Vec<_> = fastqs
            .into_par_iter()
            .map(|fq| {
                let mut fq_reader = FastqReader(seq_io::fastq::Reader::with_capacity(
                    std::fs::File::open(fq).unwrap(),
                    1024 * 1024,
                ));
                // let mut fq_reader = FqReader::new(File::open(fq).unwrap());
                fq_reader.stats(false)
            })
            .flatten()
            .collect();

        // let output_file = Path::new("/Users/aadali/multi_threads_test_with_debug.tsv");
        let output_file = Path::new("/Users/aadali/multi_threads_test_with_release.tsv");
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
        for x in &all_stats {
            writer
                .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
                .unwrap()
        }
        writer.flush().unwrap();
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }

    #[test]
    #[ignore]
    fn stats_multi_fastqs_one_thread() {
        // one thread time with debug: 119.728782667s
        // one thread time with release: 24.214846084s
        let fastqs_path = Path::new("/Users/aadali/test_data/fastqs");
        let mut fastqs = vec![];
        for fs in fastqs_path.read_dir().unwrap() {
            if let Ok(fq) = fs {
                fastqs.push(fq.path())
            }
        }
        let start = Instant::now();
        let mut all_stats = vec![];
        for fq in &fastqs {
            // let mut fq_reader = FastqReader(seq_io::fastq::Reader::from_path(fq).unwrap());
            let mut fq_reader = FastqReader(seq_io::fastq::Reader::new(std::io::BufReader::new(
                std::fs::File::open(fq).unwrap(),
            )));
            // let mut fq_reader = FqReader::new(std::fs::File::open(fq).unwrap());
            let stats_result = fq_reader.stats(false);
            stats_result.into_iter().for_each(|x| all_stats.push(x));
        }
        // let output_file = Path::new("/Users/aadali/one_thread_test_with_debug.tsv");
        let output_file = Path::new("/Users/aadali/one_thread_test_with_release.tsv");
        let mut writer = std::io::BufWriter::new(std::fs::File::create(&output_file).unwrap());
        for x in &all_stats {
            writer
                .write_fmt(format_args!("{}\t{}\t{:.6}\n", x.0, x.1, x.2))
                .unwrap()
        }
        writer.flush().unwrap();
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
        ()
    }

    #[test]
    #[ignore]
    fn test_time() {
        let start = Instant::now();
        let time = Duration::from_secs(2);
        thread::sleep(time);
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }
}
