#![allow(unused_assignments)]
#![allow(unused_mut)]
#![allow(dead_code)]
// #![allow(unused_imports)]
#![allow(unused_variables)]
mod fastq;
mod stats;
mod arguments;

use crate::fastq::FqRecord;
use fastq::{EachStats, FqReader};
use rayon::prelude::*;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};
use std::thread::JoinHandle;
use std::{io, thread};

fn main() {
    println!("Hello, world!");
}
fn read_big_fastq() -> (Receiver<Vec<FqRecord>>, JoinHandle<io::Result<()>>) {
    let fq = "/Users/aadali/test_data/big.fastq";
    let chunk = 50000usize;
    let (sender, receiver): (Sender<Vec<FqRecord>>, Receiver<Vec<FqRecord>>) = mpsc::channel();
    let mut fq_reader = FqReader::new(File::open(fq).unwrap());
    let mut this_vec_fastq  = Vec::<FqRecord>::with_capacity(chunk);
    let handle = thread::spawn(move || {
        for fq_record in fq_reader {
            this_vec_fastq.push(fq_record.unwrap());
            if this_vec_fastq.len() == chunk {
                let mut temp_vec_fastq = Vec::<FqRecord>::with_capacity(chunk);
                std::mem::swap(&mut this_vec_fastq, &mut temp_vec_fastq);
                if sender.send(temp_vec_fastq).is_err() {
                    break;
                }
            }
        }
        if sender.send(this_vec_fastq).is_err() {
            return Ok(());
        }
        Ok(())
    });
    (receiver, handle)
}

fn stats_receiver(receiver: Receiver<Vec<FqRecord>>) -> Vec<EachStats> {
    let mut all_stats: Vec<EachStats> = vec![];
    for vec_fq in receiver {
        let chunk_stats: Vec<EachStats> = vec_fq.into_par_iter().map(|x| x.stats(false)).collect();
        chunk_stats.into_iter().for_each(|x| all_stats.push(x));
    }
    all_stats
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::fastq::{EachStats, FqRecord};
    use crate::stats::get_summary;
    use std::io;
    use std::sync::mpsc::{Receiver, Sender};
    use std::thread::JoinHandle;
    use std::time::{Duration, Instant};

    #[test]
    fn test_stats() {
        let fastq = "/Users/aadali/Downloads/barcode0002.fastq";
        let mut fastq_reader = FqReader::new(File::open(fastq).expect("error when open file"));
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
    fn stats_one_fastq_by_mpsc_multi_threads() {
        // multi threads time with debug: 24.88566475s
        // multi threads time with release: 6.834277375s
        let start =Instant::now();
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
    fn stats_one_fastq_multi_threads() {
        // multi threads time with debug: 27.91277025s
        // multi threads time with release: 10.424908583s
        rayon::ThreadPoolBuilder::new().num_threads(4).build_global().unwrap();
        let start = Instant::now();
        let fq = "/Users/aadali/test_data/big.fastq";
        let chunk = 50000usize;
        let mut all_stats = vec![];
        let mut fq_reader = FqReader::new(File::open(fq).unwrap());
        'outer: loop {
            let mut this_vec_fastq: Vec<FqRecord> = vec![];
            'inner: loop {
                for record in &mut fq_reader {
                    this_vec_fastq.push(record.unwrap());
                    if this_vec_fastq.len() >= chunk {
                        let chunk_stats: Vec<EachStats> =
                            this_vec_fastq.into_par_iter().map(|f| f.stats(false)).collect();
                        chunk_stats.into_iter().for_each(|x| all_stats.push(x));
                        break 'inner;
                    }
                }
                let chunk_stats: Vec<EachStats> =
                    this_vec_fastq.into_par_iter().map(|f| f.stats(false)).collect();
                chunk_stats.into_iter().for_each(|x| all_stats.push(x));
                break 'outer;
            }
        }

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
    fn stats_one_fastq_one_thread() {
        // one thread time with debug: 112.585139916s
        // one thread time with release: 24.186320875s
        let start = Instant::now();
        let fq = "/Users/aadali/test_data/big.fastq";
        let mut fq_reader = FqReader::new(File::open(fq).unwrap());
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
    fn stats_multi_fastqs_multi_threads() {
        // multi threads time with debug: 21.9353635s
        // multi threads time with release: 6.078449167s
        rayon::ThreadPoolBuilder::new().num_threads(4).build_global().unwrap();
        let fastqs_path = Path::new("/Users/aadali/test_data/fastqs");
        let mut fastqs = vec![];
        for fs in fastqs_path.read_dir().unwrap() {
            if let Ok(fq) = fs {
                fastqs.push(fq.path())
            }
        }
        let start = Instant::now();

        let all_stats: Vec<_> = fastqs
            .par_iter()
            .map(|fq| {
                let mut fq_reader = FqReader::new(File::open(fq).unwrap());
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
            let mut fq_reader = FqReader::new(std::fs::File::open(fq).unwrap());
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
    fn test_time() {
        let start = Instant::now();
        let time = Duration::from_secs(2);
        thread::sleep(time);
        let dur = start.elapsed();
        println!("Elapsed time: {:6?}", dur);
    }
}
