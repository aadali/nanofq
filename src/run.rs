use crate::fastq::{EachStats, FastqReader, ReadStats};
use crate::stats::{ write_stats, write_summary};
use clap::ArgMatches;
use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use seq_io::fastq::RecordSet;
use std::fs::File;
use std::io::{BufReader, Stdin};
use std::path::Path;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::{thread};

fn stats_receiver(receiver: Receiver<RecordSet>, gc: bool) -> Vec<EachStats> {
    let mut all_stats: Vec<EachStats> = vec![];
    for record_set in receiver {
        let mut record_vec = vec![];
        record_set.into_iter().for_each(|x| record_vec.push(x));
        all_stats.extend(
            record_vec
                .into_par_iter()
                .map(|x| x.stats(gc))
                .collect::<Vec<EachStats>>(),
        );
    }
    all_stats
}

fn stats_stdin(thread: usize, gc: bool) -> Vec<EachStats> {
    let mut reader = FastqReader::<Stdin>::with_stdin();
    if thread == 1 {
        reader.stats(gc)
    } else {
        let (sender, receiver) = mpsc::sync_channel(1000);
        let handle = thread::spawn(move || {
            loop {
                let mut record_set = seq_io::fastq::RecordSet::default();
                if reader.read_record_set(&mut record_set).is_none() {
                    break;
                }
                if sender.send(record_set).is_err() {
                    break;
                }
            }
        });
        stats_receiver(receiver, gc)
    }
}

fn stats_file(file_path: &str, thread: usize, gc: bool) -> Vec<EachStats> {
    let mut reader = FastqReader::<File>::with_fastq(file_path);
    if thread == 1 {
        reader.stats(gc)
    } else {
        let (sender, receiver) = mpsc::sync_channel(1000);
        let handle = thread::spawn(move || {
           loop {
               let mut record_set = seq_io::fastq::RecordSet::default();
               if reader.read_record_set(&mut record_set).is_none() {
                   break;
               }
               if sender.send(record_set).is_err() {
                   break;
               }
           } 
        });
        stats_receiver(receiver, gc)
    }
}

fn stats_file_gz(file_path: &str, thread: usize, gc: bool) -> Vec<EachStats> {
    let mut reader = FastqReader::<MultiGzDecoder<BufReader<File>>>::with_fastq_gz(file_path);
    if thread == 1 {
        reader.stats(gc)
    } else {
        // let (sender, receiver) = mpsc::sync_channel(1000);
        let (sender, receiver) = mpsc::channel();
        let handle = thread::spawn(move || {
            loop {
                let mut record_set = seq_io::fastq::RecordSet::default();
                if reader.read_record_set(&mut record_set).is_none() {
                    break;
                }
                if sender.send(record_set).is_err() {
                    break;
                }
            }
        });
        stats_receiver(receiver, gc)
    }
}

fn stats_fastq_dir(path: &Path, thread: usize, gc: bool) -> Vec<EachStats> {
    assert!(path.is_dir());
    let mut all_stats = vec![];
    let mut fastqs = vec![];
    for fs in path.read_dir().expect(&format!(
        "read directory: {} failed",
        path.to_str().unwrap()
    )) {
        if let Ok(fs) = fs {
            let fs_path = fs.path();
            let fs_path_str = fs_path.to_str().unwrap();
            if fs_path_str.ends_with(".fastq")
                || fs_path_str.ends_with(".fq")
                || fs_path_str.ends_with(".fastq.gz")
                || fs_path_str.ends_with(".fq.gz")
            {
                fastqs.push(fs.path())
            }
        }
    }
    for fq in &fastqs {
        if fq.to_str().unwrap().ends_with(".gz") {
            all_stats.extend(stats_file_gz(fq.to_str().unwrap(), thread, gc))
        } else {
            all_stats.extend(stats_file(fq.to_str().unwrap(), thread, gc))
        }
    }
    all_stats
}

pub fn run_stats(stats_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
    let input = stats_cmd.get_one::<String>("input");
    let output = stats_cmd.get_one::<String>("output");
    let summary = stats_cmd.get_one::<String>("summary").unwrap();
    let topn = stats_cmd.get_one::<u16>("topn").unwrap();
    let mut quality = stats_cmd.get_one::<Vec<f64>>("quality").unwrap();
    let mut lengths = stats_cmd.get_one::<Vec<usize>>("length");
    let gc = stats_cmd.get_flag("gc");
    let thread = stats_cmd.get_one::<u16>("thread").unwrap();
    let plot = stats_cmd.get_flag("plot");
    let format = stats_cmd
        .get_many::<String>("format")
        .unwrap()
        .collect::<Vec<&String>>();
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(*thread as usize)
        .build_global()?;

    let mut stats_result = Vec::<EachStats>::new();
    match input {
        None => stats_result = stats_stdin(*thread as usize, gc),
        Some(input) => {
            let input_path = Path::new(input);
            if input_path.is_file() {
                if input_path.to_str().unwrap().ends_with(".gz") {
                    stats_result = stats_file_gz(input, *thread as usize, gc);
                } else {
                    stats_result = stats_file(input, *thread as usize, gc);
                }
            } else {
                stats_result = stats_fastq_dir(input_path, *thread as usize, gc);
            }
        }
    }
    
    match output {
        None => write_stats(&stats_result, &mut std::io::stdout(), gc)?,
        Some(output_file) => write_stats(
            &stats_result,
            &mut std::io::BufWriter::new(std::fs::File::create(output_file).unwrap()),
            gc,
        )?,
    }
    write_summary(&mut stats_result, lengths, quality, *topn as usize, summary);
    Ok(())
}
