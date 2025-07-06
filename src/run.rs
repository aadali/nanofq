use crate::fastq::{EachStats, FastqReader, FilterOption, ReadStats};
use crate::summary::{write_stats, write_summary};
use clap::ArgMatches;
use flate2::bufread::MultiGzDecoder;
use rayon::prelude::*;
use seq_io::fastq::{RecordSet, RefRecord};
use std::any::Any;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread;

fn collect_fastq_dir(path: &Path) -> Result<Vec<PathBuf>, anyhow::Error> {
    assert!(path.is_dir());
    let all_fqs = path
        .read_dir()
        .expect(&format!(
            "read directory: {} failed",
            path.to_str().unwrap()
        ))
        .filter_map(|x| {
            if let Ok(fs) = x {
                let fs_path = fs.path();
                let fs_path_str = fs_path.to_str()?;
                return if fs_path_str.ends_with(".fastq")
                    || fs_path_str.ends_with(".fq")
                    || fs_path_str.ends_with(".fastq.gz")
                    || fs_path_str.ends_with(".fq.gz")
                {
                    Some(fs_path)
                } else {
                    None
                };
            }
            None
        })
        .collect::<Vec<PathBuf>>();
    Ok(all_fqs)
}

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

fn stats<R>(reader: R, thread: usize, gc: bool) -> Vec<EachStats>
where
    R: Read + Send + Any,
{
    if thread == 1 {
        FastqReader::new(reader).stats(gc)
    } else {
        let (sender, receiver) = mpsc::sync_channel(1000);
        let handle = thread::spawn(move || {
            let mut reader = FastqReader::<R>::new(reader);
            loop {
                let mut record_set = RecordSet::default();
                if reader.read_record_set(&mut record_set).is_none() {
                    break;
                }
                if sender.send(record_set).is_err() {
                    break;
                }
            }
            Result::<(), anyhow::Error>::Ok(())
        });
        stats_receiver(receiver, gc)
    }
}

fn stats_fastq_dir(path: &Path, thread: usize, gc: bool) -> Vec<EachStats> {
    let fastqs = collect_fastq_dir(path).unwrap();
    fastqs
        .into_par_iter()
        .map(|fq| {
            if fq.to_str().unwrap().ends_with(".gz") {
                stats(
                    MultiGzDecoder::new(BufReader::new(File::open(fq).unwrap())),
                    thread,
                    gc,
                )
            } else {
                stats(BufReader::new(File::open(fq).unwrap()), thread, gc)
            }
        })
        .flatten()
        .collect::<Vec<EachStats>>()
}

fn filter_receiver(
    receiver: Receiver<RecordSet>,
    fo: &FilterOption,
    writer: &mut dyn Write,
    retain_failed: bool,
    failed_writer: &mut BufWriter<File>,
) -> Result<(), anyhow::Error> {
    for record_set in receiver {
        let mut record_vec = vec![];
        record_set.into_iter().for_each(|x| record_vec.push(x));
        let vec2 = record_vec
            .par_iter()
            .map(|x| x.is_passed(fo))
            .collect::<Vec<bool>>();
        if retain_failed {
            for (ref_record, is_passed) in record_vec.iter().zip(&vec2) {
                if *is_passed {
                    ref_record.write(writer)?;
                } else {
                    ref_record.write(failed_writer)?
                }
            }
        } else {
            for (ref_record, is_passed) in record_vec.iter().zip(&vec2) {
                if *is_passed {
                    ref_record.write(writer)?;
                }
            }
        }
    }
    Ok(())
}

fn filter<R>(
    reader: R,
    thread: usize,
    writer: &mut dyn Write,
    fo: &FilterOption,
    retain_failed: bool,
    failed_writer: &mut BufWriter<File>,
) -> Result<(), anyhow::Error>
where
    R: Read + Send + Any,
{
    if thread == 1 {
        FastqReader::new(reader).filter(writer, fo, retain_failed, failed_writer)?
    } else {
        let (sender, receiver) = mpsc::sync_channel(1000);
        let handle = thread::spawn(move || {
            let mut reader = FastqReader::<R>::new(reader);
            loop {
                let mut record_set = RecordSet::default();
                if reader.read_record_set(&mut record_set).is_none() {
                    break;
                }
                if sender.send(record_set).is_err() {
                    break;
                }
            }
        });
        filter_receiver(receiver, fo, writer, retain_failed, failed_writer)?
    }
    Ok(())
}

fn filter_fastq_dir(
    path: &Path,
    thread: usize,
    writer: &mut dyn Write,
    fo: &FilterOption,
    retain_failed: bool,
    failed_writer: &mut BufWriter<File>,
) -> Result<(), anyhow::Error> {
    let fastqs = collect_fastq_dir(path).unwrap();
    for fq in fastqs {
        if fq.to_str().unwrap().ends_with(".gz") {
            filter(
                MultiGzDecoder::new(BufReader::new(File::open(fq)?)),
                thread,
                writer,
                fo,
                retain_failed,
                failed_writer,
            )?
        } else {
            filter(
                BufReader::new(File::open(fq)?),
                thread,
                writer,
                fo,
                retain_failed,
                failed_writer,
            )?;
        }
    }
    Ok(())
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
        // None => stats_result = stats_stdin(*thread as usize, gc),
        None => stats_result = stats(std::io::stdin(), *thread as usize, gc),
        Some(input) => {
            let input_path = Path::new(input);
            if input_path.is_file() {
                if input_path.to_str().unwrap().ends_with(".gz") {
                    stats_result = stats(
                        MultiGzDecoder::new(BufReader::new(File::open(input)?)),
                        *thread as usize,
                        gc,
                    );
                } else {
                    stats_result = stats(BufReader::new(File::open(input)?), *thread as usize, gc);
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
            &mut std::io::BufWriter::new(File::create(output_file).unwrap()),
            gc,
        )?,
    }
    write_summary(&mut stats_result, lengths, quality, *topn as usize, summary);
    Ok(())
}

pub fn run_filter(filter_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
    let input = filter_cmd.get_one::<String>("input");
    let output = filter_cmd.get_one::<String>("output");
    let min_len = filter_cmd.get_one::<usize>("min_len").unwrap();
    let max_len = filter_cmd.get_one::<usize>("max_len").unwrap();
    let min_qual = filter_cmd.get_one::<f64>("min_qual").unwrap();
    let max_qual = filter_cmd.get_one::<f64>("max_qual").unwrap();
    let gc = filter_cmd.get_flag("gc");
    let min_gc = filter_cmd.get_one::<f64>("min_gc").unwrap();
    let max_gc = filter_cmd.get_one::<f64>("max_gc").unwrap();
    let thread = filter_cmd.get_one::<u16>("thread").unwrap();
    let failed_fq_path = filter_cmd.get_one::<String>("retain_failed");
    let filter_option = FilterOption {
        min_len: *min_len,
        max_len: *max_len,
        min_qual: *min_qual,
        max_qual: *max_qual,
        gc,
        min_gc: *min_gc,
        max_gc: *max_gc,
        retain_failed: failed_fq_path,
    };
    let failed_retain = if failed_fq_path.is_none() {
        false
    } else {
        true
    };
    let mut failed_writer = filter_option.set_failed_fastq_file()?.unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(*thread as usize)
        .build_global()?;

    let mut writer: Box<dyn Write> = if output.is_none() {
        Box::new(BufWriter::new(std::io::stdout()))
    } else {
        Box::new(BufWriter::new(File::create(output.unwrap())?))
    };

    match input {
        None => {
            filter(
                std::io::stdin(),
                *thread as usize,
                &mut writer,
                &filter_option,
                failed_retain,
                &mut failed_writer,
            )?;
        }
        Some(input_path) => {
            let ends_with_gz = input_path.ends_with(".gz");
            let input_path = Path::new(input_path);
            if input_path.is_file() {
                if ends_with_gz {
                    let mut reader = MultiGzDecoder::new(BufReader::new(File::open(input_path)?));
                    filter(
                        reader,
                        *thread as usize,
                        &mut writer,
                        &filter_option,
                        failed_retain,
                        &mut failed_writer,
                    )?
                } else {
                    let mut reader = BufReader::new(File::open(input_path)?);
                    filter(
                        reader,
                        *thread as usize,
                        &mut writer,
                        &filter_option,
                        failed_retain,
                        &mut failed_writer,
                    )?
                }
            } else {
                filter_fastq_dir(
                    input_path,
                    *thread as usize,
                    &mut writer,
                    &filter_option,
                    failed_retain,
                    &mut failed_writer,
                )?
            }
        }
    }
    Ok(())
}
