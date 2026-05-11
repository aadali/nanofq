use crate::bam::{BasicBamStatistics, index_bam, stats_indexed_bam, stats_xam};
use crate::fastq2::{FastqRecord, RecordEachStats, chunk_records_from_fastq};
use crate::input_type::{InputType, check_input_type};
use crate::summary::{ write_summary, make_plot};
use crate::utils::{calculate_read_q, collect_fastq_dir, gc, quit_with_error};
use clap::ArgMatches;
use rayon::prelude::*;
use std::io::Write;
use std::sync::mpsc::Receiver;
use std::time::Instant;
use needletail::{parse_fastx_file, Sequence};

fn stats_receiver(
    receiver: Receiver<Vec<FastqRecord>>,
    use_dorado_q: bool,
    use_gc: bool,
) -> Vec<RecordEachStats> {
    let mut all_stats = vec![];
    for records in receiver {
        all_stats.extend(
            records
                .into_par_iter()
                .map(|x| x.stats(use_dorado_q, use_gc))
                .collect::<Vec<RecordEachStats>>(),
        )
    }
    all_stats
}

pub fn fastq_stats(fastq_file: &str, use_dorado_q: bool, use_gc: bool) -> Vec<RecordEachStats> {
    let mut v = vec![];
    let mut records = parse_fastx_file(fastq_file).expect(&format!("Failed to read {fastq_file}"));
    let mut read_idx = 1;
    while let Some(Ok(record)) = records.next() {
        let mut headers = record.id().splitn(2, |x| x.is_ascii_whitespace());
        let name = headers
            .next()
            .expect(&format!("Parse read name failed at {read_idx}th record"));
        let name =
            str::from_utf8(name).expect(&format!("Parse read name failed at {read_idx}th record"));
        let seq = record.sequence();
        let quals = record
            .qual()
            .expect(&format!("Parse quality failed at {read_idx}th record"));
        let read_q = calculate_read_q(quals, use_dorado_q);
        v.push(RecordEachStats::new(
            name,
            record.num_bases(),
            read_q,
            if use_gc { Some(gc(seq)) } else { None },
        ));
        read_idx += 1;
    }
    v
}
fn stats_one_fastq(
    fastq_file: &str,
    thread: usize,
    use_dorado_q: bool,
    use_gc: bool,
    chunk: u32,
) -> Vec<RecordEachStats> {
    if thread == 1 {
        fastq_stats(fastq_file, use_dorado_q, use_gc)
    } else {
        let (read_handle, receiver) = chunk_records_from_fastq(fastq_file, chunk, false);
        let x = stats_receiver(receiver, use_dorado_q, use_gc);
        read_handle.join().unwrap();
        x
    }
}

fn stats_fastq_dir(
    fastq_dir: &str,
    thread: usize,
    use_dorado_q: bool,
    use_gc: bool,
) -> Vec<RecordEachStats> {
    let fastqs = collect_fastq_dir(fastq_dir);
    if thread == 1 {
        fastqs
            .into_iter()
            .map(|x| fastq_stats(x.to_str().unwrap(), use_dorado_q, use_gc))
            .flatten()
            .collect()
    } else {
        fastqs
            .into_par_iter()
            .map(|x| fastq_stats(x.to_str().unwrap(), use_dorado_q, use_gc))
            .flatten()
            .collect()
    }
}

pub fn run_stats(stats_cmd: &ArgMatches) {
    let input = stats_cmd.get_one::<String>("input");
    let output = stats_cmd.get_one::<String>("output");
    let summary = stats_cmd.get_one::<String>("summary").unwrap();
    let topn = stats_cmd.get_one::<u16>("topn").unwrap();
    let quality = stats_cmd.get_one::<Vec<f64>>("quality").unwrap();
    let use_dorado_q = stats_cmd.get_flag("use_dorado_quality");
    let lengths = stats_cmd.get_one::<Vec<u32>>("length");
    let use_gc = stats_cmd.get_flag("gc");
    // let bam = stats_cmd.get_flag("bam");
    let index = stats_cmd.get_flag("index");
    let thread = stats_cmd.get_one::<u16>("thread").unwrap();
    let chunk = stats_cmd.get_one::<u32>("chunk").unwrap();
    let plot = stats_cmd.get_one::<String>("plot");
    let python = stats_cmd.get_one::<String>("python").unwrap();
    let quantile = stats_cmd.get_one::<f64>("quantile").unwrap();
    let format = stats_cmd
        .get_many::<String>("format")
        .unwrap()
        .collect::<Vec<&String>>();
    let input_file = input.unwrap();
    let input_t = check_input_type(input_file);

    if thread != &1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(*thread as usize)
            .build_global()
            .unwrap()
    }

    let mut basic_bam_stats = BasicBamStatistics::default();
    let all_stats = match input_t {
        // InputType::FastqFromStdin => {}
        // InputType::OneBamOrSamFromStdin => {}
        InputType::DirectoryContainFastqsOrFastqsGzipped => {
            stats_fastq_dir(input_file, *thread as usize, use_dorado_q, use_gc)
        }
        InputType::OneFastqFile | InputType::OneFastqGzippedFile => {
            stats_one_fastq(input_file, *thread as usize, use_dorado_q, use_gc, *chunk)
        }
        InputType::OneSamFile | InputType::UnsortedBam | InputType::UnalignedBam => {
            let mut bam_reader = rust_htslib::bam::Reader::from_path(input_file)
                .expect(&format!("Failed to read {}", input_file));
            let (basic_bam_stats_, all_stats) =
                stats_xam(&mut bam_reader, *thread as usize, use_gc, use_dorado_q);
            basic_bam_stats = basic_bam_stats_;
            all_stats
        }
        InputType::SortedUnindexedBam => {
            let (basic_bam_stats_, all_stats) = if index {
                index_bam(input_file, *thread as usize)
                    .expect(&format!("Failed to index {}", input_file));
                // todo adjust stats bam logic for dorado quality
                stats_indexed_bam(input_file, *thread as usize, use_dorado_q, use_gc)
            } else {
                let mut bam_reader = rust_htslib::bam::Reader::from_path(input_file)
                    .expect(&format!("Failed to read {}", input_file));
                stats_xam(&mut bam_reader, *thread as usize, use_dorado_q, use_gc)
            };
            basic_bam_stats = basic_bam_stats_;
            all_stats
        }
        InputType::IndexedBam => {
            let (basic_bam_stats_, all_stats) =
                stats_indexed_bam(input.unwrap(), *thread as usize, use_dorado_q, use_gc);
            basic_bam_stats = basic_bam_stats_;
            all_stats
        }
        _ => quit_with_error("error"),
    };
    // get_summary2(all_stats, )
    let tmp_stats_outfile = format!("/tmp/NanofqStatsTmpResult_{}.tsv", uuid::Uuid::new_v4());
    match output {
        None => {
            if plot.is_some() {
                let writer = std::fs::File::create(&tmp_stats_outfile);
                let mut writer = std::io::BufWriter::new(
                    writer.expect(&format!("Failed to open {tmp_stats_outfile}")),
                );
                for x in &all_stats {
                    writeln!(
                        &mut writer,
                        "{}\t{}\t{:.2}{}",
                        x.name,
                        x.length,
                        x.qual,
                        if use_gc {
                            format!("\t{:.2}", x.gc.unwrap())
                        } else {
                            String::default()
                        }
                    )
                    .unwrap();
                }
            }
        }
        Some(output_file) => {
            let output_file =
                std::fs::File::create(output_file).expect(&format!("Failed to open {output_file}"));
            let mut writer = std::io::BufWriter::new(output_file);
            for x in &all_stats {
                writeln!(
                    &mut writer,
                    "{}\t{}\t{:.2}{}",
                    x.name,
                    x.length,
                    x.qual,
                    if use_gc {
                        format!("\t{:.2}", x.gc.unwrap())
                    } else {
                        String::default()
                    }
                )
                .unwrap();
            }
        }
    }
    let basic_stats = write_summary(
        all_stats,
        lengths.map(|x| x.as_slice()),
        quality,
        *topn as usize,
        &basic_bam_stats,
        summary,
    );
        let formats = format
            .iter()
            .map(|x| (**x).clone())
            .collect::<Vec<String>>();
        if plot.is_some() {
            if output.is_none() {
                make_plot(
                    &basic_stats,
                    *quantile,
                    plot.unwrap(),
                    &formats,
                    python,
                    &tmp_stats_outfile,
                ).unwrap();
            } else {
                make_plot(
                    &basic_stats,
                    *quantile,
                    plot.unwrap(),
                    &formats,
                    python,
                    output.unwrap(),
                ).unwrap();
            }
        }

}

pub fn t() {
    let start = Instant::now();
    let x = stats_one_fastq(
        "/Users/aadali/projects/RustProjects/nanoamp/test_data/ont-barcode05.fastq",
        4,
        false,
        true,
        50000,
    );

    let basic = BasicBamStatistics::default();
    // for a in &x {
    //     println!("{a}")
    // }
    let dur = start.elapsed();
    println!("{dur:?}");
    let start2 = Instant::now();
    write_summary(
        x,
        None,
        &[25.0, 20.0, 18.0, 15.0, 12.0, 10.0],
        5,
        &basic,
        "/Users/aadali/Documents/summary.tsv"
    );
    let dur = start2.elapsed();
    println!("{dur:?}");
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn t1() {
        t()
    }
}
