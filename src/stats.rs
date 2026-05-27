use crate::bam::{BasicBamStatistics, index_bam, stats_indexed_bam, stats_xam};
use crate::fastq2::{FastqRecord, RecordEachStats, chunk_records_from_fastq};
use crate::input_type::{InputType, check_input_type};
use crate::summary::{make_plot, write_summary};
use crate::utils::{
    calculate_quality, collect_fqs_in_dir, gc, positive_number_parse, quit_with_error,
};
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::io::Write;
use std::sync::mpsc::Receiver;

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
        let read_q = calculate_quality(quals, use_dorado_q, false);
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
    let fastqs = collect_fqs_in_dir(fastq_dir);
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
    let topn = stats_cmd.get_one::<u32>("topn").unwrap();
    let quality = stats_cmd.get_one::<Vec<f64>>("quality").unwrap();
    let use_dorado_q = stats_cmd.get_flag("use_dorado_q");
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
            )
            .unwrap();
        } else {
            make_plot(
                &basic_stats,
                *quantile,
                plot.unwrap(),
                &formats,
                python,
                output.unwrap(),
            )
            .unwrap();
        }
    }
}

pub fn stats_cmd() -> Command {
    Command::new("stats")
        .about("stats nanopore reads, output stats result, summary and optional figures")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .required(true)
                .help("the input file, could be
    1. a single fastq[.gz]
    2. a directory containing some fastq[.gz]
    3. a bam or sam file")
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("output the stats result into this tsv file if specified. it will be truncated if it exists")
        )
        .arg(
            Arg::new("summary")
                .short('s')
                .long("summary")
                .default_value("./NanofqStatsSummary.txt")
                .help("output stats summary into this file, it will be truncated if it exists")
        )
        .arg(
            Arg::new("topn")
                .short('n')
                .long("topn")
                .default_value("5")
                .value_parser(value_parser!(u32).range(0..100))
                .help("write the top  N longest reads and highest quality reads info into summary file")
        )
        .arg(
            Arg::new("use_dorado_q")
                .short('u')
                .long("use_dorado_q")
                .action(ArgAction::SetTrue)
                .help("use dorado q-score calculation. this means the leading 60 bases will be trimmed if the read length is longer than 60 when calculate the read Q-value [default: false]")
        )
        .arg(
            Arg::new("quality")
                .short('q')
                .long("quality")
                .default_value("25,20,18,15,12,10")
                .value_parser(|x:&str | {
                    let mut qualities = x.split(",")
                        .into_iter()
                        .map(|each| {
                            match each.parse::<f64>() {
                                Ok(qual) => qual,
                                Err(_) => {
                                    quit_with_error("Failed to parse f64 from --quality")
                                }
                            }
                        })
                        .collect::<Vec<f64>>();
                    qualities.sort_by(|a, b| b.partial_cmp(a).unwrap());
                    Result::<Vec<f64>, anyhow::Error>::Ok(qualities)
                })
                .help("count the reads number that whose quality is bigger than this value, multi value can be separated by comma")
        )
        .arg(
            Arg::new("length")
                .short('l')
                .long("length")
                .value_parser(|x: &str| {
                    let mut lengths = x.split(",")
                        .into_iter()
                        .map(|each| {
                            match each.parse::<u32>() {
                                Ok(len) => len,
                                Err(err) => {
                                    eprintln!("{:?}", err);
                                    quit_with_error("Failed to parse usie from --length")
                                }
                            }
                        })
                        .collect::<Vec<u32>>();
                    lengths.sort_by_key(|x| Reverse(*x));
                    Result::<Vec<u32>, anyhow::Error>::Ok(lengths)
                })
                .help("count the reads number that whose length is bigger than this value if you set this parameter, multi values can be separated by comma")
        )
        .arg(
            Arg::new("gc")
                .long("gc")
                .action(ArgAction::SetTrue)
                .help("whether to stats the gc content [default: false]")
        )
        .arg(
            Arg::new("index")
                .short('I')
                .long("index")
                .action(ArgAction::SetTrue)
                .help("build index firstly for sorted bet unindexed bam file [default: false]")
        )
        .arg(
            Arg::new("thread")
                .short('t')
                .long("thread")
                .default_value("1")
                .value_parser(value_parser!(u16).range(1..=32))
                .help("number of threads")
        )
        .arg(
            Arg::new("chunk")
                .short('c')
                .long("chunk")
                .default_value("50000")
                .value_parser(|x: &str| positive_number_parse(x, "--chunk", false, 20000, u32::MAX))
                .help("reads chunk size when multi threads used")
        )
        .arg(
            Arg::new("python")
                .long("python")
                .default_value("python3")
                .help("python3 path, matplotlib will be imported")
        )
        .arg(
            Arg::new("plot")
                .short('p')
                .long("plot")
                .help("whether to make plot, if set, it should be the prefix of figure path without filename extension")
        )
        .arg(
            Arg::new("format")
                .short('f')
                .long("format")
                .action(ArgAction::Append)
                .value_parser(["png", "pdf", "jpg", "svg"])
                .default_value("pdf")
                .help("which format figure do you want if --plot is set, this parameter can be set multi times")
        )
        .arg(
            Arg::new("quantile")
                .long("quantile")
                .default_value("0.01")
                .value_parser(|x:&str| positive_number_parse(x, "--quantile", true, 0.0f64, 0.5f64))
                .help("the shortest ratio and longest ratio of reads will not be rendered on figure, should be in range(0.0, 0.5)")
        )
}

