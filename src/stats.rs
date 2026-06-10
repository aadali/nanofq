use crate::bam::{BasicBamStatistics, index_bam, stats_indexed_bam, stats_xam};
use crate::fastq::{FastqRecord, RecordEachStats, chunk_records_from_fastq};
use crate::input_type::{InputType, check_input_type};
use crate::summary::SummaryStats;
use crate::utils::{
    calculate_quality, check_input, check_output_file, collect_fqs_in_dir, gc, positive_f64_parse,
    quit_with_error,
};
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::cmp::Reverse;
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
    let input_file = stats_cmd.get_one::<String>("input").unwrap();
    let report = stats_cmd.get_one::<String>("report").unwrap();
    let analysis_name = stats_cmd.get_one::<String>("name").unwrap();
    let output = stats_cmd.get_one::<String>("output");
    let summary = stats_cmd.get_one::<String>("summary");
    let topn = stats_cmd.get_one::<u32>("topn").unwrap();
    let quality = stats_cmd.get_one::<Vec<f64>>("quality").unwrap();
    let use_dorado_q = stats_cmd.get_flag("use_dorado_q");
    let lengths = stats_cmd.get_one::<Vec<u32>>("length");
    let use_gc = stats_cmd.get_flag("gc");
    let index = stats_cmd.get_flag("index");
    let thread = stats_cmd.get_one::<u16>("thread").unwrap();
    let chunk = stats_cmd.get_one::<u32>("chunk").unwrap();
    let bins = stats_cmd.get_one::<u32>("bins").unwrap();
    let quantile = stats_cmd.get_one::<f64>("quantile").unwrap();
    // let input_file = input.unwrap();
    let input_t = check_input_type(input_file);
    check_input(input_file);
    check_output_file(report);
    if summary.is_some() {
        check_output_file(summary.unwrap())
    }
    if output.is_some() {
        check_output_file(output.unwrap())
    }

    if thread != &1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(*thread as usize)
            .build_global()
            .unwrap()
    }

    let mut basic_bam_stats = BasicBamStatistics::default();
    let all_stats = match input_t {
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
                stats_indexed_bam(input_file, *thread as usize, use_dorado_q, use_gc);
            basic_bam_stats = basic_bam_stats_;
            all_stats
        }
    };
    let mut stats_summary = SummaryStats::new(
        all_stats,
        lengths.map(|x| x.as_slice()),
        quality,
        use_gc,
        *topn as usize,
    );
    if output.is_some() {
        stats_summary.save_all_stats(analysis_name, output.unwrap());
    }
    if summary.is_some() {
        stats_summary.write_summary_to_text(analysis_name, &basic_bam_stats, summary.unwrap());
    }
    stats_summary.write_to_html_file(analysis_name, *bins as usize, *quantile, report);
}

pub fn stats_cmd() -> Command {
    Command::new("stats")
        .about("stats nanopore reads, output html report and optional stats result and summary file")
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
            Arg::new("report")
                .short('r')
                .long("report")
                .required(true)
                .help("the output html report file")
        )
        .arg(
            Arg::new("name")
                .long("name")
                .default_value("test001")
                .help("this analysis name, will be showed in first line of output, first line of summary and title of the html report")
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
                .help("output stats summary into this file if specified, it will be truncated if it exists")
        )
        .arg(
            Arg::new("topn")
                .short('n')
                .long("topn")
                .default_value("5")
                .value_parser(value_parser!(u32).range(0..100))
                .help("write the top N longest reads and highest quality reads info into summary file")
        )
        .arg(
            Arg::new("use_dorado_q")
                .short('u')
                .long("use_dorado_q")
                .action(ArgAction::SetTrue)
                .help("use dorado q-score calculation. this means the leading 60 bases will be trimmed if the read length is longer than 60 when calculating the read Q-value")
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
                .help("count the reads whose quality is greater than this value, multiple values can be separated by comma")
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
                .help("count reads whose length is greater than this value if you set this parameter, multiple values can be separated by comma")
        )
        .arg(
            Arg::new("gc")
                .long("gc")
                .action(ArgAction::SetTrue)
                .help("whether to calculate the GC content [default: false]")
        )
        .arg(
            Arg::new("index")
                .short('I')
                .long("index")
                .action(ArgAction::SetTrue)
                .help("build index firstly for sorted but unindexed bam file [default: false]")
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
                .value_parser(value_parser!(u32).range(10000..1000001))
                // .value_parser(|x: &str| positive_int_parse(x, "--chunk",  20000, 1000000))
                .help("reads chunk size when multi threads used")
        )
        .arg(
            Arg::new("bins")
                .long("bins")
                .default_value("100")
                .value_parser(value_parser!(u32))
                .help("bins of histogram in html report")
        )
        .arg(
            Arg::new("quantile")
                .long("quantile")
                .default_value("0.01")
                .value_parser(|x:&str| positive_f64_parse(x, "--quantile",  0.0f64, 0.5f64))
                .help("the top quantile of reads lengths will be excluded from the read length distribution in html report")
        )
}
