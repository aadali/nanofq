use crate::fastq::{FastqRecord, chunk_records_from_fastq};
use crate::input_type::{InputType, check_input_type};
use crate::utils::{calculate_quality, check_output_file, collect_fqs_in_dir, gc, positive_f64_parse, quit_with_error};
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use needletail::parser::{LineEnding, write_fastq};
use needletail::{Sequence, parse_fastx_file};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::sync::mpsc::Receiver;

#[derive(Debug, Clone, Default)]
pub struct FilterOption {
    pub min_len: u32,
    pub max_len: u32,
    pub min_qual: f32,
    pub max_qual: f32,
    pub use_dorado_q: bool,
    pub use_gc: bool,
    pub min_gc: f32,
    pub max_gc: f32,
    retain_failed: Option<String>,
}

impl FilterOption {
    fn set_failed_fastq_file(&self) -> Box<dyn Write> {
        match &self.retain_failed {
            None => Box::new(std::io::sink()),
            Some(failed_fastq_file) => {
                let file = std::fs::File::create(&failed_fastq_file).expect(&format!(
                    "Failed to create failed fastq file: {}",
                    &failed_fastq_file
                ));
                Box::new(BufWriter::new(file))
            }
        }
    }
}

fn fastq_filter(fastq_file: &str, fo: &FilterOption, passed_file: &str) {
    let retain_failed = fo.retain_failed.is_some();
    let use_gc = fo.use_gc;
    let mut failed_writer = fo.set_failed_fastq_file();
    let mut passed_writer = BufWriter::new(
        std::fs::File::create(passed_file).expect(&format!("Failed to create {passed_file}")),
    );
    let mut records = parse_fastx_file(fastq_file).expect(&format!("Failed to read {fastq_file}"));
    let mut read_idx = 1;

    while let Some(Ok(record)) = records.next() {
        // check length
        let read_len = record.num_bases() as u32;
        if read_len > fo.max_len || read_len < fo.min_len {
            if retain_failed {
                write_fastq(
                    record.id(),
                    record.sequence(),
                    record.qual(),
                    &mut failed_writer,
                    LineEnding::Unix,
                )
                .expect(&format!(
                    "Failed write {}th record into failed fastq file",
                    read_idx,
                ));
            }
            continue;
        }

        // check quality
        let quals = record
            .qual()
            .expect(&format!("Parse quality failed at {read_idx}th record"));
        let read_qual = calculate_quality(quals, fo.use_dorado_q, false);
        if read_qual > fo.max_qual || read_qual < fo.min_qual {
            println!("{}", str::from_utf8(record.id()).unwrap());
            if retain_failed {
                write_fastq(
                    record.id(),
                    record.sequence(),
                    record.qual(),
                    &mut failed_writer,
                    LineEnding::Unix,
                )
                .expect(&format!(
                    "Failed write {}th record into failed fastq file",
                    read_idx,
                ));
            }
            continue;
        }

        // check gc
        if use_gc {
            let gc = gc(record.sequence());
            if gc > fo.max_gc || gc < fo.min_gc {
                if retain_failed {
                    write_fastq(
                        record.id(),
                        record.sequence(),
                        record.qual(),
                        &mut failed_writer,
                        LineEnding::Unix,
                    )
                    .expect(&format!(
                        "Failed write {}th record into failed fastq file",
                        read_idx,
                    ));
                }
                continue;
            }
        }

        write_fastq(
            record.id(),
            record.sequence(),
            record.qual(),
            &mut passed_writer,
            LineEnding::Unix,
        )
        .expect(&format!(
            "Failed write {}th record into failed fastq file",
            read_idx,
        ));
        read_idx += 1
    }
}

fn fastq_filter_out_records(
    fastq_file: &str,
    fo: &FilterOption,
) -> (Vec<FastqRecord>, Vec<FastqRecord>) {
    let retain_failed = fo.retain_failed.is_some();
    let mut failed_records = vec![];
    let mut passed_records = vec![];
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

        let mut is_passed = true;
        // check length
        let read_length = record.num_bases() as u32;
        if read_length > fo.max_len || read_length < fo.min_len {
            is_passed = false;
            if !retain_failed {
                continue;
            }
        }

        // check quality
        let quals = record
            .qual()
            .expect(&format!("Parse quality failed at {read_idx}th record"));
        let read_qual = calculate_quality(quals, fo.use_dorado_q, false);
        if read_qual > fo.max_qual || read_qual < fo.min_qual {
            is_passed = false;
            if !retain_failed {
                continue;
            }
        }

        // check gc
        if fo.use_gc {
            let gc_content = gc(seq);
            if gc_content > fo.max_gc || gc_content < fo.min_gc {
                is_passed = false;
                if !retain_failed {
                    continue;
                }
            }
        }

        let description = headers
            .next()
            .map(|x| str::from_utf8(x).unwrap_or_default());
        let record = FastqRecord::new(name, description, seq, quals);
        if is_passed {
            passed_records.push(record)
        } else {
            if retain_failed {
                failed_records.push(record)
            }
        }
        read_idx += 1;
    }
    (passed_records, failed_records)
}

fn filter_receiver(
    receiver: Receiver<Vec<FastqRecord>>,
    fo: &FilterOption,
    passed_writer: &mut dyn Write,
    failed_writer: &mut dyn Write,
) {
    for records in receiver {
        let x = records
            .par_iter()
            .map(|x| x.is_passed(&fo))
            .collect::<Vec<_>>();
        for (record, is_passed) in records.iter().zip(x) {
            if is_passed {
                record.write(passed_writer).unwrap();
            } else {
                if fo.retain_failed.is_some() {
                    record.write(failed_writer).unwrap();
                }
            }
        }
    }
}

fn filter_one_fastq(
    fastq_file: &str,
    passed_file: &str,
    thread: usize,
    fo: &FilterOption,
    chunk: u32,
) {
    if thread == 1 {
        fastq_filter(fastq_file, fo, passed_file)
    } else {
        let (read_handle, receiver) = chunk_records_from_fastq(fastq_file, chunk, true);
        let mut failed_writer = fo.set_failed_fastq_file();
        let mut passed_writer = BufWriter::new(
            std::fs::File::create(passed_file).expect(&format!("Failed to create {passed_file}")),
        );
        filter_receiver(receiver, fo, &mut passed_writer, &mut failed_writer);
        read_handle.join().unwrap();
    }
}

fn filter_fastq_dir(fastq_dir: &str, passed_file: &str, thread: usize, fo: &FilterOption) {
    let fastqs = collect_fqs_in_dir(fastq_dir);
    let mut failed_writer = fo.set_failed_fastq_file();
    let mut passed_writer = BufWriter::new(
        std::fs::File::create(passed_file).expect(&format!("Failed to create {passed_file}")),
    );
    if thread == 1 {
        for fq in fastqs {
            let (passed, failed) = fastq_filter_out_records(fq.to_str().unwrap(), fo);
            for p in passed {
                p.write(&mut passed_writer).unwrap()
            }

            for f in failed {
                f.write(&mut failed_writer).unwrap()
            }
        }
    } else {
        let classed_records = fastqs
            .into_par_iter()
            .map(|x| fastq_filter_out_records(x.to_str().unwrap(), fo))
            .collect::<Vec<_>>();
        for (passed, failed) in classed_records {
            for p in passed {
                p.write(&mut passed_writer).unwrap()
            }
            for f in failed {
                f.write(&mut failed_writer).unwrap()
            }
        }
    }
}

pub fn run_filter(filter_cmd: &ArgMatches) {
    let input = filter_cmd.get_one::<String>("input");
    let output = filter_cmd.get_one::<String>("output");
    let min_len = filter_cmd.get_one::<u32>("min_len").unwrap();
    let max_len = filter_cmd.get_one::<u32>("max_len").unwrap();
    let use_dorado_q = filter_cmd.get_flag("use_dorado_q");
    let min_qual = filter_cmd.get_one::<f64>("min_qual").unwrap();
    let max_qual = filter_cmd.get_one::<f64>("max_qual").unwrap();
    let use_gc = filter_cmd.get_flag("gc");
    let min_gc = filter_cmd.get_one::<f64>("min_gc").unwrap();
    let max_gc = filter_cmd.get_one::<f64>("max_gc").unwrap();
    let thread = filter_cmd.get_one::<u16>("thread").unwrap();
    let chunk = filter_cmd.get_one::<u32>("chunk").unwrap();
    // let max_bases = filter_cmd.get_one::<u64>("max_bases");
    let failed_fq_path = filter_cmd.get_one::<String>("retain_failed");
    if failed_fq_path.is_some(){
        check_output_file(failed_fq_path.unwrap())
    }
    let fo = FilterOption {
        min_len: *min_len,
        max_len: *max_len,
        min_qual: *min_qual as f32,
        max_qual: *max_qual as f32,
        use_dorado_q,
        use_gc,
        min_gc: *min_gc as f32,
        max_gc: *max_gc as f32,
        retain_failed: failed_fq_path.map(|x| x.clone()),
    };

    let passed_file = output.unwrap();
    check_output_file(passed_file);
    let input_path = input.unwrap();
    let input_t = check_input_type(input.unwrap());
    if thread == &1 {
        match input_t {
            InputType::DirectoryContainFastqsOrFastqsGzipped => {
                filter_fastq_dir(input_path, passed_file, 1, &fo);
            }
            InputType::OneFastqFile | InputType::OneFastqGzippedFile => {
                filter_one_fastq(input_path, passed_file, 1, &fo, *chunk);
            }
            _ => quit_with_error("error"),
        }
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(*thread as usize)
            .build_global()
            .unwrap();
        match input_t {
            InputType::DirectoryContainFastqsOrFastqsGzipped => {
                filter_fastq_dir(input_path, passed_file, *thread as usize, &fo);
            }
            InputType::OneFastqFile | InputType::OneFastqGzippedFile => {
                filter_one_fastq(input_path, passed_file, *thread as usize, &fo, *chunk);
            }
            _ => quit_with_error("error"),
        }
    }
}

pub fn filter_cmd() -> Command {
    Command::new("filter")
        .about("filter nanopore reads by length, quality or optional gc content")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .required(true)
                .help("the input fastq, a fastq[.gz] or a directory containing some fastq[.gz]"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_parser(|output: &str| {
                    if !(output.ends_with(".fq") || output.ends_with(".fastq")) {
                        quit_with_error( "Error: output should ends with .fastq or .fq. Gzipped not supported", )
                    }
                    Result::<String, anyhow::Error>::Ok(output.to_string())
                })
                .help("output the filtered fastq into this file, it will be truncated if it exists. Compressed file is not supported")
        )
        .arg(
            Arg::new("min_len")
                .short('l')
                .long("min_len")
                .default_value("1")
                .value_parser(value_parser!(u32).range(1..4294967295))
                // .value_parser(|x: &str| positive_int_parse(x, "--min_len",1, u32::MAX as usize))
                .help("min read length")
        )
        .arg(
            Arg::new("max_len")
                .short('L')
                .long("max_len")
                .default_value("4294967295")
                .value_parser(value_parser!(u32).range(1..=4294967295))
                // .value_parser(|x: &str| positive_int_parse(x, "--max_len",  1, u32::MAX as usize))
                .help("max read length")
        )
        .arg(
            Arg::new("min_qual")
                .short('q')
                .long("min_qual")
                .default_value("7.0")
                .value_parser(|x: &str| positive_f64_parse(x, "--min_qual", 0.0, 50.0f64))
                .help("min read quality")
        )
        .arg(
            Arg::new("max_qual")
                .short('Q')
                .long("max_qual")
                .default_value("50.0")
                .value_parser(|x: &str| positive_f64_parse(x, "--min_qual", 0.0, 50.0f64))
                .help("max read quality, usually, you don't need to change this")
        )
        .arg(
            Arg::new("use_dorado_q")
                .short('u')
                .long("use_dorado_q")
                .action(ArgAction::SetTrue)
                .help("use Dorado Q-score calculation: trim leading 60 bases if read length > 60 before calculate read quality")
        )
        .arg(
            Arg::new("gc")
                .long("gc")
                .action(ArgAction::SetTrue)
                .help("whether gc content is used to filter read")
        )
        .arg(
            Arg::new("min_gc")
                .short('g')
                .long("min_gc")
                .default_value("0.0")
                .value_parser(value_parser!(f64))
                .value_parser(|x: &str| positive_f64_parse(x, "--min_qual", 0.0, 1.0f64))
                .help("min gc content when --gc is set")
        )
        .arg(
            Arg::new("max_gc")
                .short('G')
                .long("max_gc")
                .default_value("1.0")
                .value_parser(value_parser!(f64))
                .value_parser(|x: &str| positive_f64_parse(x, "--min_qual", 0.0, 1.0f64))
                .help("max gc content when --gc is set")
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
            Arg::new("retain_failed")
                .long("retain_failed")
                .help("whether to save the failed records, if set, it should be path of failed fastq")
        )
}
