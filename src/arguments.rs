use crate::trim::adapter::get_seq_info;
use ansi_term::Color;
use clap::{Arg, ArgAction, ArgGroup, ArgMatches, Command, value_parser};
use std::fmt::Display;
use std::path::Path;
use std::str::FromStr;

const U32_MAX: &str = "4294967295";
fn input_value_parser(input: &str) {
    let input_path = Path::new(input);
    match input_path.try_exists() {
        Ok(ok) => {
            if !ok {
                eprintln!(
                    "{}",
                    Color::Red.paint(format!(
                        "{}: No such file or directory, check --input",
                        input
                    ))
                );
                std::process::exit(1);
            } else {
                if input_path.is_file() {
                    if !(input.ends_with(".fastq")
                        || input.ends_with(".fq")
                        || input.ends_with(".fastq.gz")
                        || input.ends_with(".fq.gz"))
                    {
                        eprintln!("{}", Color::Red.paint("bad suffix for input file, possible suffix is one of [.fastq, .fq, .fastq.gz, .fq.gz], check --input"));
                        std::process::exit(1);
                    }
                } else if input_path.is_dir() {
                    let mut count = 0;
                    for entry in input_path
                        .read_dir()
                        .expect(&format!("open directory: {:?} failed", input_path))
                    {
                        if let Ok(entry) = entry {
                            let p = entry.path();
                            let p = p.to_str().unwrap();
                            if p.ends_with(".fastq")
                                || p.ends_with(".fq")
                                || p.ends_with(".fastq.gz")
                                || p.ends_with(".fq.gz")
                            {
                                count += 1;
                            }
                        }
                    }
                    if count == 0 {
                        eprintln!(
                            "{}",
                            Color::Red.paint(format!(
                                "No fastq or fastq.gz file found in directory: {}",
                                input
                            ))
                        );
                        std::process::exit(1);
                    }
                } else {
                    eprintln!("{}", Color::Red.paint("input should be file or directory"));
                    std::process::exit(1);
                }
            }
        }
        Err(error) => {
            eprintln!("{}", Color::Red.paint(format!("{:?}", error)));
            std::process::exit(1)
        }
    }
}

fn positive_number_parse<T: FromStr + PartialOrd + Display>(
    x: &str,
    para: &str,
    float: bool,
    min: T,
    max: T,
) -> Result<T, anyhow::Error> {
    let min_length = match x.parse::<T>() {
        Ok(value) => {
            if value < min || value > max {
                eprintln!(
                    "{}",
                    Color::Red.paint(format!(
                        "Error: {} must be between {} and {}",
                        para, min, max
                    ))
                );
                std::process::exit(1);
            }
            value
        }
        Err(err) => {
            let num_type = if float { "float" } else { "int" };
            eprintln!(
                "{}",
                Color::Red.paint(format!("Error: {} must be positive {}", para, num_type))
            );
            std::process::exit(1);
        }
    };
    Ok(min_length)
}

pub fn parse_arguments() -> ArgMatches {
    let input_arg = Arg::new("input")
        .short('i')
        .long("input")
        .action(ArgAction::Set)
        .value_parser(|input:&str| {
            input_value_parser(&input);
            Result::<String, anyhow::Error>::Ok(input.to_string())
        })
        .help("the input fastq, may be a single fastq[.gz] or a directory containing some fastq[.gz], default stdin");

    let output_arg = Arg::new("output")
        .short('o')
        .long("output")
        .action(ArgAction::Set);

    let thread_arg = Arg::new("thread")
        .short('t')
        .long("thread")
        .action(ArgAction::Set)
        .default_value("1")
        .value_parser(value_parser!(u16).range(1..=32))
        .help("how many threads will be used");

    let stats_cmd = Command::new("stats")
                .about("stats fastq")
                .arg( &input_arg )
                .arg(output_arg.clone().help("output the stats result into this, a tsv file or default stdout. it will be truncated if it's a existing file"))
                .arg(Arg::new("summary")
                        .short('s')
                        .long("summary")
                        .action(ArgAction::Set)
                        .default_value("./NanofqStatsSummary.txt")
                        .help("output stats summary into this file, it will be truncated if it exists"))
                .arg(Arg::new("topn")
                        .short('n')
                        .long("topn")
                        .action(ArgAction::Set)
                        .default_value("5")
                        .value_parser(value_parser!(u16))
                        .help("write the top N longest reads and highest quality reads info into summary file"))
                .arg(Arg::new("quality")
                        .short('q')
                        .long("quality")
                        .value_parser(|x: &str| {
                            let mut qualities = x.split(",")
                                .into_iter()
                                .map(|each| {
                                    match each.parse::<f64>() {
                                        Ok(qual) => qual,
                                        Err(_) => {
                                            eprintln!("{}", Color::Red.paint("Error: parse f64 from each quality, check --quality"));
                                            std::process::exit(1);
                                        }
                                    }
                                })
                                .collect::<Vec<f64>>();
                            // decrease quality
                            qualities.sort_by(|a, b| b.partial_cmp(a).unwrap());
                            Result::<Vec<f64>, anyhow::Error>::Ok(qualities)
                        })
                        .default_value("25,20,18,15,12,10")
                        .help("count the reads number that whose quality is bigger than this value, multi value can be separated by comma"))
                .arg(Arg::new("length")
                        .short('l')
                        .long("length")
                        .value_parser(|x: &str| {
                            let mut lengths = x.split(",")
                                .into_iter()
                                .map(|each| {
                                    match each.parse::<usize>() {
                                        Ok(len) => len,
                                        Err(err) => {
                                            eprintln!("{:?}", err);
                                            eprintln!("{}", Color::Red.paint("Error: parse usize from each length, check --length"));
                                            std::process::exit(1)
                                        }
                                    }
                                })
                                .collect::<Vec<usize>>();
                            lengths.sort_by(|a,b| b.partial_cmp(a).unwrap());
                            Result::<Vec<usize>, anyhow::Error>::Ok(lengths)
                        })
                        .help("count the reads number that whose length is bigger than this value, multi values can be separated by comma"))
                .arg(Arg::new("gc")
                        // .short('g')
                        .long("gc")
                        .action(ArgAction::SetTrue)
                        .help("whether stats the gc content"))
                .arg(thread_arg.clone())
                .arg(Arg::new("plot")
                        .short('p')
                        .long("plot")
                        .action(ArgAction::SetTrue)
                        .help("whether make plot"))
                .arg(Arg::new("format")
                        .short('f')
                        .long("format")
                        .action(ArgAction::Append)
                        .value_parser(["png", "pdf", "jpg"])
                        .default_value("png")
                        .help("which format figure do you want if --plot is true, this para can be set multi times"));
    let filter_cmd = Command::new("filter")
        .about("filter fastq")
        .arg(&input_arg)
        .arg(
            output_arg
                .clone()
                .help("output the filtered fastq into this file or default stdout, it will be truncated if it's a existing file"),
        )
        .arg(
            Arg::new("min_len")
                .short('l')
                .long("min_len")
                .action(ArgAction::Set)
                .default_value("1")
                .value_parser(|x: &str| positive_number_parse::<usize>(x, "--min_len", false, 1, u32::MAX as usize))
                .help("min read length"),
        )
        .arg(
            Arg::new("max_len")
                .short('L')
                .long("max_len")
                .action(ArgAction::Set)
                .default_value(U32_MAX)
                .value_parser(|x: &str| positive_number_parse::<usize>(x, "--max_len", false, 1, u32::MAX as usize))
                .help("min read length"),
        )
        .arg(
            Arg::new("min_qual")
                .short('q')
                .long("min_qual")
                .action(ArgAction::Set)
                .default_value("1.0")
                .value_parser(|x: &str| positive_number_parse::<f64>(x, "--min_qual", true, 1.0f64, 50.0f64))
                .help("max read qual"),
        )
        .arg(
            Arg::new("max_qual")
                .short('Q')
                .long("max_qual")
                .action(ArgAction::Set)
                .default_value("50.0")
                .value_parser(|x: &str| positive_number_parse::<f64>(x, "--max_qual", true, 1.0f64, 50.0f64))
                .help("max read qual, but in most cases, you won't specify this value"),
        )
        .arg(
            Arg::new("gc")
                // .short('g')
                .long("gc")
                .action(ArgAction::SetTrue)
                .help("whether use gc content to filter read [default: false]"),
        )
        .arg(
            Arg::new("min_gc")
                .short('g')
                .long("min_gc")
                .default_value("0.0")
                .action(ArgAction::Set)
                .value_parser(|x: &str| positive_number_parse::<f64>(x, "--min_gc", true, 0.0f64, 1.0f64))
                .help("min gc content if --gc is set true"),
        )
        .arg(
            Arg::new("max_gc")
                .short('G')
                .long("max_gc")
                .default_value("1.0")
                .action(ArgAction::Set)
                .value_parser(|x: &str| positive_number_parse::<f64>(x, "--max_gc", true, 0.0f64, 1.0f64))
                .help("max gc content if --gc is set true"),
        )
        .arg(thread_arg.clone())
        .arg(
            Arg::new("retain_failed")
                .long("retain_failed")
                .action(ArgAction::Set)
                .help("whether store the failed fastq, if set, this value should be the path of failed fastq. this file will be truncated if it exists")
        );
    let trim_cmd = Command::new("trim")
        .about("trim adapter, barcode, primer that artificial sequence from long fastq")
        .arg(&input_arg)
        .arg(output_arg
            .clone()
            .help("output the trimmed fastq into this file or default stdout, it will be truncated if it's a existing file"),
        )
        .arg(
            Arg::new("log")
                .short('l')
                .long("log")
                .action(ArgAction::Set)
                .help("whether store the trimmed log, if set, this value should be the path of trimmed log file, this file will be truncated if it exists")
        )
        .arg(
            Arg::new("kit")
                .short('k')
                .long("kit")
                .help("Which kit you used. Each kit has it's own search parameter, but can be changed by [search parameter]. NBD_{number} means kit name with barcode number.")
                .default_value("LSK")
                .value_parser(|x: &str| {
                    get_seq_info().contains_key(x);
                    if get_seq_info().contains_key(x) {
                        return Result::<String, anyhow::Error>::Ok(x.to_string())
                    } else {
                        eprintln!("Error: invalid \'{}\' for \'--kit <kit>\'", x);
                        eprintln!("[possible values: LSK, RAD, ULK, RBK, PCS, PCB, NBD_1, NBD_2, ..., NBD_95, NBD_96]");
                        std::process::exit(1);
                    }
                })
        )
        .arg(
            Arg::new("primers")
                .short('p')
                .long("primers")
                .value_parser(|x: &str| {
                    let bases = [b'A', b'T', b'C', b'G', b'R', b'Y', b'M', b'K', b'S', b'W', b'H', b'B', b'V', b'D', b','];
                    x.as_bytes().iter().for_each(|base| {
                       if !bases.contains(base) {
                           eprintln!("Error: Invalid base char found in primers");
                           std::process::exit(1);
                       }
                    });
                    let comma_numbers = x.as_bytes().iter().filter(|base| *base == &b',').count();
                    if comma_numbers != 1 || x.as_bytes().last().unwrap() == &b','{
                        eprintln!("One and only one comma must be used to separate paired primers");
                        std::process::exit(1);
                    }
                    Result::<String, anyhow::Error>::Ok(x.to_string())
                })
                .help("a paired primers separated by comma, the first one is forward primer and second is reversed, the direction should from 5' end to 3' end. Degenerate bases supported")
        )
        .group(
            ArgGroup::new("seq")
                .arg("kit")
                .arg("primers")
                .required(true)
                .multiple(false)
        )
        .arg( thread_arg.clone() )
        .arg(
            Arg::new("match")
                .short('m')
                .long("match")
                .default_value("3")
                .value_parser(value_parser!(i32).range(1..=100))
                .help("match score, positive int")
        )
        .arg(
            Arg::new("mismatch")
                .short('M')
                .long("mismatch")
                .default_value("-3")
                .value_parser(value_parser!(i32).range(-100..=0))
                .help("mismatch penalty score, negative int")
        )
        .arg(
            Arg::new("gap_opened")
                .short('g')
                .long("gap_open")
                .default_value("-5")
                .value_parser(value_parser!(i32).range(-100..=0))
                .help("gap opened penalty score, negative int")
        )
        .arg(
            Arg::new("gap_extend")
                .short('G')
                .long("gap_extend")
                .default_value("-1")
                .value_parser(value_parser!(i32).range(-100..=0))
                .help("gap extend penalty score, negative int")
        )
        .arg(
            Arg::new("rev_com_not_used")
                .long("rev_com_not_used")
                .action(ArgAction::SetTrue)
                .help("whether used rev com sequences of primers to query in read if primers is used. If it's false, \
                we will assume that fwd primer is in 5'end of read and rev_com of rev primer is in 3'end of read")
        )
        .arg(
            Arg::new("end5_len")
                .long("end5_len")
                .default_value("150")
                .help("[search parameter]: search in the first N bases from 5'end of reads to find front adapter")
                .value_parser(value_parser!(u32).range(100..=1000))
        )
        .arg(
            Arg::new("end5_align_pct")
                .long("end5_align_pct")
                .help("[search parameter]: the ratio between align length of front adapter and the full length of adapter should be bigger than this value for 5' end")
                .default_value("0.8")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end5_align_ident")
                .long("end5_align_ident")
                .default_value("0.8")
                .help("[search parameter]: the ratio between the identity bases number of align and the align length of adapter should be bigger than this value for 5'end")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end3_len")
                .long("end3_len")
                .default_value("130")
                .help("[search parameter]: search in the last N bases from 3'end of reads to find rear adapter")
                .value_parser(value_parser!(u32).range(100..=1000))
        )
        .arg(
            Arg::new("end3_align_pct")
                .long("end3_align_pct")
                .help("[search parameter]: the ratio between align length of rear adapter and the full length of adapter should be bigger than this value for 3' end")
                .default_value("0.8")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end3_align_ident")
                .long("end3_align_ident")
                .default_value("0.8")
                .help("[search parameter]: the ratio between the identity bases number of align and the align length of adapter should be bigger than this value for 3'end")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end5_len_rc")
                .long("end5_len_rc")
                .default_value("150")
                .help("[search parameter]: search in the first N bases from 5'end of reads to find front adapter if this read is reverse complementary")
                .value_parser(value_parser!(u32).range(100..=1000))
        )
        .arg(
            Arg::new("end5_align_pct_rc")
                .long("end5_align_pct_rc")
                .help("[search parameter]: the ratio between align length of front adapter and the full length of adapter should be bigger than this value for 5' end if this read is reverse complementary")
                .default_value("0.8")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end5_align_ident_rc")
                .long("end5_align_ident_rc")
                .default_value("0.8")
                .help("[search parameter]: the ratio between the identity bases number of align and the align length of adapter should be bigger than this value for 5'end if this read is reverse complementary")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end3_len_rc")
                .long("end3_len_rc")
                .default_value("130")
                .help("[search parameter]: search in the last N bases from 3'end of reads to find rear adapter if this read is reverse complementary")
                .value_parser(value_parser!(u32).range(100..=1000))
        )
        .arg(
            Arg::new("end3_align_pct_rc")
                .long("end3_align_pct_rc")
                .help("[search parameter]: the ratio between align length of rear adapter and the full length of adapter should be bigger than this value for 3' end if this read is reverse complementary")
                .default_value("0.8")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )
        .arg(
            Arg::new("end3_align_ident_rc")
                .long("end3_align_ident_rc")
                .default_value("0.8")
                .help("[search parameter]: the ratio between the identity bases number of align and the align length of adapter should be bigger than this value for 3'end if this read is reverse complementary")
                .value_parser(|x: &str|positive_number_parse::<f64>(x, "--end5_align_pct", true, 0.0f64, 1.0f64))
        )

        ;

    let cmd = Command::new("nanofq")
        .version("0.0.1")
        .about("A tool for nanopore fastq file")
        .subcommand(stats_cmd)
        .subcommand(filter_cmd)
        .subcommand(trim_cmd);
    let x = cmd.get_matches();
    x
}
