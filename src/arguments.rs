use ansi_term::Color;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use std::path::Path;

fn input_value_parser(input: &str) {
    let input_path = Path::new(input);
    if !input_path.is_absolute() {
        eprintln!("{}", Color::Red.paint("input must be absolute path"));
        std::process::exit(1);
    }
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

    let cmd = Command::new("nanofq")
        .version("0.0.1")
        .about("A tool for nanopore fastq file")
        .subcommand(
            Command::new("stats")
                .about("stats fastq")
                .arg( &input_arg )
                .arg(output_arg
                         .help("output the stats result into this, a tsv file or default stdout")
                )
                .arg(
                    Arg::new("summary")
                        .short('s')
                        .long("summary")
                        .action(ArgAction::Set)
                        .default_value("./NanofqStatsSummary.txt")
                        .help("output stats summary into this file")
                )
                .arg(
                    Arg::new("topn")
                        .short('n')
                        .long("topn")
                        .action(ArgAction::Set)
                        .default_value("5")
                        .value_parser(value_parser!(u16))
                        .help("write the top N longest reads and highest quality reads info into summary file")
                )
                .arg(
                    Arg::new("quality")
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
                        .help("count the reads number that whose quality is bigger than this value, multi value can be separated by coma")
                )
                .arg(
                    Arg::new("length")
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
                        .help("count the reads number that whose length is bigger than this value, multi values can be separated by coma")
                )
                .arg(
                    Arg::new("gc")
                        .short('g')
                        .long("gc")
                        .action(ArgAction::SetTrue)
                        .help("whether stats the gc content")
                )
                .arg(
                    Arg::new("thread")
                        .short('t')
                        .long("thread")
                        .action(ArgAction::Set)
                        .default_value("1")
                        .value_parser(value_parser!(u16).range(1..16))
                        .help("how many threads used to stats fastqs")
                )
                .arg(
                    Arg::new("plot")
                        .short('p')
                        .long("plot")
                        .action(ArgAction::SetTrue)
                        .help("whether make plot")
                )
                .arg(
                    Arg::new("format")
                        .short('f')
                        .long("format")
                        .action(ArgAction::Append)
                        .value_parser(["png", "pdf", "jpg"])
                        .default_value("png")
                        .help("which format figure do you want if --plot is true, this para can be set multi times")
                )
        );
    cmd.get_matches()
}
