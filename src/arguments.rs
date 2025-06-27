use ansi_term::Color;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};

pub fn parse_arguments() -> ArgMatches {
    let input_arg = Arg::new("input")
        .short('i')
        .long("input")
        .action(ArgAction::Set)
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
                            let qualities = x.split(",")
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
                            let lengths = x.split(",")
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
                        .value_parser(value_parser!(u8).range(1..16))
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