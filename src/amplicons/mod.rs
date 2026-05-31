use crate::amplicons::preprocess::{ReadsClassifier, ReadsWithPrimer};
use crate::primer_barcode::{BARCODES, Barcode, Primer};
use crate::utils::{check_and_create_dir, quit_with_error, run_abpoa, run_minimap2_and_index};
use ahash::{HashMap, HashSet, RandomState};
use clap::parser::ValueSource;
use clap::{Arg, ArgMatches, Command, value_parser};
use log::{error, info, warn};
use rayon::prelude::*;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record};
use std::io::{BufWriter, Write};
use std::path::Path;

mod _consensus;
pub mod preprocess;

pub fn amplicon_with_known_primers(
    input_fastq: &str,
    bar_idx: usize,
    left_range: usize,
    right_range: usize,
    max_distance: u8,
    primers: HashMap<String, Primer>,
    reads_downsample: usize,
    output_dir: &str,
    min_read_quality: f64,
    length_range: f64,
    abpoa: Option<&str>,
    analysis_name: &str,
) {
    let barcode = Barcode::new(BARCODES[bar_idx].as_bytes());
    let reads_collector = preprocess::ReadsCollector::new(
        input_fastq,
        &barcode,
        left_range,
        right_range,
        max_distance,
    );
    let all_reads = reads_collector.collect_fastqs();
    let mut classifier = preprocess::ReadsClassifier::new(
        all_reads,
        21,
        left_range,
        right_range,
        max_distance,
        analysis_name.to_string(),
    );
    let primer_name2reads_idx = classifier.classify_reads_with_known_primers(&primers);
    let primer_reads = primer_name2reads_idx
        .into_iter()
        .map(|(primer_name, reads_idx)| {
            let reads = classifier.remove_reads_with_idxes(&reads_idx);
            (
                primer_name.clone(),
                ReadsWithPrimer::new(
                    reads,
                    // primers.get(&primer_name).unwrap().clone(),
                    reads_downsample,
                ),
            )
        })
        .collect::<Vec<_>>();

    check_and_create_dir(output_dir);
    primer_reads
        .into_par_iter()
        .for_each(|(primer_name, mut reads_with_primer)| {
            let save_fastq_path =
                format!("{output_dir}/{analysis_name}_{primer_name}.trimmed_filtered.fastq");
            reads_with_primer.filter(min_read_quality, length_range);
            reads_with_primer.save_fastq(&save_fastq_path);
            run_abpoa(&save_fastq_path, output_dir, &primer_name, abpoa);
        })
}

pub fn amplicon_with_unknown_primers(
    input_fastq: &str,
    bar_idx: usize,
    left_range: usize,
    right_range: usize,
    max_distance: u8,
    amplicons_number: usize,
    lead_length: usize,
    guess_reads_number: usize,
    rev_primer_found_ratio: f64,
    min_rev_match_reads_number: usize,
    reads_downsample: usize,
    output_dir: &str,
    min_read_quality: f64,
    length_range: f64,
    abpoa: Option<&str>,
    minimap2: Option<&str>,
    samtools: Option<&str>,
    min_mapq: u8,
    analysis_name: &str,
) {
    let mut detected_primers = vec![];
    let barcode = Barcode::new(BARCODES[bar_idx].as_bytes());
    let reads_collector = preprocess::ReadsCollector::new(
        input_fastq,
        &barcode,
        left_range,
        right_range,
        max_distance,
    );
    let all_reads = reads_collector.collect_fastqs();
    let mut classifier = preprocess::ReadsClassifier::new(
        all_reads,
        lead_length,
        left_range,
        right_range,
        max_distance,
        analysis_name.to_string(),
    );
    for primer_idx in 1..amplicons_number + 1 {
        info!("Analysis {primer_idx}th Amplicon......");
        // Step1.1: guess one paired primer from barcode trimmed fastq
        let primer_opt = classifier.detect_one_primer(
            guess_reads_number,
            min_rev_match_reads_number,
            rev_primer_found_ratio,
            primer_idx,
            false,
            10,
        );
        let primer = match primer_opt {
            None => {
                // Step1.2: if no primer detected, maybe because too few lead seq freq,
                // detect primer again with similar lead seqs merged
                let primer_opt2 = classifier.detect_one_primer(
                    guess_reads_number,
                    min_rev_match_reads_number,
                    rev_primer_found_ratio,
                    primer_idx,
                    true,
                    10,
                );
                warn!("No primer detect, try to merge similar lead seqs of reads");
                match primer_opt2 {
                    None => {
                        error!("No primer detect after similar lead seqs merged");
                        quit_with_error(&format!("Couldn't detect {primer_idx}th primer"))
                    }
                    Some(primer) => {
                        info!(
                            "Primers named by [{}] found with Forward: {} and Reverse: {}",
                            primer.name,
                            primer.fwd(),
                            primer.rev()
                        );
                        primer
                    }
                }
            }
            Some(primer) => {
                info!(
                    "Primers named by [{}] found with Forward: {} and Reverse: {}",
                    primer.name,
                    primer.fwd(),
                    primer.rev()
                );
                primer
            }
        };
        detected_primers.push(primer.clone());
        classifier = draft_consensus_with_one_known_primer(
            primer,
            classifier,
            min_read_quality,
            length_range,
            reads_downsample,
            output_dir,
            abpoa,
            minimap2,
            samtools,
            min_mapq,
        );
    }
    let mut output_primers = BufWriter::new(
        std::fs::File::create(format!("{output_dir}/{analysis_name}_detected_primers.txt"))
            .unwrap(),
    );
    for primer in &detected_primers {
        writeln!(
            &mut output_primers,
            "{}\t{}\t{}",
            primer.name,
            primer.fwd(),
            primer.rev()
        )
        .unwrap()
    }
}

fn draft_consensus_with_one_known_primer(
    primer: Primer,
    classifier: ReadsClassifier,
    min_read_quality: f64,
    length_range: f64,
    reads_downsample: usize,
    output_dir: &str,
    abpoa: Option<&str>,
    minimap2: Option<&str>,
    samtools: Option<&str>,
    min_mapq: u8,
) -> ReadsClassifier {
    let primers = primer.name2primer();
    let mut classifier = classifier;

    // Step2: Trim sequence outside this known primer of amplicon,
    // try to get clean reads indexes that startswith fwd primer and ends with rev primer rc
    let primer_name2reads_idx = classifier.classify_reads_with_known_primers(&primers);

    for (primer_name, reads_idx) in primer_name2reads_idx.iter() {
        info!(
            "{} reads with paired primers [{}] found",
            reads_idx.len(),
            primer_name,
        );
    }

    // assert_eq!(primer_name2reads_idx.len(), 1);
    let (primer_name, reads_idx) = primer_name2reads_idx
        .into_iter()
        .collect::<Vec<_>>()
        .pop()
        .unwrap();

    // Step3: Remove this clean reads from classifier and save them in ReadsWithPrimer
    let reads = classifier.remove_reads_with_idxes(&reads_idx);
    info!(
        "{} reads with paired primer [{}] used to construct ReadsWithPrimer",
        reads.len(),
        primer_name
    );
    let mut reads_with_primer = ReadsWithPrimer::new(reads, reads_downsample);

    // Step4: Filter reads depending on read quality and length
    let bad_records_size = reads_with_primer.filter(min_read_quality, length_range);
    info!("{bad_records_size} reads dropped after filter depending on read quality and length");

    // Step5: Save this clean fastq and remaining barcode trimmed fastq and redundant records
    check_and_create_dir(output_dir);
    let redundant_fastq = format!("{output_dir}/{primer_name}.redundant.trimmed_filtered.fastq");
    reads_with_primer.save_redundant_fastq(&redundant_fastq);
    info!(
        "{} redundant reads hadn't be used to run abpoa found, saved into {redundant_fastq}",
        reads_with_primer.redundant_reads.len()
    );

    let this_primer_fastq = format!("{output_dir}/{primer_name}.trimmed_filtered.fastq",);
    reads_with_primer.save_fastq(&this_primer_fastq);
    info!(
        "{} reads used to run abpoa and construct draft consensus, saved into {this_primer_fastq}",
        reads_with_primer.reads.len()
    );

    let barcode_trimmed_remaining_fastq = format!( "{output_dir}/{primer_name}_remaining.fastq" );
    classifier.save_fastq(&barcode_trimmed_remaining_fastq);
    info!("{primer_name} remaining records saved into {barcode_trimmed_remaining_fastq}");

    // Step6: Run abpoa to construct draft consensus
    info!("Running abpoa for {primer_name}...");
    let draft_consensus = run_abpoa(&this_primer_fastq, output_dir, &primer_name, abpoa);

    /*
    Step7 and Step8 are used to remove reads that mapped to this amplicon draft consensus
    avoiding them to disturb next classifier.guess_one_primer.
    */
    // Step7: Map all remaining barcod trimmed fastq records to this draft consensus.
    info!(
        "Map remaining reads to draft_consensus to search this amplicon's reads that has no paired primers detected on read"
    );
    let sorted_bam = run_minimap2_and_index(
        output_dir,
        &barcode_trimmed_remaining_fastq,
        &draft_consensus,
        &primer_name,
        minimap2,
        samtools,
    );

    // Step8: Remove mapped reads from classifier's all_reads
    let mapped_read_names = get_reads_name_from_bam(&sorted_bam, min_mapq);
    classifier.remove_reads_with_names(&mapped_read_names);
    info!(
        "{} reads found and remove them from ReadsClassifier and go to find next amplicon",
        mapped_read_names.len()
    );
    /*
    Theoretically, all reads from this amplicon have been removed from classifier here.
    Some were removed by classifier.remove_reads_with_indexes. These reads used to construct draft_consensus
    Some were removed by classifier.remove_reads_with_names. These reads just be removed.
     */

    // Step9: For next loop
    classifier
}

fn get_reads_name_from_bam(bam_file: &str, min_mapq: u8) -> HashSet<String> {
    let mut read_names = HashSet::with_hasher(RandomState::new());
    let mut bam_reader =
        IndexedReader::from_path(bam_file).expect(&format!("Failed to read {bam_file}"));
    bam_reader.fetch(FetchDefinition::All).unwrap();
    let mut bam_record = Record::new();
    while let Some(Ok(_)) = bam_reader.read(&mut bam_record) {
        if bam_record.mapq() >= min_mapq {
            let read_name = str::from_utf8(bam_record.qname()).unwrap();
            read_names.insert(read_name.to_string());
        }
    }
    read_names
}

fn parse_primers_from_cli(primers: &str, analysis_name: &str) -> HashMap<String, Primer> {
    let pattern =
        r"^([A-Za-z0-9_-]+,[ATCGatcg]+,[ATCGatcg]+)(?:;([A-Za-z0-9_-]+,[ATCGatcg]+,[ATCGatcg]+))*$";
    let pattern = regex::Regex::new(pattern).unwrap();
    if pattern.is_match(primers) {
        Primer::parse_primers_from_str(primers, analysis_name)
    } else if Path::new(primers).is_file() {
        Primer::parse_primer_from_file(primers, analysis_name)
    } else {
        quit_with_error("Bad --primers format found")
    }
}

pub fn run_amplicons(amp_cmd: &ArgMatches) {
    let input = amp_cmd.get_one::<String>("input").unwrap();
    let output = amp_cmd.get_one::<String>("output").unwrap();
    let barcode = amp_cmd.get_one::<u32>("barcode").unwrap();
    let left = amp_cmd.get_one::<usize>("left").unwrap();
    let right = amp_cmd.get_one::<usize>("right").unwrap();
    let distance = amp_cmd.get_one::<u8>("distance").unwrap();
    let primers_opt = amp_cmd.get_one::<String>("primers");
    let downsample = amp_cmd.get_one::<usize>("downsample").unwrap();
    let min_qual = amp_cmd.get_one::<f64>("min_qual").unwrap();
    let len_range = amp_cmd.get_one::<f64>("len_range").unwrap();
    let amp_numbers = amp_cmd.get_one::<usize>("number").unwrap();
    let lead_length = amp_cmd.get_one::<usize>("lead").unwrap();
    let detect_rev_primer_reads_number =
        amp_cmd.get_one::<usize>("detect_rev_primer_reads").unwrap();
    let rev_primer_found_ratio = amp_cmd.get_one::<f64>("rev_primer_found_ratio").unwrap();
    let min_rev_match_reads_number = amp_cmd.get_one::<usize>("min_rev_match_reads").unwrap();
    let min_mapq = amp_cmd.get_one::<u8>("min_mapq").unwrap();
    let abpoa = amp_cmd.get_one::<String>("abpoa").map(|x| x.as_ref());
    let minimap2 = amp_cmd.get_one::<String>("minimap2").map(|x| x.as_ref());
    let samtools = amp_cmd.get_one::<String>("samtools").map(|x| x.as_ref());
    let analysis_name = amp_cmd.get_one::<String>("prefix").unwrap();

    let primers_is_set = amp_cmd.value_source("primers") == Some(ValueSource::CommandLine);
    let amplicon_number_is_set = amp_cmd.value_source("number") == Some(ValueSource::CommandLine);
    if primers_is_set && amplicon_number_is_set {
        quit_with_error("--primers and number couldn't be specified together")
    }

    // set output absolute path
    let output_path = std::env::current_dir()
        .expect("Failed to get current directory path")
        .join(Path::new(output).to_owned());
    let output = output_path
        .to_str()
        .expect("Failed to set output directory");

    match primers_opt {
        None => amplicon_with_unknown_primers(
            input,
            *barcode as usize,
            *left,
            *right,
            *distance,
            *amp_numbers,
            *lead_length,
            *detect_rev_primer_reads_number,
            *rev_primer_found_ratio,
            *min_rev_match_reads_number,
            *downsample,
            output,
            *min_qual,
            *len_range,
            abpoa,
            minimap2,
            samtools,
            *min_mapq,
            analysis_name,
        ),
        Some(primers) => {
            let primers = parse_primers_from_cli(primers, analysis_name);
            amplicon_with_known_primers(
                input,
                *barcode as usize,
                *left,
                *right,
                *distance,
                primers,
                *downsample,
                output,
                *min_qual,
                *len_range,
                abpoa,
                analysis_name,
            )
        }
    }
}

pub fn amplicons_cmd() -> Command {
    Command::new("amplicon")
        .about("get draft consensuses from mixed Ligation Nanopore Long amplicons reads with known or unknown primers")
        .arg(
            Arg::new("input")
            .short('i')
            .long("input")
            .required(true)
            .help("the input fastq[.gz] file")
    ).arg(
        Arg::new("output")
            .short('o')
            .long("output")
            .required(true)
            .help("output directory")
    ).arg(
        Arg::new("primers")
            .short('p')
            .long("primers")
            .help("known primers string with format: <PrimerName>,<FwdPrimer>,<RevPrimer>[;PrimerName,FwdPrimer,RevPrimer...] or primers file with each line format: PrimerName\\tFwdPrimer\\tRevPrimer")
    ).arg(
        Arg::new("number")
            .short('n')
            .long("number")
            .default_value("1")
            .value_parser(value_parser!(usize))
            .help("when no known primers specified, how many amplicons mixed in this sample, used with [UnknownPrimers] parameters")
    ).arg(
        Arg::new("barcode")
            .short('b')
            .long("barcode")
            .default_value("0")
            .value_parser(value_parser!(u32).range(0..=96))
            .help("which barcode used, chose from 0-96. 0 means Ligation Sequence Kit used without barcode")
    ).arg(
        Arg::new("left")
            .short('l')
            .long("left")
            .required(false)
            .default_value("150")
            .value_parser(value_parser!(usize))
            .help("the first this value bases of read used to detect barcode or primer")
    ).arg(
        Arg::new("right")
            .short('r')
            .long("right")
            .required(false)
            .default_value("150")
            .value_parser(value_parser!(usize))
            .help("the last this value bases of read used to detect barcode or primer")
    ).arg(
        Arg::new("distance")
            .short('d')
            .long("distance")
            .default_value("3")
            .value_parser(value_parser!(u8))
            .help("min distance between barcode or primer with sequence in read")
    ).arg(
        Arg::new("downsample")
            .long("downsample")
            .default_value("10000")
            .value_parser(value_parser!(usize))
            .help("max reads number that with paired primer at dual reads will be used to construct consensus")
    ).arg(
        Arg::new("min_qual")
            .long("min_qual")
            .default_value("15")
            .value_parser(value_parser!(f64))
            .help("min read quality that with paired primers at dual reads")
    ).arg(
        Arg::new("len_range")
            .long("len_range")
            .default_value("0.05")
            .value_parser(value_parser!(f64))
            .help("length of reads that with paired primer should be in mean_length * (1-len_range) and mean_length * (1+len_range)")
    ).arg(
        Arg::new("lead")
            .long("lead")
            .default_value("21")
            .value_parser(value_parser!(usize))
            .help("when no known primers specified, use first lead bases as candidate fwd primer after barcode trimmed")
    ).arg(
        Arg::new("detect_rev_primer_reads")
            .long("detect_rev_primer_reads")
            .default_value("1000")
            .value_parser(value_parser!(usize))
            .help("[UnknownPrimers]: how many reads used to detect rev primer for reads with candidate fwd primer")
    ).arg(
        Arg::new("rev_primer_found_ratio")
            .long("rev_primer_found_ratio")
            .default_value("0.5")
            .value_parser(value_parser!(f64))
            .help("[UnknownPrimers]: if rev primer can be detected in rear of this proportion reads, consider a paired primers found")
    ).arg(
        Arg::new("min_rev_match_reads")
            .long("min_rev_match_reads")
            .default_value("100")
            .value_parser(value_parser!(usize))
            .help("[UnknownPrimers]: if reads number with rev primer detected in rear is less than this value, ignore it")
    ).arg(
        Arg::new("min_mapq")
            .long("min_mapq")
            .default_value("50")
            .value_parser(value_parser!(u8))
            .help("[UnknownPrimers]: the min mapq used to collect reads that with no paired primer detected but can be mapped to draft consensus")
    ).arg(
        Arg::new("abpoa")
            .long("abpoa")
            .default_value("abpoa")
            .help("abpoa path")
    ).arg(
        Arg::new("minimap2")
            .long("minimap2")
            .default_value("minimap2")
            .help("[UnknownPrimers]: minimap2 path")
    ).arg(
        Arg::new("samtools")
            .long("samtools")
            .default_value("samtools")
            .help("[UnknownPrimers]: samtools path")
    ).arg(
        Arg::new("prefix")
            .long("prefix")
            .default_value("test001")
            .help("the prefix of output file name")
    )
}
