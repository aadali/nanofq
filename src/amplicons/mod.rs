use crate::amplicons::preprocess::{ReadsClassifier, ReadsWithPairedPrimers};
use crate::primer_barcode::{BARCODES, Barcode, Primer};
use crate::utils::{
    check_and_create_dir, init_log, quit_with_error, run_abpoa, run_minimap2_and_index,
};
use ahash::{HashMap, HashSet, RandomState};
use clap::parser::ValueSource;
use clap::{Arg, ArgAction, ArgMatches, Command, value_parser};
use log::info;
use rayon::prelude::*;
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record};
use std::io::BufWriter;
use std::path::Path;

mod _consensus;
pub mod preprocess;

struct FileNameSuffix {
    clean: &'static str,
    paired_primers_reads: &'static str,
    good_reads: &'static str,
    redundant_reads: &'static str,
    bad_reads: &'static str,
    detected_primers: &'static str,
    remaining_reads: &'static str,
}

impl Default for FileNameSuffix {
    fn default() -> Self {
        FileNameSuffix {
            clean: ".clean.fastq",
            paired_primers_reads: ".with_paired_primers.fastq",
            good_reads: ".with_paired_primers.good.fastq",
            redundant_reads: ".with_paired_primers.redundant.fastq",
            bad_reads: ".with_paired_primers.bad.fastq",
            detected_primers: ".detected_primers.tsv",
            remaining_reads: ".remaining.fastq",
        }
    }
}

fn amplicon_with_known_primers(
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
    save_failed: bool,
    file_name_suffix: &FileNameSuffix,
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
    let mut classifier = ReadsClassifier::new(
        all_reads,
        21,
        left_range,
        right_range,
        max_distance,
        analysis_name.to_string(),
    );
    check_and_create_dir(output_dir);
    let clean_fastq = format!("{output_dir}/{analysis_name}{}", file_name_suffix.clean);
    classifier.save_clean_fastq(&clean_fastq);
    let primer_name2reads_idx = classifier.classify_reads_with_known_primers(&primers);
    primer_name2reads_idx
        .par_iter()
        .for_each(|(primer_name, reads_idx)| {
            let fq = format!(
                "{output_dir}/{analysis_name}_{primer_name}{}",
                file_name_suffix.paired_primers_reads
            );
            let mut file = BufWriter::new(std::fs::File::create(fq).unwrap());
            for read_idx in reads_idx {
                let record = classifier.all_reads.get(read_idx).unwrap();
                record.write(&mut file).unwrap();
            }
            info!(
                "{primer_name}: {} reads with paired primers at dual ends found",
                reads_idx.len()
            );
        });

    let primer_reads = primer_name2reads_idx
        .into_iter()
        .map(|(primer_name, reads_idx)| {
            let reads = classifier.remove_reads_with_idxes(&reads_idx);
            (
                primer_name.clone(),
                ReadsWithPairedPrimers::new(
                    reads,
                    reads_downsample,
                ),
            )
        })
        .collect::<Vec<_>>();

    check_and_create_dir(output_dir);
    primer_reads
        .into_par_iter()
        .for_each(|(primer_name, mut reads_with_primer)| {
            let total_reads_with_paired_primers = reads_with_primer.reads.len();
            let this_primer_fastq = format!(
                "{output_dir}/{analysis_name}_{primer_name}{}",
                file_name_suffix.good_reads
            );
            let failed_fastq = format!(
                "{output_dir}/{analysis_name}_{primer_name}{}",
                file_name_suffix.bad_reads
            );
            let bad_reads_size = reads_with_primer.filter(
                min_read_quality,
                length_range,
                if save_failed {
                    Some(&failed_fastq)
                } else {
                    None
                },
            );

            info!(
                "{primer_name}: {bad_reads_size}/{total_reads_with_paired_primers} bad reads dropped{}",
                if save_failed {
                    format!(", saved into {}", failed_fastq)
                } else {
                    "".to_string()
                }
            );

            if !reads_with_primer.redundant_reads.is_empty() {
                let redundant_fastq = format!(
                    "{output_dir}/{analysis_name}_{primer_name}{}",
                    file_name_suffix.redundant_reads
                );
                reads_with_primer.save_redundant_fastq(&redundant_fastq);

                info!(
                    "{primer_name}: {}/{total_reads_with_paired_primers} redundant reads found, saved into {}",
                    reads_with_primer.redundant_reads.len(),
                    &redundant_fastq
                )
            }
            reads_with_primer.save_fastq(&this_primer_fastq);
            info!(
                "{primer_name}: {}/{total_reads_with_paired_primers} good reads used as input of abpoa, saved into {this_primer_fastq}",
                reads_with_primer.reads.len()
            );
            run_abpoa(&this_primer_fastq, output_dir, &primer_name, abpoa);
        })
}

fn amplicon_with_unknown_primers(
    input_fastq: &str,
    bar_idx: usize,
    left_range: usize,
    right_range: usize,
    max_distance: u8,
    amplicons_number: usize,
    lead_length: usize,
    guess_reads_number: usize,
    reads_downsample: usize,
    output_dir: &str,
    min_read_quality: f64,
    length_range: f64,
    abpoa: Option<&str>,
    minimap2: Option<&str>,
    samtools: Option<&str>,
    min_mapq: u8,
    analysis_name: &str,
    save_failed: bool,
    file_name_suffix: &FileNameSuffix,
    thread: usize,
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
    let mut classifier = ReadsClassifier::new(
        all_reads,
        lead_length,
        left_range,
        right_range,
        max_distance,
        analysis_name.to_string(),
    );
    check_and_create_dir(output_dir);
    classifier.save_clean_fastq(&format!(
        "{output_dir}/{analysis_name}{}",
        file_name_suffix.clean
    ));
    for primer_idx in 1..amplicons_number + 1 {
        {
            info!("Analysis {primer_idx}th Amplicon......");
        }
        // Step1.1: guess one paired primer from barcode trimmed fastq
        let detected_paired_primer =
            classifier.detect_one_primer(guess_reads_number, amplicons_number, true, 5);

        let primer =
            detected_paired_primer.generate_primer(&format!("{analysis_name}_primer{primer_idx}"));

        {
            info!(
                "Primer named by {} detected: Fwd: {} >> {} and Rev: {} >> {}",
                primer.name,
                primer.fwd(),
                detected_paired_primer.fwd_primer_reads.len(),
                primer.rev(),
                detected_paired_primer.rev_primer_reads.len(),
            );
        }

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
            save_failed,
            file_name_suffix,
            thread,
        );
    }

    let mut primers_contents = detected_primers
        .iter()
        .map(|primer| format!("{}\t{}\t{}", primer.name, primer.fwd(), primer.rev()))
        .collect::<Vec<_>>()
        .join("\n");
    primers_contents.push('\n');
    std::fs::write(
        format!(
            "{output_dir}/{analysis_name}{}",
            file_name_suffix.detected_primers
        ),
        primers_contents.as_bytes(),
    )
    .expect("Failed to save detected primers");
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
    save_failed: bool,
    file_name_suffix: &FileNameSuffix,
    thread: usize,
) -> ReadsClassifier {
    let primers = primer.name2primer();
    let mut classifier = classifier;
    let total_reads = classifier.all_reads.len();

    // Step2: Trim sequence outside this known primer of amplicon,
    // try to get clean reads indexes that startswith fwd primer and ends with rev primer rc
    let primer_name2reads_idx = classifier.classify_reads_with_known_primers(&primers);

    // assert_eq!(primer_name2reads_idx.len(), 1);
    let (primer_name, reads_idx) = primer_name2reads_idx
        .into_iter()
        .collect::<Vec<_>>()
        .pop()
        .unwrap();

    // Step3: Move this clean reads from classifier into ReadsWithPrimer
    let reads = classifier.remove_reads_with_idxes(&reads_idx);
    let total_paired_reads = reads.len();
    info!(
        "{} reads with paired primer [{}] and used to construct ReadsWithPrimer",
        total_paired_reads, primer_name
    );
    let mut reads_with_primer = ReadsWithPairedPrimers::new(reads, reads_downsample);

    // Step4: Filter reads depending on read quality and length
    let failed_fastq = format!(
        "{output_dir}/{}_{primer_name}{}",
        classifier.analysis_name, file_name_suffix.bad_reads
    );
    let bad_reads_size = reads_with_primer.filter(
        min_read_quality,
        length_range,
        if save_failed {
            Some(&failed_fastq)
        } else {
            None
        },
    );
    info!(
        "{bad_reads_size}/{total_paired_reads} bad reads dropped{}",
        if save_failed {
            format!(", saved them into {}", failed_fastq)
        } else {
            "".to_string()
        }
    );

    // Step5: Save good reads used to run abpoa of this primer and clean remaining reads and redundant records
    check_and_create_dir(output_dir);
    let redundant_fastq = format!(
        "{output_dir}/{primer_name}{}",
        file_name_suffix.redundant_reads
    );
    reads_with_primer.save_redundant_fastq(&redundant_fastq);
    info!(
        "{}/{total_paired_reads} redundant reads found, saved into {redundant_fastq}",
        reads_with_primer.redundant_reads.len()
    );

    let this_primer_fastq = format!("{output_dir}/{primer_name}{}", file_name_suffix.good_reads);
    reads_with_primer.save_fastq(&this_primer_fastq);
    info!(
        "{}/{total_paired_reads} good reads used as input of abpoa, saved into {this_primer_fastq}",
        reads_with_primer.reads.len()
    );

    let clean_remaining_fastq = format!(
        "{output_dir}/{primer_name}{}",
        file_name_suffix.remaining_reads
    );
    classifier.save_fastq(&clean_remaining_fastq);
    info!(
        "{}={total_reads}-{total_paired_reads} remaining records, saved into {clean_remaining_fastq}",
        total_reads - total_paired_reads
    );

    // Step6: Run abpoa to construct draft consensus
    info!("Running abpoa for {primer_name}...");
    let draft_consensus = run_abpoa(&this_primer_fastq, output_dir, &primer_name, abpoa);

    /*
    Step7 and Step8 are used to remove reads that mapped to this amplicon draft consensus
    avoiding them to disturb next classifier.guess_one_primer.
    */
    // Step7: Map all remaining barcod trimmed fastq records to this draft consensus.
    info!(
        "Map remaining reads to draft_consensus to search reads that without paired primer detected at ends"
    );
    let sorted_bam = run_minimap2_and_index(
        output_dir,
        &clean_remaining_fastq,
        &draft_consensus,
        &primer_name,
        minimap2,
        samtools,
        thread,
    );

    // Step8: Remove mapped reads from classifier's all_reads
    let mapped_read_names = get_reads_name_from_bam(&sorted_bam, min_mapq);
    classifier.remove_reads_with_names(&mapped_read_names);
    info!(
        "{} reads found and remove them from ReadsClassifier. GO to next loop for next amplicon\n\n",
        mapped_read_names.len()
    );

    /*
    Theoretically, all reads from this amplicon have been removed from classifier here.
    1.  Some were removed by classifier.remove_reads_with_indexes.
        These reads have paired primers detected at dual ends.
        Some of these reads used to construct draft_consensus
    2.  Some were removed by classifier.remove_reads_with_names.
        Paired primers can't be detected at ends.
        I collect them by align them to draft_consensus. These reads just be removed.
     */

    // Step9: For next loop
    classifier
}

/*
TODO
In some cases, multi amplicons share same forward primer and different reverse primers. How extract one of them in many many similar amplicons?.
Consider mapq and alignment ratio

amplicon1: 5'----------------------------------------------------------------------------------------------------3'
amplicon2: 5'-------------------------------------------------------------------------3'
amplicon3: 5'--------------------------------------3'
if use amplicon2 as reference, when the real and good reads from amplicon2 mapped to this reference,
the start position of alignment on reference should be near 0 position and the end alignment position be near reference.len() position.
Also, alignment start and alignment end of this read should be near of read's 0 position and read.len() position
 */

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
    let save_failed = amp_cmd.get_one::<bool>("retain_failed").unwrap();
    let len_range = amp_cmd.get_one::<f64>("len_range").unwrap();
    let amp_numbers = amp_cmd.get_one::<usize>("number").unwrap();
    let lead_length = amp_cmd.get_one::<usize>("lead").unwrap();
    let detect_rev_primer_reads_number =
        amp_cmd.get_one::<usize>("detect_rev_primer_reads").unwrap();
    let min_mapq = amp_cmd.get_one::<u8>("min_mapq").unwrap();
    let abpoa = amp_cmd.get_one::<String>("abpoa").map(|x| x.as_ref());
    let minimap2 = amp_cmd.get_one::<String>("minimap2").map(|x| x.as_ref());
    let samtools = amp_cmd.get_one::<String>("samtools").map(|x| x.as_ref());
    let analysis_name = amp_cmd.get_one::<String>("prefix").unwrap();
    let thread = amp_cmd.get_one::<u16>("thread").unwrap();

    let primers_is_set = amp_cmd.value_source("primers") == Some(ValueSource::CommandLine);
    let amplicon_number_is_set = amp_cmd.value_source("number") == Some(ValueSource::CommandLine);
    if primers_is_set && amplicon_number_is_set {
        quit_with_error("--primers and number couldn't be specified together")
    }
    init_log();
    {
        eprintln!("Args...........");
        eprintln!("\t--input\t{input}");
        eprintln!("\t--output\t{output}");
        eprintln!("\t--primers\t{primers_opt:?}");
        eprintln!("\t--number\t{amp_numbers}");
        eprintln!("\t--barcode\t{barcode}");
        eprintln!("\t--left\t{left}");
        eprintln!("\t--right\t{right}");
        eprintln!("\t--distance\t{distance}");
        eprintln!("\t--downsample\t{downsample}");
        eprintln!("\t--min_qual\t{min_qual}");
        eprintln!("\t--len_range\t{len_range}");
        eprintln!("\t--retain_failed\t{save_failed}");
        eprintln!("\t--lead\t{lead_length}");
        eprintln!("\t--prefix\t{analysis_name}");
        eprintln!("\t--detect_rev_primer_reads\t{detect_rev_primer_reads_number}");
        eprintln!("\t--min_mapq\t{min_mapq}");
        eprintln!("\t--abpoa\t{abpoa:?}");
        eprintln!("\t--thread\t{thread}");
        eprintln!("\t--minimap2\t{minimap2:?}");
        eprintln!("\t--samtools\t{samtools:?}");
        eprintln!("Args...........");
    }
    let output_path = std::env::current_dir()
        .expect("Failed to get current directory path")
        .join(Path::new(output).to_owned());
    let output = output_path
        .to_str()
        .expect("Failed to set output directory");

    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(*thread as usize)
        .build_global();
    let file_name_suffix = FileNameSuffix::default();

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
            *downsample,
            output,
            *min_qual,
            *len_range,
            abpoa,
            minimap2,
            samtools,
            *min_mapq,
            analysis_name,
            *save_failed,
            &file_name_suffix,
            *thread as usize,
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
                *save_failed,
                &file_name_suffix,
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
            .default_value("5000")
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
        Arg::new("retain_failed")
            .long("retain_failed")
            .action(ArgAction::SetTrue)
            .help("whether to save reads with paired primers found at ends but with low quality or abnormal length")
    ).arg(
        Arg::new("lead")
            .long("lead")
            .default_value("21")
            .value_parser(value_parser!(usize))
            .help("when no known primers specified, use first lead bases as candidate fwd primer after barcode trimmed")
    ).arg(
        Arg::new("prefix")
            .long("prefix")
            .default_value("test001")
            .help("the prefix of output file name")
    ).arg(
        Arg::new("detect_rev_primer_reads")
            .long("detect_rev_primer_reads")
            .default_value("500")
            .value_parser(value_parser!(usize))
            .help("[UnknownPrimers]: how many reads used to detect rev primer for reads with candidate fwd primer")
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
        Arg::new("thread")
            .short('t')
            .long("thread")
            .default_value("4")
            .value_parser(value_parser!(u16).range(1..=32))
            .help("minimap2 thread")
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
    )
}
