use crate::amplicons::preprocess::{ReadsClassifier, ReadsWithPrimer};
use crate::primer_barcode::{BARCODES, Barcode, Primer};
use crate::utils::{check_and_create_dir, quit_with_error, run_abpoa, run_minimap2_and_index};
use ahash::{HashMap, HashSet, RandomState};
use rayon::prelude::*;
use rust_htslib::bam::{IndexedReader, Read, Record};
use std::io::{BufWriter, Write};

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
    let mut classifier =
        preprocess::ReadsClassifier::new(all_reads, 21, left_range, right_range, max_distance);
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
            let save_fastq_path = format!("{output_dir}/{primer_name}.trimmed_filtered.fastq");
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
    );

    for primer_idx in 1..amplicons_number + 1 {
        // Step1.1: guess one paired primer from barcode trimmed fastq
        let primer_opt = classifier.guess_one_primer(
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
                let primer_opt2 = classifier.guess_one_primer(
                    guess_reads_number,
                    min_rev_match_reads_number,
                    rev_primer_found_ratio,
                    primer_idx,
                    true,
                    10,
                );
                match primer_opt2 {
                    None => quit_with_error(&format!("Couldn't detect {primer_idx}th primer")),
                    Some(primer) => primer,
                }
            }
            Some(primer) => primer,
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
            primer_idx,
        );
    }
    let mut output_primers = BufWriter::new(
        std::fs::File::create(format!("{output_dir}/detected_primers.txt")).unwrap(),
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
    remain_idx: usize,
) -> ReadsClassifier {
    let primers = primer.name2primer();
    let mut classifier = classifier;

    // Step2: Trim sequence outside this known primer of amplicon,
    // try to get clean reads indexes that startswith fwd primer and ends with rev primer rc
    let primer_name2reads_idx = classifier.classify_reads_with_known_primers(&primers);

    assert_eq!(primer_name2reads_idx.len(), 1);
    let (primer_name, reads_idx) = primer_name2reads_idx
        .into_iter()
        .collect::<Vec<_>>()
        .pop()
        .unwrap();

    // Step3: Remove this clean reads from classifier and save them in ReadsWithPrimer
    let reads = classifier.remove_reads_with_idxes(&reads_idx);
    let mut reads_with_primer = ReadsWithPrimer::new(reads, reads_downsample);

    /*
    TODO for reads_with_primer.filter, if read is filtered, insert this read into classifier again.
     This read will be mapped to draft_consensus and not be used for next loop
     */
    // Step4: Filter reads depending on read quality and length
    reads_with_primer.filter(min_read_quality, length_range);

    // Step5: Save this clean fastq and remaining barcode trimmed fastq
    check_and_create_dir(output_dir);
    let save_fastq_path = format!("{output_dir}/{primer_name}.trimmed_filtered.fastq");
    reads_with_primer.save_fastq(&save_fastq_path);
    let barcode_trimmed_fastq = format!("{output_dir}/remaining_{remain_idx}.fastq");
    classifier.save_fastq(&barcode_trimmed_fastq);

    // Step6: Run abpoa to construct draft consensus
    let draft_consensus = run_abpoa(&save_fastq_path, output_dir, &primer_name, abpoa);

    /*
    Step7 and Step8 are used to remove reads that mapped to this amplicon draft consensus
    avoiding them to disturb next classifier.guess_one_primer.
    */
    // Step7: Map all remaining barcod trimmed fastq records to this draft consensus.
    let sorted_bam = run_minimap2_and_index(
        output_dir,
        &barcode_trimmed_fastq,
        &draft_consensus,
        &primer_name,
        minimap2,
        samtools,
    );

    // Step8: Remove mapped reads from classifier's all_reads
    let mapped_read_names = get_reads_name_from_bam(&sorted_bam, min_mapq);
    classifier.remove_reads_with_names(&mapped_read_names);
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
    let mut bam_records =
        IndexedReader::from_path(bam_file).expect(&format!("Failed to read {bam_file}"));
    let mut bam_record = Record::new();
    while let Some(Ok(_)) = bam_records.read(&mut bam_record) {
        if bam_record.mapq() >= min_mapq {
            let read_name = str::from_utf8(bam_record.qname()).unwrap();
            if !read_names.contains(read_name) {
                read_names.insert(read_name.to_string());
            }
        }
    }
    read_names
}
