use crate::amplicons::preprocess::ReadsWithPrimer;
use crate::primer_barcode::{BARCODES, Barcode, Primer};
use crate::utils::{check_and_create_dir, run_abpoa};
use ahash::{HashMap, RandomState};
use rayon::prelude::*;

mod consensus;
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
                    primers.get(&primer_name).unwrap().clone(),
                    reads_downsample,
                ),
            )
        })
        .collect::<Vec<_>>();

    check_and_create_dir(output_dir);
    primer_reads
        .into_par_iter()
        .for_each(|(primer_name, mut reads_with_primer)| {
            let save_fastq_path = format!("{output_dir}/{primer_name}.fastq");
            reads_with_primer.filter(min_read_quality, length_range);
            reads_with_primer.save(&save_fastq_path);
            run_abpoa(&save_fastq_path, output_dir, &primer_name, abpoa)
        })
}

fn amplicon_with_unknown_primers(input_fastq: &str, bar_idx: usize) {}
