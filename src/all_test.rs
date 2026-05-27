#![allow(unused)]
use crate::amplicons::{amplicon_with_known_primers, amplicon_with_unknown_primers, preprocess};
use crate::primer_barcode::{BARCODES, Barcode, Primer};
use bio::alphabets::dna::revcomp;

pub fn test_amplicon() {
    let primer_str = "py_bar04_1600,GCAACAACAACCTTTCATCCT,TACACATAAATCCTGTCAAAT;py_bar03_1200,CAACAACAACCTTTCATCCTA,TTGAGGACTGGAAAGATCAAT";
    // let fwd_primer = "CAGCACCTCCAGGGTTTTCC";
    // let fwd_primer = "CAGCACCTCCAGGGTTTTCCC";
    // let fwd_primer = "CCAGGGTTTTCCCAGTCACGA";
    // let fwd_primer = "CCAGGGTTTTCCCAGTCACGA";
    // let rev_primer = "CCAGGGTTTTCCCAGTCACGA";
    // let rev_primer_rc = "CACACAACATACGAGCCGGAA";

    // let raw_fastq = "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode04-1600.fastq";
    let raw_fastq = "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode03-barcode04-merged.fastq";
    let primers = Primer::parse_primer_from_cli(primer_str);
    amplicon_with_known_primers(
        raw_fastq,
        4,
        150,
        150,
        3,
        primers,
        10000,
        "/Users/aadali/projects/RustProjects/nanoamp/test_data/output",
        15.0,
        0.05,
        Some("/Users/aadali/biotools/abPOA-v1.5.6_arm64-macos/abpoa")
    );
    // let barcode = Barcode::new(BARCODES[4].as_bytes());
    // let start_primer_str = "GCAACAACAACCTTTCATCCT".to_owned();
    // // let start_primer_str = "GCAACAACAACCTTTCATCCTAATTCTGG".to_owned();
    // let end_primer_str = "ATTTGACAGGATTTATGTGTA".to_owned();
    // let end_primer_str = "TACACATAAATCCTGTCAAAT".to_owned();
    // println!(
    //     "{}",
    //     str::from_utf8(revcomp(end_primer_str.as_bytes()).as_slice()).unwrap()
    // );
    // let primer = Primer::new(
    //     "py1600",
    //     start_primer_str.as_bytes(),
    //     end_primer_str.as_bytes(),
    // );
    // let reads_collector = preprocess::ReadsCollector::new(
    //     "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode04-1600.fastq",
    //     &barcode,
    //     150,
    //     150,
    //     3,
    // );
    //
    // let all_reads = reads_collector.collect_fastqs();
    // println!("{}", all_reads.len());
    // let bar_trimmed_fastq = "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode04-1600-barcode-trimmed.fastq";
    // let mut trimmed_output = std::fs::File::create(bar_trimmed_fastq).unwrap();
    // for (x, y) in all_reads.iter() {
    //     y.write(&mut trimmed_output).unwrap();
    // }




    // let mut classifier = preprocess::ReadsClassifier::new(all_reads, 150, 150, 3);
    // let reads = classifier.classify_reads_with_known_primers(&primers);
    // for (primer_name, reads_idx) in reads {
    //     let reads = classifier.remove_reads(&reads_idx);
    //     let mut reads_with_primer = preprocess::ReadsWithPrimer::new(
    //         reads, primers.get(&primer_name).unwrap().clone()
    //     );
    //     reads_with_primer.filter(12.0, 0.05);
    //     reads_with_primer.save("/Users/aadali/projects/RustProjects/nanofq/test_data/py-barcode04-1600-filtered.fastq");
    // }




    // let primer_trimmed_fastq = "/Users/aadali/projects/RustProjects/nanofq/test_data/py-barcode04-1600-primer-trimmed.fastq";
    // let mut primer_trimmed_output = std::fs::File::create(primer_trimmed_fastq).unwrap();
    // let mut count = 0;
    // for x in reads.iter() {
    //     for read_idx in x.1 {
    //         let read = classifier.all_reads.get(read_idx).unwrap();
    //
    //         if read.seq.starts_with("GCAACAACAACCTTTCATCCT".as_bytes()) {
    //             count += 1;
    //         }
    //         classifier
    //             .all_reads
    //             .get(read_idx)
    //             .unwrap()
    //             .write(&mut primer_trimmed_output)
    //             .unwrap();
    //     }
    // }
    println!("hello world");
}

pub fn test_amplicons_with_unknown_primer() {
    let input_fastq = "/Users/aadali/projects/RustProjects/nanoamp/test_data/py-barcode03-barcode04-merged.fastq";
    amplicon_with_unknown_primers(
        input_fastq,
        3,
        150,
        150,
        3,
        2,
        21,
        1000,
        0.5,
        100,
        10000,
        "/Users/aadali/projects/RustProjects/nanoamp/test_data/output2",
        14.0,
        0.05,
        Some("abpoa"),
        Some("minimap2"),
        Some("samtools"),
        50
    )
}
