mod collect_reads;
mod classsify_reads;
mod draft_consensus;

use crate::fastq2::{FastqRecord, read_fastq};
use crate::primer_barcode::{Barcode, PO, Primer};
use ahash::{HashMap, RandomState};
use bio::pattern_matching::myers::Myers;

//     let start_primer_str = "GCAACAACAACCTTTCATCCT".to_owned();
//     let start_primer_str = "GCAACAACAACCTTTCATCCTAATTCTGG".to_owned();
//     let start_primer_str = "AGGGAAACACGGTAAGATCCGAACAGCACCTAGGGAAACACGATAGAATCCGAACAGCACCTGCAACAACAACCTTTCATCCT".to_owned();

//     let end_primer_str = "ATTTGACAGGATTTATGTGTA".to_owned();


