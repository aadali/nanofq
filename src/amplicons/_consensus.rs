use crate::fastq2::FastqRecord;
use crate::primer_barcode::Primer;
use crate::utils::quit_with_error;
use ahash::{HashMap, RandomState};
use std::sync::OnceLock;

const BIT2BASE: [u8; 4] = [b'A', b'C', b'G', b'T'];
// const A: u64 = 0b00;
// const C: u64 = 0b01;
// const G: u64 = 0b10;
// const T: u64 = 0b11;
const BASE2BIT: [u64; 128] = {
    let mut table = [0u64; 128];
    table[b'A' as usize] = 0b00;
    // table[b'a' as usize] = 0;
    table[b'C' as usize] = 0b01;
    // table[b'c' as usize] = 1;
    table[b'G' as usize] = 0b10;
    // table[b'g' as usize] = 2;
    table[b'T' as usize] = 0b11;
    // table[b't' as usize] = 3;
    table
};

pub type BitKmer = u64;

static K_SIZE: OnceLock<usize> = OnceLock::new();

pub trait KTrait {
    fn from<T: AsRef<[u8]>>(seq: T) -> BitKmer;

    fn to_string(&self) -> String;

    fn next_kmers(&self) -> [BitKmer; 4];

    fn prev_kmers(&self) -> [BitKmer; 4];
}

impl KTrait for BitKmer {
    fn from<T: AsRef<[u8]>>(seq: T) -> BitKmer {
        debug_assert_eq!(seq.as_ref().len(), *K_SIZE.get().unwrap());
        debug_assert!(*K_SIZE.get().unwrap() < 33);
        let bit_kmer = seq
            .as_ref()
            .iter()
            .fold(0u64, |acc, base| (acc << 2) | (BASE2BIT[*base as usize]));
        bit_kmer
    }

    fn to_string(&self) -> String {
        let mut seq = String::with_capacity(*K_SIZE.get().unwrap());
        for x in 0..*K_SIZE.get().unwrap() {
            let shift = (*K_SIZE.get().unwrap() - x - 1) * 2;
            let base = BIT2BASE[(*self >> shift & 0b11) as usize];
            seq.push(base as char);
        }
        seq
    }
    fn next_kmers(&self) -> [BitKmer; 4] {
        [
            self << 2 | 0b00,
            self << 2 | 0b01,
            self << 2 | 0b10,
            self << 2 | 0b11,
        ]
    }

    fn prev_kmers(&self) -> [BitKmer; 4] {
        [
            *self >> 2 | 0,
            (*self >> 2) | (0b01 << (*K_SIZE.get().unwrap() * 2 - 2)),
            (*self >> 2) | (0b10 << (*K_SIZE.get().unwrap() * 2 - 2)),
            (*self >> 2) | (0b11 << (*K_SIZE.get().unwrap() * 2 - 2)),
        ]
    }
}

pub struct ReadsWithPairedPrimers<'a> {
    good_reads: HashMap<usize, FastqRecord>,
    primer: &'a Primer,
}

impl<'a> ReadsWithPairedPrimers<'a> {
    pub fn new_with_known_primer(
        all_reads: &mut HashMap<usize, FastqRecord>,
        good_reads_idxes: &Vec<usize>,
        primer: &'a Primer,
    ) -> Self {
        let mut good_reads =
            HashMap::with_capacity_and_hasher(good_reads_idxes.len(), RandomState::new());
        for idx in good_reads_idxes {
            let read = all_reads.remove(idx);
            if read.is_none() {
                quit_with_error(
                    &format!("Failed to get read with index {idx} from all reads, is this read a good read?"),
                )
            }
            good_reads.insert(*idx, read.unwrap()).unwrap();
        }
        ReadsWithPairedPrimers { good_reads, primer }
    }

    pub fn new_with_unknown_primer(all_reads: &mut HashMap<usize, FastqRecord>) {}

    fn filter(&mut self, read_q: f64, length_range: f64) {
        let mean_length = (self.good_reads.iter().fold(0usize, |acc, x| acc + x.0) as f64)
            / (self.good_reads.len() as f64);
        self.good_reads.retain(|_, read| {
            read.qual(false) > (read_q as f32)
                && (read.len() as f64) < (mean_length * (1.0 + length_range))
                && (read.len() as f64) > (mean_length * (1.0 - length_range))
        });
    }

    fn write_to_file(&self, mode: &str) {
        println!("hello world");
        println!("hello java");
        println!("hello python");
        println!("hello cpp");
        println!("hello R");
        println!("hello groovy");
        println!("hello rust");
        println!("hello nextflow");
        println!("hello world");
        println!("hello javascript");
        println!("hello typescript");
        println!("hello ruby");
    }
}
