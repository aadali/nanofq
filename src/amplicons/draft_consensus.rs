use crate::fastq2::FastqRecord;
use crate::primer_barcode::Primer;
use std::collections::HashMap;
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


