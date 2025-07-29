use crate::fastq::{FastqReader, NanoRead};
use crate::utils::rev_com;
use seq_io::fastq::{OwnedRecord, Record, RefRecord};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::io::{Read, Write};

struct SubRecord<T> {
    ref_record: T,
    start: usize,
    end: usize,
}

impl<T> Record for SubRecord<T>
where
    T: Record,
{
    fn head(&self) -> &[u8] {
        self.ref_record.head()
    }

    fn seq(&self) -> &[u8] {
        &self.ref_record.seq()[self.start..=self.end]
    }

    fn qual(&self) -> &[u8] {
        &self.ref_record.seq()[self.start..=self.end]
    }
}
impl<'a> SubRecord<RefRecord<'a>> {
    fn write_rev_com(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
        write!(
            writer,
            "@{}\n{}\n+\n{}\n",
            format!(
                "{}_rev_com {}",
                self.ref_record.id()?,
                self.desc().unwrap_or(Ok("")).unwrap()
            ),
            rev_com(unsafe { std::str::from_utf8_unchecked(self.seq()) }),
            unsafe {
                String::from_utf8_unchecked(
                    self.qual().iter().rev().map(|x| *x).collect::<Vec<u8>>(),
                )
            }
        )?;
        Ok(())
    }

    fn to_owned(&self) -> SubRecord<OwnedRecord> {
        SubRecord {
            ref_record: self.ref_record.to_owned_record(),
            start: self.start,
            end: self.end,
        }
    }
}

struct QualOwnedRecord {
    qual: f64,
    sub_record: SubRecord<OwnedRecord>,
}

impl Ord for QualOwnedRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Eq for QualOwnedRecord {}
impl PartialOrd for QualOwnedRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        other.qual.partial_cmp(&self.qual)
    }
}

impl PartialEq for QualOwnedRecord {
    fn eq(&self, other: &Self) -> bool {
        (self.qual * 1_000_000.0) as usize == (other.qual * 1_000_000.00) as usize
    }
}
fn push_heap(
    ref_record: RefRecord,
    trim_from: usize,
    trim_to: usize,
    best_n_reads: &mut BinaryHeap<QualOwnedRecord>,
    end3: &str,
    min_len: usize,
    n: usize,
) {
    let trim_to = trim_to + end3.len() - 1;
    if !(trim_to <= trim_from || trim_to - trim_from < min_len) {
        let sub_record = SubRecord {
            ref_record: ref_record.to_owned_record(),
            start: trim_from,
            end: trim_to,
        };
        let this_qual = sub_record.calculate_read_quality().1;
        let heap_item = QualOwnedRecord {
            qual: this_qual,
            sub_record,
        };
        if best_n_reads.len() < n {
            best_n_reads.push(heap_item)
        } else if heap_item > *best_n_reads.peek().unwrap() {
            best_n_reads.pop();
            best_n_reads.push(heap_item);
        }
    }
}
pub fn clean_amplicon<R: Read>(
    fastq_reader: &mut FastqReader<R>,
    fwd_primer: &str,
    rev_primer: &str,
    min_len: usize,
    n: usize,
    output_writer: &mut dyn Write,
) {
    let end5 = fwd_primer;
    let end3 = rev_com(rev_primer);
    let rev_com_end5 = rev_primer;
    let rev_com_end3 = rev_com(fwd_primer);
    let mut best_n_reads = BinaryHeap::<QualOwnedRecord>::new();
    loop {
        if let Some(Ok(ref_record)) = fastq_reader.next() {
            if ref_record.seq().len() < min_len {
                continue;
            }
            let seq_str = unsafe { std::str::from_utf8_unchecked(ref_record.seq()) };
            if let Some(trim_from) = seq_str.find(end5) {
                if let Some(trim_to_) = seq_str.find(&end3) {
                    push_heap(
                        ref_record,
                        trim_from,
                        trim_to_,
                        &mut best_n_reads,
                        &end3,
                        min_len,
                        n,
                    );
                }
            } else {
                if let Some(trim_from) = seq_str.find(rev_com_end5) {
                    if let Some(trim_to_) = seq_str.find(&rev_com_end3) {
                        push_heap(
                            ref_record,
                            trim_from,
                            trim_to_,
                            &mut best_n_reads,
                            &rev_com_end3,
                            min_len,
                            n,
                        );
                    }
                }
            }
        } else {
            break;
        }
    }
    let best_n_reads = best_n_reads.into_sorted_vec();
    for x in &best_n_reads {
        NanoRead::write(&x.sub_record, output_writer).unwrap();
    }
}

pub fn run_amplicon() {
    let fwd= "CAAGATCGTGCCACGGTACTCCAGCC";
    let rev = "CAGACCAGAGGATTTCGGAGGGTCTG";
    let raw_fastq = "/Users/aadali/test_data/amplicon/barcode03.fastq";
    let mut reader = FastqReader::new(std::fs::File::open(raw_fastq).unwrap());
    let min_len = 1200;
    let n = 2000;
    let mut bar03_trimmed_fq = std::io::BufWriter::new(
        std::fs::File::create("/Users/aadali/test_data/amplicon/barcode03_primer.trimmed.fastq")
            .unwrap(),
    );
    clean_amplicon(&mut reader, fwd, rev, min_len, n, &mut bar03_trimmed_fq);
}