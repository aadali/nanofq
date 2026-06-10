use crate::filter::FilterOption;
use crate::utils::{calculate_quality, complement, find_most_left_rear, find_most_right_front, gc};
use bio::pattern_matching::myers::Myers;
use needletail::{Sequence, parse_fastx_file};
use std::fmt::{Display, Formatter};
use std::io::Write;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread::JoinHandle;
use std::{io, str, thread};

#[derive(Debug)]
pub struct FastqRecord {
    pub name: String,
    description: Option<String>,
    pub seq: Vec<u8>,
    quality: Vec<u8>,
}
impl FastqRecord {
    pub fn new<T, U>(name: T, description: Option<T>, seq: U, quality: U) -> Self
    where
        T: Into<String>,
        U: Into<Vec<u8>>,
    {
        FastqRecord {
            name: name.into(),
            description: description.map(|x| x.into()),
            seq: seq.into(),
            quality: quality.into(),
        }
    }

    pub fn stats(self, use_dorado_q: bool, use_gc: bool) -> RecordEachStats {
        RecordEachStats {
            name: self.name,
            length: self.seq.len() as u32,
            qual: calculate_quality(self.quality, use_dorado_q, false),
            gc: if use_gc { Some(gc(self.seq)) } else { None },
        }
    }

    pub fn is_passed(&self, fo: &FilterOption) -> bool {
        let read_length = self.seq.len() as u32;
        if read_length > fo.max_len || read_length < fo.min_len {
            return false;
        }
        let read_qual = self.qual(fo.use_dorado_q);
        if read_qual > fo.max_qual || read_qual < fo.min_qual {
            return false;
        }
        if fo.use_gc {
            let gc = self.gc();
            if gc > fo.max_gc || gc < fo.min_gc {
                return false;
            }
        }
        true
    }

    pub fn len(&self) -> u32 {
        self.seq.len() as u32
    }

    pub fn qual(&self, use_dorado_q: bool) -> f32 {
        calculate_quality(&self.quality, use_dorado_q, false)
    }

    pub fn gc(&self) -> f32 {
        gc(&self.seq)
    }

    pub fn reversed(&mut self) {
        self.name.push_str("_rc");
        self.seq
            .iter_mut()
            .for_each(|base| *base = complement(*base));
        self.seq.reverse();
        self.quality.reverse();
    }

    pub fn write(&self, writer: &mut dyn Write) -> Result<(), io::Error> {
        let (sep, description) = if self.description.is_some() {
            (" ", self.description.as_ref().unwrap().as_ref())
        } else {
            ("", "")
        };
        write!(
            writer,
            "@{}{}{}\n{}\n+\n{}\n",
            self.name,
            sep,
            description,
            unsafe { str::from_utf8_unchecked(&self.seq) },
            unsafe { str::from_utf8_unchecked(&self.quality) }
        )?;
        Ok(())
    }

    pub fn split_off_front_barcode_end(
        &mut self,
        front_bar_par: &mut Myers,
        left_range: usize,
        max_distance: u8,
    ) -> (bool, bool) {
        let mut is_trimmed = false ;
        let search_seq = if left_range < self.seq.len() {
            &self.seq[..left_range]
        } else {
            self.seq.as_slice()
        };
        let matches = front_bar_par
            .find_all(search_seq, max_distance)
            .collect::<Vec<_>>();
        if let Some((_, idx, _)) = find_most_right_front(matches.clone(), max_distance) {
            self.quality = self.quality.split_off(idx);
            self.seq = self.seq.split_off(idx);
            debug_assert_eq!(self.seq.len(), self.quality.len());
            is_trimmed = true
        }
        (self.len() ==0, is_trimmed)
    }

    pub fn truncate_at_rear_barcode_start(
        &mut self,
        rear_bar_par: &mut Myers,
        right_range: usize,
        max_distance: u8,
    ) -> (bool, bool) {
        let mut is_truncated = false;
        let (search_seq, real_right_range) = if right_range < self.seq.len() {
            (&self.seq[self.seq.len() - right_range..], right_range)
        } else {
            (&self.seq[..], self.seq.len())
        };
        let matches = rear_bar_par
            .find_all(search_seq, max_distance)
            .collect::<Vec<_>>();
        if let Some((idx, _, _)) = find_most_left_rear(matches, max_distance) {
            self.quality
                .truncate(self.seq.len() - real_right_range + idx);
            self.seq.truncate(self.seq.len() - real_right_range + idx);
            debug_assert_eq!(self.seq.len(), self.quality.len());
            is_truncated = true;
        }
        (self.len() == 0, is_truncated)
    }

    pub fn _find_fwd_primer(
        &self,
        fwd_primer_pat: &mut Myers,
        left_range: usize,
        max_distance: u8,
    ) -> bool {
        let search_seq = if left_range < self.seq.len() {
            &self.seq[..left_range]
        } else {
            self.seq.as_slice()
        };
        fwd_primer_pat
            .find_all(search_seq, max_distance)
            .next()
            .is_some()
    }

    pub fn find_rev_primer(
        &self,
        rev_primer_pat: &mut Myers,
        right_range: usize,
        max_distance: u8,
    ) -> bool {
        let search_seq = if right_range < self.seq.len() {
            &self.seq[self.seq.len() - right_range..]
        } else {
            &self.seq[..]
        };
        rev_primer_pat
            .find_all(search_seq, max_distance)
            .next()
            .is_some()
    }

    pub fn split_off_fwd_primer(
        &mut self,
        fwd_primer_pat: &mut Myers,
        left_range: usize,
        max_distance: u8,
    ) -> bool {
        let search_seq = if left_range < self.seq.len() {
            &self.seq[..left_range]
        } else {
            self.seq.as_slice()
        };
        let match_position = fwd_primer_pat
            .find_all(search_seq, max_distance)
            .min_by_key(|x| x.2);
        // println!("{}, {}, {}", self.name, )
        if let Some(position) = match_position {
            let fwd_primer_start_idx = position.0;
            self.seq = self.seq.split_off(fwd_primer_start_idx);
            self.quality = self.quality.split_off(fwd_primer_start_idx);
            debug_assert_eq!(self.seq.len(), self.quality.len());
            true
        } else {
            false
        }
    }

    pub fn truncate_at_rev_primer_start(
        &mut self,
        rev_primer_pat: &mut Myers,
        right_range: usize,
        max_distance: u8,
    ) -> bool {
        let (search_seq, real_right_range) = if right_range < self.seq.len() {
            (&self.seq[self.seq.len() - right_range..], right_range)
        } else {
            (&self.seq[..], self.seq.len())
        };
        let match_position = rev_primer_pat
            .find_all(search_seq, max_distance)
            .min_by_key(|x| x.2);
        if let Some(position) = match_position {
            let rev_primer_end_idx = position.1;
            self.quality
                .truncate(self.seq.len() - real_right_range + rev_primer_end_idx);
            self.seq
                .truncate(self.seq.len() - real_right_range + rev_primer_end_idx);
            debug_assert_eq!(self.seq.len(), self.quality.len());
            true
        } else {
            false
        }
    }
}

#[derive(Clone)]
pub struct RecordEachStats {
    pub name: String,
    pub length: u32,
    pub qual: f32,
    pub gc: Option<f32>,
}

impl Display for RecordEachStats {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}", self.name, self.length, self.qual)
    }
}
impl RecordEachStats {
    pub fn new<T: Into<String>>(name: T, length: usize, qual: f32, gc: Option<f32>) -> Self {
        RecordEachStats {
            name: name.into(),
            length: length as u32,
            qual,
            gc,
        }
    }
}

impl Display for FastqRecord {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "@{}\n{}\n+\n{}",
            self.name,
            str::from_utf8(&self.seq).unwrap(),
            str::from_utf8(&self.quality).unwrap()
        )
    }
}

pub fn read_fastq(fastq_file: &str, need_description: bool) -> Vec<FastqRecord> {
    let mut fastq_records = vec![];
    let mut records = parse_fastx_file(fastq_file).expect(&format!("Failed to read {fastq_file}"));
    let mut read_idx = 1u32;
    while let Some(Ok(record)) = records.next() {
        let mut headers = record.id().splitn(2, |x| x.is_ascii_whitespace());
        let name = headers
            .next()
            .expect(&format!("Parse read name failed at {read_idx}th record"));

        let name =
            str::from_utf8(name).expect(&format!("Parse read name failed at {read_idx}th record"));

        let description;
        if need_description {
            description = headers
                .next()
                .map(|x| str::from_utf8(x).unwrap_or_default());
        } else {
            description = None
        }

        let seq = record.sequence();
        let quals = record
            .qual()
            .expect(&format!("Parse quality failed at {read_idx}th record"));
        fastq_records.push(FastqRecord::new(name, description, seq, quals));
        read_idx += 1;
    }
    fastq_records
}

pub fn chunk_records_from_fastq(
    fastq_file: &str,
    chunk: u32,
    need_description: bool,
) -> (JoinHandle<()>, Receiver<Vec<FastqRecord>>) {
    let (sender, receiver) = mpsc::sync_channel::<Vec<FastqRecord>>(1000);
    let mut records = parse_fastx_file(fastq_file).expect(&format!("Failed to read {fastq_file}"));
    let mut fastq_records = Vec::with_capacity(chunk as usize);
    let mut read_idx = 1;
    let handle = thread::spawn(move || {
        while let Some(Ok(record)) = records.next() {
            let mut headers = record.id().splitn(2, |x| x.is_ascii_whitespace());
            let name = headers
                .next()
                .expect(&format!("Parse read name failed at {read_idx}th record"));

            let name = str::from_utf8(name)
                .expect(&format!("Parse read name failed at {read_idx}th record"));

            let description;
            if need_description {
                description = headers
                    .next()
                    .map(|x| str::from_utf8(x).unwrap_or_default());
            } else {
                description = None
            }

            let seq = record.sequence();
            let quals = record
                .qual()
                .expect(&format!("Parse quality failed at {read_idx}th record"));
            fastq_records.push(FastqRecord::new(name, description, seq, quals));
            if fastq_records.len() as u32 >= chunk {
                let mut empty_records = Vec::with_capacity(chunk as usize);
                std::mem::swap(&mut empty_records, &mut fastq_records);
                if sender.send(empty_records).is_err() {
                    break;
                };
            }
            read_idx += 1;
        }
        sender
            .send(fastq_records)
            .expect("Failed to send last chunk");
    });
    (handle, receiver)
}
