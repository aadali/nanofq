use crate::filter::FilterOption;
use crate::utils::{calculate_read_q, gc};
use needletail::{Sequence, parse_fastx_file};
use std::fmt::{Display, Formatter};
use std::io::Write;
use std::sync::mpsc;
use std::sync::mpsc::Receiver;
use std::thread::JoinHandle;
use std::{io, str, thread};

#[derive(Debug)]
pub struct FastqRecord {
    name: String,
    description: Option<String>,
    seq: Vec<u8>,
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
            qual: calculate_read_q(self.quality, use_dorado_q),
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
        calculate_read_q(&self.quality, use_dorado_q)
    }

    pub fn gc(&self) -> f32 {
        gc(&self.seq)
    }

    pub fn write(&self, writer: &mut dyn Write) -> Result<(), io::Error> {
        write!(
            writer,
            "@{}{}\n{}\n+\n{}\n",
            self.name,
            if self.description.is_some() {
                self.description.as_ref().unwrap().as_ref()
            } else {
                ""
            },
            unsafe { str::from_utf8_unchecked(&self.seq) },
            unsafe { str::from_utf8_unchecked(&self.quality) }
        )?;
        Ok(())
    }
}

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

fn read_fastq(fastq_file: &str, need_description: bool) -> Vec<FastqRecord> {
    let mut fastq_records = vec![];
    let mut records = parse_fastx_file(fastq_file).expect(&format!("Failed to read {fastq_file}"));
    let mut read_idx = 1;
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
