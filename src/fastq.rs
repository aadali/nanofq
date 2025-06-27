use seq_io::fastq;
use seq_io::fastq::{Record,  RefRecord};
use std::io::{Read};
use std::ops::{Deref, DerefMut};
use std::path::Path;
use anyhow::anyhow;

const BUFF: usize = 1024 * 1024;
pub type EachStats = (Box<String>, usize, f64, Option<f64>);

pub trait ReadStats {
    fn gc_count(&self) -> f64;
    fn calculate_read_quality(&self) -> f64;

    fn stats(&self, gc: bool) -> EachStats;

    fn is_passed(
        &self,
        min_len: usize,
        max_len: usize,
        min_read_qvalue: f64,
        gc: bool,
        min_gc: f64,
        max_gc: f64,
    ) -> bool;
}

impl<'a> ReadStats for RefRecord<'a> {
    #[inline]
    fn gc_count(&self) -> f64 {
        let seq_len = self.seq().len() as f64;
        let gc_number: usize = self
            .seq()
            .iter()
            .map(|x| {
                if x == &b'G' || x == &b'C' || x == &b'g' || x == &b'c' {
                    1usize
                } else {
                    0usize
                }
            })
            .sum();
        gc_number as f64 / seq_len
    }

    #[inline]
    fn calculate_read_quality(&self) -> f64 {
        let seq_len = self.seq().len() as f64;
        (self
            .qual()
            .iter()
            .map(|x| 10.0f64.powf((x - 33) as f64 / -10.0))
            .sum::<f64>()
            / seq_len)
            .log10()
            * -10.0
    }

    #[inline]
    fn stats(&self, gc: bool) -> EachStats {
        let len = self.seq().len();
        let read_quality = self.calculate_read_quality();
        let gc = if gc { Some(self.gc_count()) } else { None };

        (
            Box::new(self.id().expect("parse id to str error").to_string()),
            len,
            read_quality,
            gc,
        )
    }

    #[inline]
    fn is_passed(
        &self,
        min_len: usize,
        max_len: usize,
        min_read_qvalue: f64,
        gc: bool,
        min_gc: f64,
        max_gc: f64,
    ) -> bool {
        let seq_len = self.seq().len();
        let gc_passed = if gc {
            let gc = self.gc_count();
            gc > min_gc && gc < max_gc
        } else {
            true
        };
        seq_len >= min_len
            && seq_len <= max_len
            && self.calculate_read_quality() > min_read_qvalue
            && gc_passed
    }
}

pub struct FastqReader<R: Read>(pub fastq::Reader<R>);

trait IsFastqReader {}

impl<R: Read> IsFastqReader for FastqReader<R> {}

impl<R: Read> FastqReader<R> {
    fn new<P: AsRef<Path>>(path: Option<P>) -> Result<Box<dyn IsFastqReader>, anyhow::Error> {
        if path.is_none() {
            Ok(Box::new(FastqReader(fastq::Reader::with_capacity(
                std::io::stdin(),
                BUFF,
            ))))
        } else {
            let path = path.unwrap();
            if path.as_ref().ends_with(".fastq") || path.as_ref().ends_with(".fq") {
                Ok(Box::new(FastqReader(fastq::Reader::with_capacity(
                    std::fs::File::open(path)?,
                    BUFF,
                ))))
            } else if path.as_ref().ends_with(".fastq.gz") || path.as_ref().ends_with(".fq.gz") {
                Ok(Box::new(FastqReader(fastq::Reader::with_capacity(
                    flate2::bufread::MultiGzDecoder::new(std::io::BufReader::new(
                        std::fs::File::open(path)?,
                    )),
                    BUFF,
                ))))
            } else {
                Err(anyhow!("bad file suffix"))
            }
        }
    }
    
    #[inline]
    pub(crate) fn stats(&mut self, gc: bool) -> Vec<EachStats> {
        let mut stats_result = vec![];
        loop {
            let ref_record = self.next();
            if ref_record.is_none() {
                break;
            }
            stats_result.push(ref_record.unwrap().expect("Error").stats(gc))
        }
        stats_result
    }
}

impl<R: Read> Deref for FastqReader<R> {
    type Target = fastq::Reader<R>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<R: Read> DerefMut for FastqReader<R> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
