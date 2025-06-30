use ansi_term::Color;
use flate2::bufread::MultiGzDecoder;
use seq_io::fastq;
use seq_io::fastq::{Record, RefRecord};
use std::fs::File;
use std::io::{BufReader, Read, Stdin};
use std::ops::{Deref, DerefMut};

const BUFF: usize = 1024 * 1024;
pub type EachStats = (Box<String>, usize, (f64, f64), Option<f64>); // (f64, f64): (this_read_average_error_pro, this_read_quality)

pub trait ReadStats {
    fn gc_count(&self) -> f64;
    fn calculate_read_quality(&self) -> (f64, f64);

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
    fn calculate_read_quality(&self) -> (f64, f64) {
        let seq_len = self.seq().len() as f64;
        let avg_err_prob = self
            .qual()
            .iter()
            .map(|x| 10.0f64.powf((x - 33) as f64 / -10.0))
            .sum::<f64>()
            / seq_len;
        let quality = avg_err_prob .log10() * -10.0;
        (avg_err_prob, quality)
    }

    #[inline]
    fn stats(&self, gc: bool) -> EachStats {
        let len = self.seq().len();
        let read_quality = self.calculate_read_quality();
        let gc = if gc { Some(self.gc_count()) } else { None };

        (
            Box::new(
                self.id()
                    .expect(&Color::Red.paint("parse id to str error").to_string())
                    .to_string(),
            ),
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
            && self.calculate_read_quality().1 > min_read_qvalue
            && gc_passed
    }
}

pub struct FastqReader<R: Read>(pub fastq::Reader<R>);

impl<R: Read> FastqReader<R> {
    pub fn with_stdin() -> FastqReader<Stdin> {
        FastqReader(fastq::Reader::with_capacity(std::io::stdin(), BUFF))
    }

    pub fn with_fastq(path: &str) -> FastqReader<File> {
        FastqReader(fastq::Reader::with_capacity(
            File::open(path).expect(&format!(
                "{}: {}",
                Color::Red.paint("Open failed failed: "),
                path
            )),
            BUFF,
        ))
    }

    pub fn with_fastq_gz(path: &str) -> FastqReader<MultiGzDecoder<BufReader<File>>> {
        FastqReader(fastq::Reader::with_capacity(
            MultiGzDecoder::new(BufReader::new(File::open(path).expect(&format!(
                "{}: {}",
                Color::Red.paint("Open failed failed: "),
                path
            )))),
            BUFF,
        ))
    }

    #[inline]
    pub(crate) fn stats(&mut self, gc: bool) -> Vec<EachStats> {
        let mut stats_result = vec![];
        loop {
            let ref_record = self.next();
            if let Some(record) = ref_record {
                stats_result.push(
                    record
                        .expect(&Color::Red.paint("Error: ").to_string())
                        .stats(gc),
                )
            } else {
                break;
            }
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
