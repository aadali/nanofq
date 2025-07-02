use ansi_term::Color;
use seq_io::fastq;
use seq_io::fastq::{Record, RefRecord};
use std::fmt::{Display};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::ops::{Deref, DerefMut};

const BUFF: usize = 1024 * 1024;
pub type EachStats = (Box<String>, usize, (f64, f64), Option<f64>); // (f64, f64): (this_read_average_error_pro, this_read_quality)

#[derive(Clone)]
pub(crate) struct FilterOption {
    pub min_len: usize,
    pub max_len: usize,
    pub min_qual: f64,
    pub max_qual: f64,
    pub gc: bool,
    pub min_gc: f64,
    pub max_gc: f64,
    // pub retain_failed: Option<String>,
}
impl FilterOption {
    pub(crate) fn set_failed_fastq_file(file_path: Option<String>) -> Result<Option<BufWriter<File>>, anyhow::Error> {
        if file_path.is_none() {
            Ok(None)
        } else {
            Ok(Some(BufWriter::new(File::create_new(file_path.unwrap())?)))
        }
    }
}
pub trait ReadStats {
    fn gc_count(&self) -> f64;
    fn calculate_read_quality(&self) -> (f64, f64);

    fn stats(&self, gc: bool) -> EachStats;

    fn is_passed(&self, fo: &FilterOption) -> bool;
    
    fn write<W: Write>(&self, writer: &mut W) -> Result<(), anyhow::Error>;
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
        let quality = avg_err_prob.log10() * -10.0;
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
    fn is_passed(&self, fo: &FilterOption) -> bool {
        let seq_len = self.seq().len();
        let gc_passed = if fo.gc {
            let gc = self.gc_count();
            gc > fo.min_gc && gc < fo.max_gc
        } else {
            true
        };
        let this_read_qual = self.calculate_read_quality().1;
        seq_len >= fo.min_len
            && seq_len <= fo.max_len
            && this_read_qual > fo.min_qual
            && this_read_qual < fo.max_qual
            && gc_passed
    }
    
    fn write<W:Write>(&self, writer: &mut W) -> Result<(), anyhow::Error>{
        unsafe {
            write!(
                writer,
                "@{}\n{}\n+\n{}\n",
                std::str::from_utf8_unchecked(self.head()),
                std::str::from_utf8_unchecked(self.seq()),
                std::str::from_utf8_unchecked(self.qual())
            )?;
        }
        Ok(())
    }
}

pub struct FastqReader<R: Read>(pub fastq::Reader<R>);

impl<R: Read> FastqReader<R> {
    pub fn new(reader: R) -> FastqReader<R> {
        FastqReader::<R>(fastq::Reader::with_capacity(reader, BUFF))
    }

    #[inline]
    pub(crate) fn stats(&mut self, gc: bool) -> Vec<EachStats> {
        let mut stats_result = vec![];
        loop {
            let ref_record_opt = self.next();
            if let Some(ref_record) = ref_record_opt {
                stats_result.push(
                    ref_record
                        .expect(
                            &Color::Red
                                .paint("Error: failed to get fastq record")
                                .to_string(),
                        )
                        .stats(gc),
                )
            } else {
                break;
            }
        }
        stats_result
    }

    pub(crate) fn filter<W: Write>(
        &mut self,
        writer: &mut BufWriter<W>,
        fo: &FilterOption,
    ) -> Result<(), anyhow::Error> {
        loop {
            match self.next() {
                Some(ref_record_res) => {
                    let ref_record = ref_record_res.expect(
                        &Color::Red
                            .paint("Error: failed to get fastq record")
                            .to_string(),
                    );
                    if ref_record.is_passed(fo) {
                        ReadStats::write(&ref_record, writer)?
                    }
                }
                None => {
                    return Ok(());
                }
            }
        }
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
