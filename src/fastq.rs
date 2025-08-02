use crate::alignment::LocalAligner;
use crate::trim::adapter::TrimConfig;
use crate::trim::trim_seq;
use crate::utils::get_q2p_table;
use ansi_term::Color;
use seq_io::fastq;
use seq_io::fastq::{Record, RefRecord};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::ops::{Deref, DerefMut};

const BUFF: usize = 1024 * 1024;

// (f64, f64): (this_read_average_error_pro, this_read_quality)
// (ReadID, Length, (ReadAverageErrProb, ReadQuality), Option<GCContent>)
pub type EachStats = (Box<String>, usize, (f64, f64), Option<f64>);

#[derive(Clone)]
pub(crate) struct FilterOption<'a> {
    pub min_len: usize,
    pub max_len: usize,
    pub min_qual: f64,
    pub max_qual: f64,
    pub gc: bool,
    pub min_gc: f64,
    pub max_gc: f64,
    pub retain_failed: Option<&'a String>,
}
impl<'a> FilterOption<'a> {
    pub(crate) fn set_failed_fastq_file(&self) -> Result<Option<BufWriter<File>>, anyhow::Error> {
        match self.retain_failed {
            None => Ok(Some(BufWriter::new(File::create(
                "/tmp/NanoFqFailed.fastq",
            )?))),
            Some(failed_fastq_file) => Ok(Some(BufWriter::new(File::create(failed_fastq_file)?))),
        }
    }
}
pub trait NanoRead {
    fn gc_count(&self) -> f64;
    fn calculate_read_quality(&self) -> (f64, f64);

    fn stats(&self, gc: bool) -> EachStats;

    fn is_passed(&self, fo: &FilterOption) -> bool;

    fn write(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error>;

    fn trim(
        &self,
        trim_cfg: &TrimConfig,
        aligner: &mut LocalAligner,
        min_len: usize,
        pretty_log: bool,
        trim_primer: bool,
    ) -> (Option<(&[u8], &[u8])>, Option<String>);
}
impl<'a, T> NanoRead for T
where
    T: Record,
{
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
            .map(|x| get_q2p_table()[*x as usize])
            .sum::<f64>()
            / seq_len;
        if avg_err_prob.is_finite() {
            let quality = avg_err_prob.log10() * -10.0;
            (avg_err_prob, quality)
        } else {
            (0.0, 0.0)
        }
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

    fn write(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
            write!(
                writer,
                "@{}\n{}\n+\n{}\n",
                unsafe {std::str::from_utf8_unchecked(self.head())},
                unsafe {std::str::from_utf8_unchecked(self.seq())},
                unsafe {std::str::from_utf8_unchecked(self.qual())}
            )?;
        Ok(())
    }

    fn trim(
        &self,
        trim_cfg: &TrimConfig,
        aligner: &mut LocalAligner,
        min_len: usize,
        pretty_log: bool,
        trim_primer: bool,
    ) -> (Option<(&[u8], &[u8])>, Option<String>) {
        let (trim_from, trim_to, log_string, _) = trim_seq(
            trim_cfg,
            self.seq(),
            &format!(
                "{}: {}",
                self.id().expect("parse into read id error"),
                self.seq().len()
            ),
            aligner,
            pretty_log,
            min_len,
            trim_primer,
        );
        if trim_from == 0 && trim_to == 0 {
            (None, log_string)
        } else {
            (
                Some((
                    &self.seq()[trim_from..trim_to],
                    &self.qual()[trim_from..trim_to],
                )),
                log_string,
            )
        }
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

    pub(crate) fn filter(
        &mut self,
        writer: &mut dyn Write,
        fo: &FilterOption,
        retain_failed: bool,
        failed_writer: &mut BufWriter<File>,
    ) -> Result<(), anyhow::Error> {
        if retain_failed {
            loop {
                match self.next() {
                    Some(ref_record_res) => {
                        let ref_record = ref_record_res.expect(
                            &Color::Red
                                .paint("Error: failed to get fastq record")
                                .to_string(),
                        );
                        if ref_record.is_passed(fo) {
                            NanoRead::write(&ref_record, writer)?
                        }
                    }
                    None => {
                        return Ok(());
                    }
                }
            }
        } else {
            loop {
                match self.next() {
                    Some(ref_record_res) => {
                        let ref_record = ref_record_res.expect(
                            &Color::Red
                                .paint("Error: failed to get fastq record")
                                .to_string(),
                        );
                        if ref_record.is_passed(fo) {
                            NanoRead::write(&ref_record, writer)?
                        } else {
                            NanoRead::write(&ref_record, failed_writer)?
                        }
                    }
                    None => {
                        return Ok(());
                    }
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
