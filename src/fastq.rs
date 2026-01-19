use crate::alignment::LocalAligner;
use crate::trim::adapter::TrimConfig;
use crate::trim::trim_seq;
use crate::utils::get_q2p_table;
use ansi_term::Color;
use seq_io::fastq;
use seq_io::fastq::Record;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::ops::{Deref, DerefMut};

const BUFF: usize = 1024 * 1024;
pub const DORADO_TRIM_LEADING_BASE_NUMBER: usize = 60;

// (ReadID, Length,  ReadQuality, Option<GCContent>)
pub type EachStats = (String, u32, f32, Option<f32>);

#[derive(Clone)]
pub(crate) struct FilterOption<'a> {
    pub min_len: u32,
    pub max_len: u32,
    pub dont_use_dorado_quality: bool,
    pub min_qual: f32,
    pub max_qual: f32,
    pub gc: bool,
    pub min_gc: f32,
    pub max_gc: f32,
    pub retain_failed: Option<&'a String>,
}
impl<'a> FilterOption<'a> {
    pub(crate) fn set_failed_fastq_file(&self) -> Result<Option<Box<dyn Write>>, anyhow::Error> {
        match self.retain_failed {
            None => Ok(Some(Box::new(std::io::sink()))),
            Some(failed_fastq_file) => Ok(Some(Box::new(BufWriter::new(File::create(
                failed_fastq_file,
            )?)))),
        }
    }
}
pub trait NanoRead {
    fn gc_count(&self) -> Option<f32>;
    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> Option<f32>;

    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> Option<EachStats>;

    fn is_passed(&self, fo: &FilterOption) -> (bool, String, u32, f32);

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
    fn gc_count(&self) -> Option<f32> {
        let seq_len = self.seq().len();
        if seq_len == 0 {
            eprintln!(
                "Empty sequence found for: {:?}",
                self.id().unwrap_or("UnknowReadName")
            );
            return None;
        }
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
        Some(gc_number as f32 / seq_len as f32)
    }

    #[inline]
    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> Option<f32> {
        let qual_real_len = self.qual().len();
        if qual_real_len == 0 {
            eprintln!(
                "Empty quality found for: {:?}",
                self.id().unwrap_or("UnknowReadName")
            );
            return None;
        }
        let trim_leading =
            (!dont_use_dorado_quality) && qual_real_len > DORADO_TRIM_LEADING_BASE_NUMBER;
        let seq_len = if trim_leading {
            qual_real_len - DORADO_TRIM_LEADING_BASE_NUMBER
        } else {
            qual_real_len
        };
        let avg_err_prob = self
            .qual()
            .iter()
            .skip(if trim_leading {
                DORADO_TRIM_LEADING_BASE_NUMBER
            } else {
                0
            })
            .map(|x| get_q2p_table()[*x as usize])
            .sum::<f64>()
            / seq_len as f64;
        let quality = avg_err_prob.log10() * -10.0;
        Some(quality as f32)
    }

    #[inline]
    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> Option<EachStats> {
        let len = self.seq().len();
        if len == 0 {
            return None;
        }
        let read_quality_opt = self.calculate_read_quality(dont_use_dorado_quality);
        let gc = if gc { self.gc_count() } else { None };
        match read_quality_opt {
            Some(read_quality) => Some((
                str::from_utf8(
                    self.head()
                        .split(|b| *b == b' ' || *b == b'\t') // the header of dorado basecaller output fastq separated by tab
                        .next()
                        .unwrap(),
                )
                .expect(&Color::Red.paint(format!(
                    "parse read id failed in record header: {}",
                    str::from_utf8(self.head()).unwrap_or("unknow header")
                )))
                .to_string(),
                len as u32,
                read_quality,
                gc,
            )),
            None => None,
        }
    }

    #[inline]
    fn is_passed(&self, fo: &FilterOption) -> (bool, String, u32, f32) {
        let seq_len = self.seq().len() as u32;
        let gc_passed = if fo.gc {
            let gc_opt = self.gc_count();
            if gc_opt.is_some() {
                let gc = gc_opt.unwrap();
                gc > fo.min_gc && gc < fo.max_gc
            } else {
                false
            }
        } else {
            true
        };
        let this_read_qual = self
            .calculate_read_quality(fo.dont_use_dorado_quality)
            .unwrap_or(0.0);
        let is_passed = seq_len >= fo.min_len
            && seq_len <= fo.max_len
            && this_read_qual > fo.min_qual
            && this_read_qual < fo.max_qual
            && gc_passed;
        (
            is_passed,
            str::from_utf8(
                self.head()
                    .split(|b| *b == b' ' || *b == b'\t')
                    .next()
                    .unwrap(),
            )
            .expect(&Color::Red.paint(format!(
                "parse read id failed in record header: {}",
                str::from_utf8(self.head()).unwrap_or("unknow header")
            )))
            .to_string(),
            seq_len,
            this_read_qual,
        )
    }

    fn write(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
        write!(
            writer,
            "@{}\n{}\n+\n{}\n",
            unsafe { std::str::from_utf8_unchecked(self.head()) },
            unsafe { std::str::from_utf8_unchecked(self.seq()) },
            unsafe { std::str::from_utf8_unchecked(self.qual()) }
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
    pub(crate) fn stats(&mut self, gc: bool, dont_use_dorado_quality: bool) -> Vec<EachStats> {
        let mut stats_result = vec![];
        loop {
            let ref_record_opt = self.next();
            if let Some(ref_record) = ref_record_opt {
                let each_stats = ref_record
                    .expect(
                        &Color::Red
                            .paint("Error: failed to get fastq record")
                            .to_string(),
                    )
                    .stats(gc, dont_use_dorado_quality);
                if each_stats.is_some() {
                    stats_result.push(each_stats.unwrap());
                }
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
        failed_writer: &mut dyn Write,
    ) -> Result<Vec<(String, u32, f32)>, anyhow::Error> {
        let mut stats_results = Vec::new();
        if retain_failed {
            loop {
                match self.next() {
                    Some(ref_record_res) => {
                        let ref_record = ref_record_res.expect(
                            &Color::Red
                                .paint("Error: failed to get fastq record")
                                .to_string(),
                        );
                        let read_filter_res = ref_record.is_passed(fo);
                        if read_filter_res.0 {
                            stats_results.push((
                                read_filter_res.1,
                                read_filter_res.2,
                                read_filter_res.3,
                            ));
                            NanoRead::write(&ref_record, writer)?
                        }
                    }
                    None => {
                        return Ok(stats_results);
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
                        let read_filter_res = ref_record.is_passed(fo);
                        if read_filter_res.0 {
                            stats_results.push((
                                read_filter_res.1,
                                read_filter_res.2,
                                read_filter_res.3,
                            ));
                            NanoRead::write(&ref_record, writer)?
                        } else {
                            NanoRead::write(&ref_record, failed_writer)?
                        }
                    }
                    None => {
                        return Ok(stats_results);
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
