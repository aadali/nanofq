use crate::fastq::{DORADO_TRIM_LEADING_BASE_NUMBER, EachStats};
use crate::utils::{get_q2p_table, quit_with_error};
use bio::bio_types::genome::AbstractInterval;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar};
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::sync::OnceLock;

static ENCODED_BASES_GC_COUNT: OnceLock<HashMap<u8, usize>> = OnceLock::new();
pub fn get_encoded_bases_gc_count_table() -> &'static HashMap<u8, usize> {
    ENCODED_BASES_GC_COUNT.get_or_init(|| {
        HashMap::from([
            (18, 1),  // AC: 18
            (20, 0),  // AG: 20
            (24, 0),  // AT: 24
            (17, 0),  // AA: 17
            (129, 0), // TA: 129
            (136, 0), // TT: 136
            (132, 1), // TG: 132
            (130, 1), // TC: 130
            (33, 1),  // CA: 33
            (40, 1),  // CT: 40
            (36, 2),  // CG: 36
            (34, 2),  // CC: 34
            (65, 1),  // GA: 65
            (66, 2),  // GC: 66
            (68, 2),  // GG: 68
            (72, 1),  // GT: 72
            (16, 0),  // A : 16
            (128, 0), // T : 128
            (32, 1),  // C : 32
            (64, 1),  // G : 64
        ])
    })
}
pub trait BamRecordStats {
    fn gc_count(&self) -> f64;
    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> (f64, f64);
    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> EachStats;
}
impl BamRecordStats for rust_htslib::bam::Record {
    fn gc_count(&self) -> f64 {
        let seq_len = self.qual().len();
        let gc_number: usize = self
            .seq()
            .encoded
            .iter()
            .map(|x| *get_encoded_bases_gc_count_table().get(x).unwrap_or(&0usize))
            .sum();
        gc_number as f64 / seq_len as f64
    }

    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> (f64, f64) {
        let quals = self.qual();
        let real_seq_len = quals.len();
        if dont_use_dorado_quality {
            let avg_err_prob = quals
                .iter()
                .map(|x| get_q2p_table()[*x as usize + 33])
                .sum::<f64>()
                / real_seq_len as f64;
            if avg_err_prob.is_finite() {
                let read_quality = avg_err_prob.log10() * -10.0;
                (avg_err_prob, read_quality)
            } else {
                (0.0, 0.0)
            }
        } else {
            let quality_tag_res = self.aux(b"qs");
            if quality_tag_res.is_ok() {
                let quality_tag = quality_tag_res.unwrap();
                let read_quality = match quality_tag {
                    Aux::Float(quality) => quality as f64,
                    Aux::Double(quality) => quality,
                    _ => {
                        quit_with_error(&format!(
                            "Parse qs tag for {:?} in chr: {}, at position: {}",
                            str::from_utf8(self.qname()).unwrap(),
                            self.contig(),
                            self.reference_start(),
                        ));
                        0.0
                    }
                };
                (read_quality.powf(read_quality / -10.0), read_quality)
            } else {
                let (seq_len, skip) = if real_seq_len > DORADO_TRIM_LEADING_BASE_NUMBER {
                    (
                        real_seq_len - DORADO_TRIM_LEADING_BASE_NUMBER,
                        DORADO_TRIM_LEADING_BASE_NUMBER,
                    )
                } else {
                    (real_seq_len, 0)
                };
                let avg_err_prob = quals
                    .iter()
                    .skip(skip)
                    .map(|x| get_q2p_table()[*x as usize + 33])
                    .sum::<f64>()
                    / seq_len as f64;
                if avg_err_prob.is_finite() {
                    let quality = avg_err_prob.log10() * -10.0;
                    (avg_err_prob, quality)
                } else {
                    (0.0, 0.0)
                }
            }
        }
    }

    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> EachStats {
        let len = self.qual().len();
        let read_quality = self.calculate_read_quality(dont_use_dorado_quality);
        let gc = if gc { Some(self.gc_count()) } else { None };
        (
            Box::new(str::from_utf8(self.qname()).unwrap().to_string()),
            len,
            read_quality,
            gc,
        )
    }
}

#[derive(Default, Debug)]
pub struct BasicBamStatistics {
    // only for Aligned Bam
    reads_mapped: usize,            // primary alignments reads number
    reads_unmapped: usize,          // unmapped reads number
    reads_mq0: usize,               // primary alignment reads with map quality == 0
    primary_alignment: usize,       // primary alignments
    supplementary_alignment: usize, // supplementary alignments
    secondary_alignment: usize,     // secondary alignments
    bases_mapped: usize,            // primary_alignment seq_len()
    bases_mapped_cigar: usize, //  Match + Ins of (primary_alignment and supplementary_alignment), more accurate
    mismatches: usize,         // NM tag
    error_rate: f64,           // mismatches / bases_mapped_cigar
    map_rate: f64,             // reads_mapped / (reads_mapped + reads_unmapped)
}
impl BasicBamStatistics {
    pub fn is_empty(&self) -> bool {
        self.reads_mapped == 0
            && self.reads_unmapped == 0
            && self.reads_mq0 == 0
            && self.reads_mq0 == 0
            && self.primary_alignment == 0
            && self.supplementary_alignment == 0
            && self.secondary_alignment == 0
            && self.bases_mapped == 0
            && self.bases_mapped_cigar == 0
            && self.bases_mapped_cigar == 0
            && self.mismatches == 0
    }
}
impl Display for BasicBamStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "#BamStatsSummary\n\
            ReadsMapped\t{}\n\
            ReadsUnmapped\t{}\n\
            ReadsMapQ0\t{}\t# primary alignments with map quality==0\n\
            PrimaryAlignment\t{}\n\
            SupplementaryAlignment\t{}\n\
            SecondaryAlignment\t{}\n\
            BasesMapped\t{}\t# primary alignments seq_len, ignore clipping\n\
            BasesMappedCigar\t{}\t# Match + Ins of (primary alignments and supplementary alignments), more accurate\n\
            Mismatches\t{}\t# NM tag\n\
            ErrorRate\t{:.6}\t# Mismatches / BasesMappedCigar\n\
            MapRate\t{:.6}\t# ReadsMapped / (ReadsMapped + ReadsUnmapped)\n",
            self.reads_mapped,
            self.reads_unmapped,
            self.reads_mq0,
            self.primary_alignment,
            self.supplementary_alignment,
            self.secondary_alignment,
            self.bases_mapped,
            self.bases_mapped_cigar,
            self.mismatches,
            self.error_rate,
            self.map_rate
        )
    }
}

fn get_nm_aux(record: &bam::Record) -> usize {
    let nm_tag = record.aux(b"NM").unwrap_or(Aux::I32(0));
    match nm_tag {
        // NM tag in bam file is I32, but NM tag in sam file is U8, U16,
        // it maybe a bug in rust-htslib
        Aux::I8(nm) => nm as usize,
        Aux::U8(nm) => nm as usize,
        Aux::I16(nm) => nm as usize,
        Aux::U16(nm) => nm as usize,
        Aux::I32(nm) => nm as usize,
        Aux::U32(nm) => nm as usize,
        _ => {
            quit_with_error(&format!(
                "Parse NM tag failed for {}",
                str::from_utf8(record.qname()).unwrap()
            ));
            0usize
        }
    }
}

fn get_cigar_base_length(cigar: &Cigar) -> usize {
    match cigar {
        Cigar::Match(m) => *m as usize,
        Cigar::Ins(i) => *i as usize,
        _ => 0usize,
    }
}
pub fn stats_xam(
    bam_reader: &mut bam::Reader,
    thread: usize,
    gc: bool,
    dont_use_dorado_quality: bool,
) -> (BasicBamStatistics, Vec<EachStats>) {
    debug_assert!(thread > 0);
    let mut basic_bam_stats = BasicBamStatistics::default();
    let mut all_stats = vec![];
    if thread > 0 {
        bam_reader.set_threads(thread).unwrap();
    }
    let mut record = bam::Record::new();
    record.set_qname(b"InitRecord");
    while let Some(x) = bam_reader.read(&mut record) {
        if x.is_err() {
            quit_with_error(&format!(
                "Parse record failed: {}",
                str::from_utf8(record.qname()).unwrap()
            ));
        }
        if record.is_unmapped() {
            basic_bam_stats.reads_unmapped += 1;
            all_stats.push(record.stats(gc, dont_use_dorado_quality));
        } else {
            if record.flags() & 0x900 == 0 {
                all_stats.push(record.stats(gc, dont_use_dorado_quality));
                basic_bam_stats.bases_mapped += record.seq_len();
                basic_bam_stats.primary_alignment += 1;
                basic_bam_stats.reads_mapped += 1;
                if record.mapq() == 0 {
                    basic_bam_stats.reads_mq0 += 1;
                }
                let nm = get_nm_aux(&record);
                basic_bam_stats.mismatches += nm;
                for cigar in record.cigar().iter() {
                    basic_bam_stats.bases_mapped_cigar += get_cigar_base_length(cigar);
                }
            } else {
                if record.is_supplementary() {
                    basic_bam_stats.supplementary_alignment += 1;
                    let nm = get_nm_aux(&record);
                    basic_bam_stats.mismatches += nm;
                    for cigar in record.cigar().iter() {
                        basic_bam_stats.bases_mapped_cigar += get_cigar_base_length(cigar);
                    }
                }
                if record.is_secondary() {
                    basic_bam_stats.secondary_alignment += 1;
                }
            }
        }
    }
    basic_bam_stats.error_rate =
        basic_bam_stats.mismatches as f64 / basic_bam_stats.bases_mapped_cigar as f64;
    basic_bam_stats.map_rate = basic_bam_stats.reads_mapped as f64
        / (basic_bam_stats.reads_mapped + basic_bam_stats.reads_unmapped) as f64;
    (basic_bam_stats, all_stats)
}
