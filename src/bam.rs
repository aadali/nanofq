use crate::fastq::{DORADO_TRIM_LEADING_BASE_NUMBER, EachStats};
use crate::utils::{get_q2p_table, quit_with_error};
use bio::bio_types::genome::AbstractInterval;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
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
        let trim_leading =
            (!dont_use_dorado_quality) && quals.len() > DORADO_TRIM_LEADING_BASE_NUMBER;
        let (seq_len, skip) = if trim_leading {
            (
                real_seq_len - DORADO_TRIM_LEADING_BASE_NUMBER,
                DORADO_TRIM_LEADING_BASE_NUMBER,
            )
        } else {
            (real_seq_len, 0)
        };
        let quality_tag_res = self.aux(b"qs");
        let read_quality = if quality_tag_res.is_ok() {
            let quality_tag = quality_tag_res.unwrap();
            match quality_tag {
                Aux::Float(quality) => quality as f64,
                _ => {
                    quit_with_error(&format!(
                        "error qs tag for {:?} in chr: {}, at position: {}",
                        str::from_utf8(self.qname()).unwrap(),
                        self.contig(),
                        self.reference_start(),
                    ));
                    0.0
                }
            }
        } else {
            let avg_err_prob = quals
                .iter()
                .skip(skip)
                .map(|x| get_q2p_table()[*x as usize + 33])
                .sum::<f64>()
                / seq_len as f64;
            if avg_err_prob.is_finite() {
                let quality = avg_err_prob.log10() * -10.0;
                quality
            } else {
                0.0
            }
        };
        (0.0, read_quality)
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
