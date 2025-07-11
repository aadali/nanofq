pub mod adapter;

use crate::utils::{get_dege_bases, get_bases, IS_MATCHED};
use std::cell::Ref;
use crate::trim::adapter::{SequenceInfo, get_seq_info};
use adapter::End;
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::*;
use clap::builder::OsStringValueParser;
use seq_io::fastq::{Record, RefRecord};
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::iter::repeat_n;
use std::pin::pin;
use std::sync::OnceLock;
use std::time::Instant;
use flate2::Status::BufError;
use crate::fastq::NanoRead;

static DEGE_BASES: OnceLock<HashMap<u8, HashSet<u8>>> = OnceLock::new();
static BASES: OnceLock<HashMap<u8, u8>> = OnceLock::new();


struct TrimResult<'a> {
    id: &'a str,
    align_end5_end: Option<usize>,
    align_end3_start: Option<usize>,
    align_end5_str: Option<String>,
    align_end3_str: Option<String>,
}
impl<'a> TrimResult<'a> {
    fn is_empty(&self) -> bool {
        self.align_end5_end.is_none() && self.align_end3_start.is_none()
    }
    
    fn get_trim_log_and_trimmed_fastq(&self, ref_record: &RefRecord) -> (String, Option<String>) {
        let trimmed_record_string_opt = ref_record.sub_record_to_string(self.align_end5_end, self.align_end3_start);
        let mut log_string = String::new();
        log_string.push_str(&format!("{}: {}\nSearchIn5'End:\n", ref_record.id().unwrap(), ref_record.seq().len()));
        if self.align_end5_str.is_some() {
            log_string.push_str(self.align_end5_str.as_deref().unwrap());
        } else {
            log_string.push_str("NoTrimmed\n");
        }
        log_string.push_str("SearchIn3'End:\n");
        if self.align_end3_str.is_some() {
            log_string.push_str(self.align_end3_str.as_deref().unwrap());
        } else {
            log_string.push_str("NoTrimmed\n");
        }
        if trimmed_record_string_opt.is_none() {
            log_string.push_str(&format!("the full length of {} was trimmed\n", ref_record.id().unwrap()))
        }
        log_string.push_str(&(repeat_n('-', 100).collect::<String>() + "\n"));
        (log_string, trimmed_record_string_opt)
    }
}
impl<'a> Default for TrimResult<'a> {
    fn default() -> Self {
        TrimResult {
            id: "",
            align_end5_end: None,
            align_end3_start: None,
            align_end5_str: None,
            align_end3_str: None,
        }
    }
}

fn pct_and_ident_passed(end: &(&str, End), operations: &Vec<AlignmentOperation>) ->  bool {
    let match_len = operations.len();
    let ident_len = operations
        .iter()
        .filter(|x| *x == &AlignmentOperation::Match)
        .count();
    // todo!("the Match in bio::Alignment means equal, it's not for dege bases");
     match_len as f64 / end.0.len() as f64 > end.1.1 && ident_len as f64 / match_len as f64 > end.1.2
}

fn trim_end5<F: Fn(u8,u8) -> i32>(
    end: Option<(&str, End)>,
    read_seq: &[u8],
    aligner: &mut Aligner<F>,
) -> Option<(usize, String)> {
    if let Some(end5) = end {
        let end5_query = end5.0.as_bytes();
        let end5_cfg = end5.1;
        let end5_target = if read_seq.len() > end5_cfg.0 {
            &read_seq[..end5_cfg.0]
        } else {
            read_seq
        };
        let end5_alignment = aligner.local(end5_target, end5_query);
        if pct_and_ident_passed(&end5, &end5_alignment.operations){
            let align_end5_end = end5_alignment.xend;
            let align_end5_str = format!(
                "ReadAlignRange: {}, [{}, {})\nAdapterAlignRange: {}, [{}, {})\n{}\n",
                end5_target.len(),
                end5_alignment.xstart,
                end5_alignment.xend,
                end5_query.len(),
                end5_alignment.ystart,
                end5_alignment.yend,
                end5_alignment.pretty(end5_target, end5_query, 200)
                
            );
            Some((align_end5_end, align_end5_str))
        } else {
            None
        }
    } else {
        None
    }
}

fn trim_end3<F: Fn(u8,u8)->i32>(
    end: Option<(&str, End)>,
    read_seq: &[u8],
    aligner: &mut Aligner<F>,
) -> Option<(usize, String)> {
    if let Some(end3) = end {
        let end3_query = end3.0.as_bytes();
        let end3_cfg = end3.1;
        let end3_target = if read_seq.len() > end3_cfg.0 {
            &read_seq[read_seq.len() - end3_cfg.0..]
        } else {
            read_seq
        };
        let end3_alignment = aligner.local(end3_target, end3_query);
        println!("{}", end3_alignment.pretty(end3_target, end3_query, 200));
        if  pct_and_ident_passed(&end3, &end3_alignment.operations){
            let align_end3_start = read_seq.len() - end3_target.len() + end3_alignment.xstart;
            let align_end3_str = format!(
                "ReadAlignRange: {}, [{}, {})\nAdapterAlignRange: {}, [{}, {})\n{}\n",
                end3_target.len(),
                align_end3_start + 1,
                align_end3_start + end3_alignment.xend - end3_alignment.xstart,
                end3_query.len(),
                end3_alignment.ystart,
                end3_alignment.yend,
                end3_alignment.pretty(end3_target, end3_query, 200)
                
            );
            Some((align_end3_start, align_end3_str))
        } else {
            None
        }
    } else {
        None
    }
}

fn get_trim_result<'a>(
    trim_end5_opt: Option<(usize, String)>,
    trim_end3_opt: Option<(usize, String)>,
    ref_record: &'a RefRecord,
) -> TrimResult<'a> {
 if let Some(trim_fwd_end5) = trim_end5_opt {
        if let Some(trim_fwd_end3) = trim_end3_opt {
            TrimResult {
                id: ref_record.id().unwrap(),
                align_end5_end: Some(trim_fwd_end5.0),
                align_end3_start: Some(trim_fwd_end3.0),
                align_end5_str: Some(trim_fwd_end5.1),
                align_end3_str: Some(trim_fwd_end3.1),
            }
        } else {
            TrimResult {
                id: ref_record.id().unwrap(),
                align_end5_end: Some(trim_fwd_end5.0),
                align_end3_start: None,
                align_end5_str: Some(trim_fwd_end5.1),
                align_end3_str: None,
            }
        }
    } else {
        if let Some(trim_fwd_end3) = trim_end3_opt {
            TrimResult {
                id: ref_record.id().unwrap(),
                align_end5_end: None,
                align_end3_start: Some(trim_fwd_end3.0),
                align_end5_str: None,
                align_end3_str: Some(trim_fwd_end3.1),
            }
        } else {
            TrimResult {
                id: ref_record.id().unwrap(),
                align_end5_end: None,
                align_end3_start: None,
                align_end5_str: None,
                align_end3_str: None,
            }
        }
    }
    
}

fn trim_seq<'a, F: Fn(u8, u8) -> i32>(
    seq_info: &SequenceInfo,
    ref_record: &'a RefRecord,
    aligner: &mut Aligner<F>,
) -> (String, Option<String>) {
    let mut trim_result = TrimResult::default();
    let trim_fwd_end5_opt = trim_end5(seq_info.end5, ref_record.seq(), aligner );
    let trim_fwd_end3_opt = trim_end3(seq_info.end3, ref_record.seq(), aligner );
    
    trim_result = get_trim_result(trim_fwd_end5_opt, trim_fwd_end3_opt, ref_record);
    if trim_result.is_empty(){
        let trim_rev_end5_opt = trim_end5(
            seq_info.rev_com_end5,
            ref_record.seq(),
            aligner,
        );
        let trim_rev_end3_opt = trim_end3(
            seq_info.rev_com_end3,
            ref_record.seq(),
            aligner);
        trim_result = get_trim_result(trim_rev_end5_opt, trim_rev_end3_opt, ref_record);
        
    }
    trim_result.get_trim_log_and_trimmed_fastq(ref_record)
}


#[test]
pub fn test_local_align() {
    let start = Instant::now();
    // let align_config = AlignConfig {
    //     match_score: 3,
    //     mismatch_score: -3,
    //     gap_open_score: -7,
    //     gap_extend_score: -1,
    // };
    let query = b"CCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTT";
    // let target = b"TGTGCACGTACTTGGTTCGTTACGTATTGCTAAGGTTAACACAAAGACACGACAACTTTCTCAGCACCTGCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTCAATTTTTTTTTTGGAATTTTTTTTTTG";
    // let target = b"TTTGTTAGTCTACTTCGTTCAGTTACGTATTGCTAAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGCTTACGGTTCACTACTCACGACGATGTTTTTTTTGGTACCTTTTTTTCACCGGAGAGGACCCCGTAAAGTGATAATGATTATCATCTACATATCACAACGTGCGTGGAGGCCATCAAACCACGTCAAATAATCAATTATGACGCAGGTATCGTATTA";
    let target = b"GTTTTGTTTAACCTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCGCAGCACCTGTGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATCAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTCAAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAACCTAAGGAGCGTCACACTTGTCCTGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACAGTCAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTGGTGCAGTTGGTTGAAGCCACCACGCCCTTCCGGCGTGGATAATGCAGCCTCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCTCAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAAT";
    // let target = b"TTTAGCCTGTGCTTCGTTTAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCAGCACCTAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACGTCAAGCGGTCTGTAGTCGTGCGTAAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCGACTCGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATGGTATCTGTCAAGAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCAT";
    // let mut aligner = Aligner::new(-5, -1, MatchParams::new(3, -3));
    let mut aligner = Aligner::with_capacity(180, 70, -5, -1, MatchParams::new(3, -3));
    let alignment = aligner.local(target, query);
    let x = alignment.pretty(&target[..180], &query[..60], 200);
    // let y = alignment.pretty(b"ATCGGGTACAGATCAGAG", b"ATCGGGTACAGATCAGAGTGCAD", 100);
    let x = get_seq_info();
    let kit = x.get("LSK").unwrap();
    let x = kit.end5;
    let y = kit.end5;
    println!("target_start: {}", alignment.xstart);
    println!("target_end: {}", alignment.xend);
    println!("query_start: {}", alignment.ystart);
    println!("query_end:{}", alignment.yend);
    println!("query_len: {}", query.len());
    // println!("{}", x);
    println!("{:?}", alignment.operations);
}

// #[test]
pub fn test_trim() {
    use crate::fastq::FastqReader;
    use crate::fastq::NanoRead;
    // let mut reader = FastqReader::new(std::fs::File::open("/Users/aadali/test_data/nbd114.24/barcode01.fastq").unwrap());
    let mut reader = FastqReader::new(std::fs::File::open("/Users/aadali/test_data/pcb114.24/SRR30594249.fastq").unwrap());
    // let nbd_1 = get_seq_info()["NBD_1"];
    let pcb_1 = get_seq_info()["PCB"];
    let mut writer = BufWriter::new(std::fs::File::create("/Users/aadali/test_data/nbd114.24/trimmed.barcode01.fastq").unwrap());
    let mut log_writer = BufWriter::new(std::fs::File::create("/Users/aadali/test_data/nbd114.24/trimmed3.log").unwrap());
    log_writer.write(pcb_1.get_info().as_bytes()).unwrap();
    let score = |a:u8, b:u8| {if IS_MATCHED(&a, &b) {3} else {-3}};
    // let mut aligner = Aligner::with_capacity(200, 100, -5, -1, MatchParams {match_score:3, mismatch_score:-3});
    let mut aligner = Aligner::with_capacity(200, 100, -5, -1, score);
    // let mut aligner = Aligner::
    loop {
        let ref_record = reader.next();
        if ref_record.is_none() {
            break
        }
        let ref_record = ref_record.unwrap().unwrap();
        if ref_record.id().unwrap() != "SRR30594249.73" {
            continue
        }
        let (log, trimmed_fq) = trim_seq(&pcb_1, &ref_record, &mut aligner);
        log_writer.write(log.as_bytes()).unwrap();
        if trimmed_fq.is_some(){
            writer.write(trimmed_fq.unwrap().as_bytes()).unwrap();
        }
    }
}

#[test]
fn a() {
    let score = |a:u8, b:u8| {if IS_MATCHED(&b, &a ) {3} else {-3}};
    // let mut aligner = Aligner::with_capacity(200, 100, -5, -1, MatchParams {match_score:3, mismatch_score:-3});
    let mut aligner = Aligner::with_capacity(200, 200, -5, -1, score);
    let target ="CTGTGCATGATTATTTACTGGTTCAGTTATCCAGCCGATATTGCAGCCTGGCGCTGGCGCCGTTGACAAAGTTGTCGGTGTCTTTGTGACTTGCCTGCTCGCTCTCTTTCAGAGGAAGTCCGCCGCCCGCAAGTTTTTTTTTTTTTTTTTTTTTTTTTGT";
    let query = "AAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG";
    let alignment = aligner.local(target.as_bytes(), query.as_bytes());
    println!("{}", alignment.pretty(target.as_bytes(), query.as_bytes(), 200));
}