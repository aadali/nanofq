pub mod adapter;
use crate::trim::adapter::SequenceInfo;
use adapter::End;
use bio::alignment::AlignmentOperation;
use bio::alignment::pairwise::*;
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use std::time::Instant;

static DEGE_BASES: OnceLock<HashMap<u8, HashSet<u8>>> = OnceLock::new();
static BASES: OnceLock<HashMap<u8, u8>> = OnceLock::new();

fn get_dege_bases() -> &'static HashMap<u8, HashSet<u8>> {
    DEGE_BASES.get_or_init(|| {
        HashMap::from([
            (b'R', HashSet::from([b'A', b'G'])),
            (b'Y', HashSet::from([b'C', b'T'])),
            (b'M', HashSet::from([b'C', b'A'])),
            (b'K', HashSet::from([b'G', b'T'])),
            (b'S', HashSet::from([b'C', b'G'])),
            (b'W', HashSet::from([b'A', b'T'])),
            (b'H', HashSet::from([b'A', b'T', b'C'])),
            (b'B', HashSet::from([b'G', b'T', b'C'])),
            (b'V', HashSet::from([b'G', b'A', b'C'])),
            (b'D', HashSet::from([b'G', b'A', b'T'])),
            (b'N', HashSet::from([b'G', b'A', b'T', b'C'])),
        ])
    })
}

fn get_bases() -> &'static HashMap<u8, u8> {
    BASES.get_or_init(|| {
        HashMap::from([
            (b'A', b'T'),
            (b'T', b'A'),
            (b'G', b'C'),
            (b'C', b'G'),
            (b'a', b'T'),
            (b't', b'A'),
            (b'g', b'C'),
            (b'c', b'G'),
            (b'N', b'N'),
            (b'n', b'N'),
        ])
    })
}

#[inline]
fn is_matched(
    query_base: &u8,
    target_base: &u8,
    dege_bases: &'static HashMap<u8, HashSet<u8>>,
) -> bool {
    query_base == target_base || {
        if let Some(target_bases) = dege_bases.get(query_base) {
            target_bases.contains(target_base)
        } else {
            false
        }
    }
}

struct QueryInfo {
    query_seq: Vec<u8>,
}

fn trim(seq_info: &SequenceInfo, target: &[u8], aligner: &mut Aligner<MatchParams>) {
    let end5_alignment = aligner.local(seq_info.end5.0.as_bytes(), target);
    ()
}

fn pct_and_ident_passed(end: &End, operations: &Vec<AlignmentOperation>) -> bool {
    let match_len = operations.len();
    let ident_len = operations
        .iter()
        .filter(|x| *x == &AlignmentOperation::Match)
        .count();
    match_len as f64 / end.0 as f64 > end.1 && ident_len as f64 / match_len as f64 > end.2
}

fn trim_end5_with_forward(
    seq_info: &SequenceInfo,
    target: &[u8],
    aligner: &mut Aligner<MatchParams>,
    output_pretty: bool,
) -> (Option<usize>, Option<String>) {
    let end5 = seq_info.end5;
    let alignment = aligner.local(target, end5.0.as_bytes());
    if pct_and_ident_passed(&end5.1, &alignment.operations) {
        if output_pretty {
            let mut align_str = String::new();
            align_str.push_str(&format!(
                "ReadAlignRange: ({}, {})\nAdapterAlignRange: ({}, {})\n{}\n",
                alignment.xstart,
                alignment.xend,
                alignment.ystart,
                alignment.yend,
                alignment.pretty(target, end5.0.as_bytes(), 200)
            ));
            (Some(alignment.xend), Some(align_str))
        } else {
            (Some(alignment.xend), None)
        }
    } else {
        (None, None)
    }
}

fn trim_end3_with_forward(
    seq_info: &SequenceInfo,
    target: &[u8],
    aligner: &mut Aligner<MatchParams>,
    output_pretty: bool,
) -> (Option<usize>, Option<String>) {
    let end3 = seq_info.end3.unwrap();
    let alignment = aligner.local(target, end3.0.as_bytes());
    if pct_and_ident_passed(&end3.1, &alignment.operations) {
        if output_pretty {
            let align_str = format!(
                "ReadAlignRange: (-{}, -{})\nAdapterAlignRange: ({}, {})\n{}\n",
                target.len() - alignment.xstart + 1,
                target.len() - alignment.xend + 1,
                alignment.ystart,
                alignment.yend,
                alignment.pretty(target, end3.0.as_bytes(), 200)
            );
            (Some(target.len() - alignment.xstart + 1), Some(align_str))
        } else {
            (Some(target.len() - alignment.xstart + 1), None)
        }
    } else {
        (None, None)
    }
}

fn trim_end5_with_rev_com(
    seq_info: &SequenceInfo,
    target: &[u8],
    aligner: &mut Aligner<MatchParams>,
    output_pretty: bool,
) -> (Option<usize>, Option<String>) {
    let rev_com_end5 = seq_info.rev_com_end5.unwrap();
    let alignment = aligner.local(target, rev_com_end5.0.as_bytes());
    if pct_and_ident_passed(&rev_com_end5.1, &alignment.operations) {
        if output_pretty {
            let align_str = format!(
                "ReadAlignRange: ({}, {})\nAdapterAlignRange: ({}, {})\n{}\n",
                alignment.xstart,
                alignment.xend,
                alignment.ystart,
                alignment.yend,
                alignment.pretty(target, rev_com_end5.0.as_bytes(), 200)
            );
            (Some(alignment.xend), Some(align_str))
        } else {
            (Some(alignment.xend), None)
        }
    } else {
        (None, None)
    }
}

fn trim_end3_with_rev_com(
    seq_info: &SequenceInfo,
    target: &[u8],
    aligner: &mut Aligner<MatchParams>,
    output_pretty: bool,
) -> (Option<usize>, Option<String>) {
    let rev_com_end3 = seq_info.rev_com_end3.unwrap();
    let alignment = aligner.local(target, rev_com_end3.0.as_bytes());
    if pct_and_ident_passed(&rev_com_end3.1, &alignment.operations) {
        if output_pretty {
            let align_str = format!(
                "ReadAlignRange: (-{}, -{})\nAdapterAlignRange: ({}, {})\n{}\n",
                target.len() - alignment.xstart + 1,
                target.len() - alignment.xend + 1,
                alignment.ystart,
                alignment.yend,
                alignment.pretty(target, rev_com_end3.0.as_bytes(), 200)
            );
            (Some(target.len() - alignment.xstart + 1), Some(align_str))
        } else {
            (Some(target.len() - alignment.xstart + 1), None)
        }
    } else {
        (None, None)
    }
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
    let target = b"TGTGCACGTACTTGGTTCGTTACGTATTGCTAAGGTTAACACAAAGACACGACAACTTTCTCAGCACCTGCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTCAATTTTTTTTTTGGAATTTTTTTTTTG";
    // let target = b"TTTGTTAGTCTACTTCGTTCAGTTACGTATTGCTAAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGCTTACGGTTCACTACTCACGACGATGTTTTTTTTGGTACCTTTTTTTCACCGGAGAGGACCCCGTAAAGTGATAATGATTATCATCTACATATCACAACGTGCGTGGAGGCCATCAAACCACGTCAAATAATCAATTATGACGCAGGTATCGTATTA";
    // let target = b"GTTTTGTTTAACCTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCGCAGCACCTGTGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATCAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTCAAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAACCTAAGGAGCGTCACACTTGTCCTGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACAGTCAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTGGTGCAGTTGGTTGAAGCCACCACGCCCTTCCGGCGTGGATAATGCAGCCTCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCTCAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAAT";
    // let target = b"TTTAGCCTGTGCTTCGTTTAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCAGCACCTAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACGTCAAGCGGTCTGTAGTCGTGCGTAAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCGACTCGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATGGTATCTGTCAAGAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCAT";
    // let mut aligner = Aligner::new(-5, -1, MatchParams::new(3, -3));
    let mut aligner =
        Aligner::with_capacity(query.len(), target.len(), -5, -1, MatchParams::new(3, -3));
    let alignment = aligner.local(target, query);
    let x = alignment.pretty(target, query, 200);
    // let y = alignment.pretty(b"ATCGGGTACAGATCAGAG", b"ATCGGGTACAGATCAGAGTGCAD", 100);
    println!("target_start: {}", alignment.xstart);
    println!("target_end: {}", alignment.xend);
    println!("query_start: {}", alignment.ystart);
    println!("query_end:{}", alignment.yend);
    println!("query_len: {}", query.len());
    println!("{}", x);
    println!("{:?}", alignment.operations);
}
