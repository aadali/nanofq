pub mod adapter;

use crate::alignment::{LocalAligner, LocalAlignment, ReadEnd};
use crate::trim::adapter::{EndConfig, SequenceInfo};
use crate::utils::SEP_LINE;
use seq_io::fastq::{Record};
use std::io::{Write};

fn trim_end<'a>(
    end_config: &'a EndConfig,
    read_seq: &'a [u8],
    aligner: &mut LocalAligner,
    end: ReadEnd,
) -> Option<(usize, usize, LocalAlignment<'a>)> {
    if let Some(end_config) = end_config {
        let end_reference = end_config.0.as_bytes();
        let end_align_para = end_config.1;
        let read_end_seq = if read_seq.len() > end_align_para.0 {
            if end == ReadEnd::End5 {
                &read_seq[..end_align_para.0]
            } else {
                &read_seq[read_seq.len() - end_align_para.0..]
            }
        } else {
            read_seq
        };
        let mut alignment = aligner.align(end_reference, read_end_seq);
        let (ident, ident_pct) = alignment.get_ident();
        let align_pct = alignment.get_percent();
        if align_pct > end_align_para.1 && ident_pct > end_align_para.2 {
            Some((read_end_seq.len(), ident, alignment))
        } else {
            None
        }
    } else {
        None
    }
}

pub fn trim_seq(
    seq_info: &SequenceInfo,
    seq: &[u8],
    id: &str,
    aligner: &mut LocalAligner,
    log: bool,
) -> (usize, usize, Option<String>) {
    let read_seq = seq;
    let mut fwd_trim_from = 0;
    let mut fwd_trim_to = read_seq.len();
    let mut end3_used_len = 0;
    let mut end5_alignment = LocalAlignment::default();
    let mut end3_alignment = LocalAlignment::default();
    let mut trim_end5_success = false;
    let mut trim_end3_success = false;
    let mut fwd_ident_score = 0;
    let mut pretty_log = if log {
        Some(format!("{}\n", id,))
    } else {
        None
    };
    // actually, the forward end5 must be used to search, this means seq_info.end5.is_some() must be true. The following expr must be true
    if seq_info.may_trim_end5() {
        // Step1. consider to align end5
        if let Some((_, end5_ident, end5_align)) =
            trim_end(&seq_info.end5, read_seq, aligner, ReadEnd::End5)
        {
            end5_alignment = end5_align;
            fwd_trim_from = end5_alignment.read_range.0 - 1;
            trim_end5_success = true;
            fwd_ident_score += end5_ident;
        }
    }
    if seq_info.may_trim_end3() {
        // Step2. consider to align end3
        if let Some((end3_len, end3_ident, end3_align)) =
            trim_end(&seq_info.end3, read_seq, aligner, ReadEnd::End3)
        {
            end3_used_len = end3_len;
            end3_alignment = end3_align;
            fwd_trim_to = read_seq.len() - end3_used_len + end3_alignment.read_range.0 - 1;
            trim_end3_success = true;
            fwd_ident_score += end3_ident;
        }
    }
    // if rev_com_end5 is used, so the rev_com_end3 must be used as well
    debug_assert_eq!(
        seq_info.may_trim_rev_com_end5(),
        seq_info.may_trim_rev_com_end3(),
        "rev_com end5 and rev_com end3 must be fit"
    );
    /*
    For some kits, such as RAD/NBD/RBK/ULK/LSK, only forward read would be checked. For these kits, I don't care the rev com read.

    For other kits, such as PCS/PCB, the both ends of forward and rev com reads would be checked. For these kits, the rev com read
    will be checked if no adapter was found in both ends of forward simultaneously [CASE1] or just in one end(end5 or end3)[CASE2] of forward.

    if CASE1 we thought the right alignments already be found and just use the trim_from and trim_to index of forward read to trim original sequence.

    else if CASE2 the rev com read will be checked as well. And the total alignment identity bases number (ident base number in end5 + ident base number in end3)
    will be calculated for forward read and rev com read. More identity bases number, Better the alignment are.
    Finally, we will use the trim_from and trim_to index of read (forward or rev com) that has more identity base
    to trim the original sequence
     */
    if !seq_info.may_trim_rev_com_end5() {
        // Step3. if for this seq_info, rev_com align is not needed, then just use trim info from Step1 and Step2
        if trim_end5_success {
            pretty_log
                .as_mut()
                .map(|x| x.push_str(&end5_alignment.pretty(ReadEnd::End5)));
        }
        if trim_end3_success {
            pretty_log
                .as_mut()
                .map(|x| x.push_str(&end3_alignment.pretty(ReadEnd::End3)));
        }
        if fwd_trim_to >= fwd_trim_from {
            pretty_log.as_mut().map(|x| x.push_str(SEP_LINE));
            return (fwd_trim_from, fwd_trim_to, pretty_log);
        } else {
            pretty_log
                .as_mut()
                .map(|x| x.push_str("Full seq was trimmed, drop it\n"));
            return (0, 0, pretty_log); // if the original sequence is too short, maybe the align start of end3 is less than the align end of end5
        }
    } else {
        // Step4. if the rev_com read should be also detected
        if trim_end5_success && trim_end3_success {
            // Step5. if the check of both ends of forward passed, then just use trim info from Step1 and Step2
            pretty_log
                .as_mut()
                .map(|x| x.push_str(&end5_alignment.pretty(ReadEnd::End5)));
            pretty_log
                .as_mut()
                .map(|x| x.push_str(&end3_alignment.pretty(ReadEnd::End3)));
            if fwd_trim_to >= fwd_trim_from {
                pretty_log.as_mut().map(|x| x.push_str(SEP_LINE));
                return (fwd_trim_from, fwd_trim_to, pretty_log);
            } else {
                pretty_log
                    .as_mut()
                    .map(|x| x.push_str("Full seq was trimmed, drop it\n"));
                return (0, 0, pretty_log); // if the original sequence is too short, maybe the align start of end3 is less than the align end of end5
            }
        } else {
            // Step6. if just one end of forward passed, then consider the both ends of rev_com and do Step7
            let mut rev_ident_score = 0;
            let mut rev_trim_from = 0;
            let mut trim_rev_com_end5_success = false;
            let mut trim_rev_com_end3_success = false;
            let mut rev_trim_to = read_seq.len();
            let mut rev_com_end5_alignment = LocalAlignment::default();
            let mut rev_com_end3_alignment = LocalAlignment::default();
            // Step7. check end5 of rev_com
            if let Some((_, rev_com_end5_ident, rev_com_end5_align)) =
                trim_end(&seq_info.rev_com_end5, read_seq, aligner, ReadEnd::End5)
            {
                rev_com_end5_alignment = rev_com_end5_align;
                rev_trim_from = rev_com_end5_alignment.read_range.0 - 1;
                rev_ident_score += rev_com_end5_ident;
                trim_rev_com_end5_success = true;
            }
            // Step8. check end3 of rev_com
            if let Some((rev_com_end3_len, rev_com_end3_ident, rev_com_end3_align)) =
                trim_end(&seq_info.rev_com_end3, read_seq, aligner, ReadEnd::End3)
            {
                end3_used_len = rev_com_end3_len;
                rev_com_end3_alignment = rev_com_end3_align;
                rev_trim_to =
                    read_seq.len() - end3_used_len + rev_com_end3_alignment.read_range.0 - 1;
                rev_ident_score += rev_com_end3_ident;
                trim_rev_com_end3_success = true;
            }
            // Step9. determine which read (forward or rev_com) will be used depends on the total identity bases number of each direction
            if fwd_ident_score > rev_ident_score {
                if fwd_trim_to >= fwd_trim_from {
                    // Step10. if identity bases numbers of forward is more, just use trim info from Step1 and Step2
                    if trim_end5_success {
                        pretty_log
                            .as_mut()
                            .map(|x| x.push_str(&end5_alignment.pretty(ReadEnd::End5)));
                    }
                    if trim_end3_success {
                        pretty_log
                            .as_mut()
                            .map(|x| x.push_str(&end3_alignment.pretty(ReadEnd::End3)));
                    }
                    pretty_log.as_mut().map(|x| x.push_str(SEP_LINE));
                    return (fwd_trim_from, fwd_trim_to, pretty_log);
                } else {
                    pretty_log
                        .as_mut()
                        .map(|x| x.push_str("Full seq was trimmed, drop it\n"));
                    return (0, 0, pretty_log);
                }
            } else {
                // Step11. if identity bases numbers of rev_com is more, just use trim info from Step7 and Step8
                if rev_trim_to >= rev_trim_from {
                    if trim_rev_com_end5_success {
                        pretty_log
                            .as_mut()
                            .map(|x| x.push_str(&rev_com_end5_alignment.pretty(ReadEnd::End5)));
                    }
                    if trim_rev_com_end3_success {
                        pretty_log
                            .as_mut()
                            .map(|x| x.push_str(&rev_com_end3_alignment.pretty(ReadEnd::End3)));
                    }
                    pretty_log.as_mut().map(|x| x.push_str(SEP_LINE));
                    return (rev_trim_from, rev_trim_to, pretty_log);
                } else {
                    pretty_log
                        .as_mut()
                        .map(|x| x.push_str("Full seq was trimmed, drop it\n"));
                    return (0, 0, pretty_log);
                }
            }
        }
    }
}

/*
#[test]
pub fn test_trim() {
    use crate::fastq::FastqReader;
    use crate::fastq::NanoRead;
    let mut reader = FastqReader::new(std::fs::File::open("/Users/aadali/test_data/nbd114.24/barcode01.fastq").unwrap());
    // let mut reader = FastqReader::new(
    //     // std::fs::File::open("/Users/aadali/test_data/pcb114.24/SRR30594249.fastq").unwrap(),
    //     std::fs::File::open("/Users/aadali/test_data/pcb114.24/sub_SRR30594249.fastq").unwrap(),
    // );
    // let nbd_1 = get_seq_info()["NBD_1"];
    // let pcb_1 = get_seq_info()["PCB"];
    let nbd_1= get_seq_info()["NBD_1"];
    // let mut writer = BufWriter::new(
    //     // std::fs::File::create("/Users/aadali/test_data/nbd114.24/trimmed.barcode01.fastq").unwrap(),
    //     std::fs::File::create("/Users/aadali/test_data/pcb114.24/trimmed.barcode01.fastq").unwrap(),
    // );
    let mut log_writer = BufWriter::new(
        std::fs::File::create("/Users/aadali/test_data/nbd114.24/nbd1_trimmed.log").unwrap(),
        // std::fs::File::create("/Users/aadali/test_data/pcb114.24/trimmed2.log").unwrap(),
    );
    // log_writer.write(pcb_1.get_info().as_bytes()).unwrap();
    log_writer.write(nbd_1.get_info().as_bytes()).unwrap();
    // let mut aligner = Aligner::with_capacity(200, 100, -5, -1, MatchParams {match_score:3, mismatch_score:-3});
    // let mut aligner = Aligner::with_capacity(200, 100, -5, -1, score);
    let mut aligner = LocalAligner::new(
        90,
        180,
        Scores {
            match_: 3,
            mismatch: -3,
            gap_open: -5,
            gap_extend: -1,
        },
    );
    // let mut aligner = Aligner::
    loop {
        let ref_record = reader.next();
        if ref_record.is_none() {
            break;
        }
        let ref_record = ref_record.unwrap().unwrap();
        // if ref_record.id().unwrap() != "SRR30594249.10" {
        //     continue;
        // }
         let (trim_from, trim_to, pretty_log) =trim_seq(nbd_1, ref_record.seq(), ref_record.id().unwrap(), &mut aligner, true);
        {
            // println!("{:?}", pretty_log);
            if pretty_log.is_some() {
                log_writer.write(pretty_log.unwrap().as_bytes()).unwrap();
            }
        }
        // let x = trim_seq(pcb_1, &ref_record, &mut aligner);
        // log_writer.write(log.as_bytes()).unwrap();
        // if trimmed_fq.is_some() {
        //     writer.write(trimmed_fq.unwrap().as_bytes()).unwrap();
        // }
    }
}

#[test]
fn a() {
    use crate::utils::IS_MATCHED;
    use bio::alignment::pairwise::Aligner;
    let score = |a: u8, b: u8| {
        if IS_MATCHED(&b, &a) { 3 } else { -3 }
    };
    // let mut aligner = Aligner::with_capacity(200, 100, -5, -1, MatchParams {match_score:3, mismatch_score:-3});
    let mut aligner = Aligner::with_capacity(200, 200, -5, -1, score);
    let target = "CTGTGCATGATTATTTACTGGTTCAGTTATCCAGCCGATATTGCAGCCTGGCGCTGGCGCCGTTGACAAAGTTGTCGGTGTCTTTGTGACTTGCCTGCTCGCTCTCTTTCAGAGGAAGTCCGCCGCCCGCAAGTTTTTTTTTTTTTTTTTTTTTTTTTGT";
    let query = "AAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG";
    let alignment = aligner.local(target.as_bytes(), query.as_bytes());
    println!(
        "{}",
        alignment.pretty(target.as_bytes(), query.as_bytes(), 200)
    );
}


 */
