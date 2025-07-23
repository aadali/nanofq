pub mod adapter;

use crate::alignment::{LocalAligner, LocalAlignment, ReadEnd};
use crate::trim::adapter::{EndConfig, TrimConfig};
use crate::utils::SEP_LINE;

fn trim_end<'a>(
    end_cfg: &'a EndConfig,
    read_seq: &'a [u8],
    aligner: &mut LocalAligner,
    end: ReadEnd,
) -> Option<(usize, usize, LocalAlignment<'a>)> {
    if let Some(end_config) = end_cfg {
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
        let alignment = aligner.align(end_reference, read_end_seq);
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

fn get_trim_return(
    trim_from: usize,
    trim_to: usize,
    min_len: usize,
    log: &mut Option<String>,
) -> (usize, usize) {
    if trim_to > trim_from {
        if trim_to - trim_from > min_len {
            log.as_mut().map(|x| x.push_str(SEP_LINE));
            (trim_from, trim_to)
        } else {
            log.as_mut().map(|x| {
                x.push_str(&format!(
                    "{} is too short after trimming, drop it\n{}",
                    trim_to - trim_from,
                    SEP_LINE
                ))
            });
            (0, 0)
        }
    } else {
        log.as_mut()
            .map(|x| x.push_str(&format!("Full length is trimmed, drop it\n{}", SEP_LINE)));
        (0, 0)
    }
}
pub fn trim_seq(
    trim_cfg: &TrimConfig,
    seq: &[u8],
    id: &str,
    aligner: &mut LocalAligner,
    log: bool,
    min_len: usize,
    trim_primer: bool,
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
    // actually, the forward end5 will always be used to search, this means trim_cfg.end5.is_some() must be true. The following expr always be true
    if trim_cfg.may_trim_end5() {
        // Step1. consider to align end5
        if let Some((_, end5_ident, end5_align)) =
            trim_end(&trim_cfg.end5, read_seq, aligner, ReadEnd::End5)
        {
            end5_alignment = end5_align;
            fwd_trim_from = if trim_primer {end5_alignment.read_range.1} else {end5_alignment.read_range.0};
            trim_end5_success = true;
            fwd_ident_score += end5_ident;
        }
    }
    if trim_cfg.may_trim_end3() {
        // Step2. consider to align end3
        if let Some((end3_len, end3_ident, end3_align)) =
            trim_end(&trim_cfg.end3, read_seq, aligner, ReadEnd::End3)
        {
            end3_used_len = end3_len;
            end3_alignment = end3_align;
            fwd_trim_to = if trim_primer {
                read_seq.len() - end3_used_len + end3_alignment.read_range.0 - 1
            } else {
                read_seq.len() - end3_used_len + end3_alignment.read_range.1 - 1
            };
            trim_end3_success = true;
            fwd_ident_score += end3_ident;
        }
    }
    // if rev_com_end5 is used, so the rev_com_end3 must be used as well
    debug_assert_eq!(
        trim_cfg.may_trim_rev_com_end5(),
        trim_cfg.may_trim_rev_com_end3(),
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
    if !trim_cfg.may_trim_rev_com_end5() {
        // Step3. if for this trim_cfg, rev_com align is not needed, then just use trim info from Step1 and Step2
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
        let (a, b) = get_trim_return(fwd_trim_from, fwd_trim_to, min_len, &mut pretty_log);
        (a, b, pretty_log)
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
            let (a, b) = get_trim_return(fwd_trim_from, fwd_trim_to, min_len, &mut pretty_log);
            (a, b, pretty_log)
        } else {
            // Step6. if just one end of forward passed, then consider the both ends of rev_com and do Step7
            let mut rev_ident_score = 0;
            let mut rev_trim_from = 0;
            let mut rev_trim_to = read_seq.len();
            let mut trim_rev_com_end5_success = false;
            let mut trim_rev_com_end3_success = false;
            let mut rev_com_end5_alignment = LocalAlignment::default();
            let mut rev_com_end3_alignment = LocalAlignment::default();
            // Step7. check end5 of rev_com
            if let Some((_, rev_com_end5_ident, rev_com_end5_align)) =
                trim_end(&trim_cfg.rev_com_end5, read_seq, aligner, ReadEnd::End5)
            {
                rev_com_end5_alignment = rev_com_end5_align;
                rev_trim_from = if trim_primer {rev_com_end5_alignment.read_range.1} else {rev_com_end5_alignment.read_range.0};
                rev_ident_score += rev_com_end5_ident;
                trim_rev_com_end5_success = true;
            }
            // Step8. check end3 of rev_com
            if let Some((rev_com_end3_len, rev_com_end3_ident, rev_com_end3_align)) =
                trim_end(&trim_cfg.rev_com_end3, read_seq, aligner, ReadEnd::End3)
            {
                end3_used_len = rev_com_end3_len;
                rev_com_end3_alignment = rev_com_end3_align;
                rev_trim_to = if trim_primer {
                    read_seq.len() - end3_used_len + rev_com_end3_alignment.read_range.0 - 1
                } else {
                    read_seq.len() - end3_used_len + rev_com_end3_alignment.read_range.1 - 1
                };
                rev_ident_score += rev_com_end3_ident;
                trim_rev_com_end3_success = true;
            }
            // Step9. determine which read (forward or rev_com) will be used depends on the total identity bases number of each direction
            if fwd_ident_score > rev_ident_score {
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
                // Step10. if identity bases numbers of forward is more, just use trim info from Step1 and Step2
                let (a, b) = get_trim_return(fwd_trim_from, fwd_trim_to, min_len, &mut pretty_log);
                (a, b, pretty_log)
            } else {
                // Step11. if identity bases numbers of rev_com is more, just use trim info from Step7 and Step8
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
                let (a, b) = get_trim_return(rev_trim_from, rev_trim_to, min_len, &mut pretty_log);
                (a, b, pretty_log)
            }
        }
    }
}
