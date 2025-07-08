use crate::alignment::{AlignConfig, AlignMatrix, AlignResult, Direction, GlobalAlignMatrix, LocalAlignMatrix};
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use std::time::Instant;


static DEGE_BASES: OnceLock<HashMap<u8, HashSet<u8>>> = OnceLock::new();
static BASES: OnceLock<HashMap<u8, u8>> = OnceLock::new();

fn get_dege_bases() -> &'static HashMap<u8, HashSet<u8>>{
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
            (b'n', b'N')
        ])
    })
}

#[inline]
fn is_matched(query_base: &u8, target_base: &u8) -> bool {
    query_base == target_base || {
        if let Some(target_bases) = get_dege_bases().get(query_base) {
            target_bases.contains(target_base)
        } else {
            false
        }
    }
}
pub fn smith_waterman(
    query: &[u8],
    target: &[u8],
    align_config: &AlignConfig,
    align_matrix: &mut LocalAlignMatrix,
    align_result: &mut AlignResult,
) -> Result<(), anyhow::Error> {
    let query = if align_matrix.query_len >= query.len() {
        query
    } else {
        &query[..align_matrix.query_len]
    };
    let target = if align_matrix.target_len >= target.len() {
        target
    } else {
        &target[..align_matrix.target_len]
    };

    for query_row in 1..query.len() + 1 {
        for target_col in 1..target.len() + 1 {
            // calculate score from diag
            let prev_score_from_diag = align_matrix.get_score(query_row - 1, target_col - 1);
            let bases_align_score = if is_matched(&query[query_row - 1], &target[target_col - 1]) {
                align_config.match_score
            } else {
                align_config.mismatch_score
            };
            let score_from_diag = prev_score_from_diag + bases_align_score;

            //calculate score from up
            let prev_score_from_up = align_matrix.get_score(query_row - 1, target_col);
            let gap_score_from_up =
                if align_matrix.get_direction(query_row - 1, target_col) == Direction::Up {
                    align_config.gap_extend_score
                } else {
                    align_config.gap_open_score
                };
            let score_from_up = prev_score_from_up + gap_score_from_up;

            // calculate score from left
            let prev_score_from_left = align_matrix.get_score(query_row, target_col - 1);
            let gap_score_from_left =
                if align_matrix.get_direction(query_row, target_col - 1) == Direction::Left {
                    align_config.gap_extend_score
                } else {
                    align_config.gap_open_score
                };
            let score_from_left = prev_score_from_left + gap_score_from_left;

            let this_max_score = *[0, score_from_diag, score_from_up, score_from_left]
                .iter()
                .max()
                .unwrap();

            // save this cell's max score
            align_matrix.set_score(query_row, target_col, this_max_score)?;

            // save the direction of this cell in direction_matrix
            if this_max_score == score_from_diag {
                align_matrix.set_direction(query_row, target_col, Direction::Diag)?;
            } else if this_max_score == score_from_up {
                align_matrix.set_direction(query_row, target_col, Direction::Up)?
            } else if this_max_score == score_from_left {
                align_matrix.set_direction(query_row, target_col, Direction::Left)?
            } else {
                // do nothing
            }
            align_matrix.set_score(query_row, target_col, this_max_score)?;

            // update the max score of the entire alignment and its row index and col index in matrix
            if this_max_score > align_result.max_score {
                align_result.idx_of_max_score = (query_row, target_col);
                align_result.max_score = this_max_score;
            }
        }
    }
    let (mut trace_query_row, mut trace_target_col) = align_result.idx_of_max_score;
    while trace_query_row > 0
        && trace_target_col > 0
        && align_matrix.get_score(trace_query_row, trace_target_col) > 0
    {
        let this_direction = align_matrix.get_direction(trace_query_row, trace_target_col);
        if this_direction == Direction::Diag {
            align_result.align_query.push(query[trace_query_row - 1]);
            align_result.align_target.push(target[trace_target_col - 1]);
            if is_matched(&query[trace_query_row - 1], &target[trace_target_col-1]) {
                align_result.align_line.push(b'|');
            } else {
                align_result.align_line.push(b':')
            }
            trace_target_col -= 1;
            trace_query_row -= 1;
        } else if this_direction == Direction::Up {
            align_result.align_target.push(b'-');
            align_result.align_query.push(query[trace_query_row - 1]);
            trace_query_row -= 1;
            align_result.align_line.push(b':');
        } else if this_direction == Direction::Left {
            align_result.align_query.push(b'-');
            align_result.align_target.push(target[trace_target_col - 1]);
            trace_target_col -= 1;
            align_result.align_line.push(b' ');
        } else {
        }
    }
    align_result.idx_of_start = (trace_query_row, trace_target_col);
    align_result.reverse();
    Ok(())
}

// global
fn needleman_wunsch(query: &[u8], target: &[u8], align_config: AlignConfig) {}

// #[cfg(test)]
#[test]
#[ignore]
pub fn test_local_align() {
    let start = Instant::now();
    let align_config = AlignConfig {
        match_score: 3,
        mismatch_score: -3,
        gap_open_score: -7,
        gap_extend_score: -1,
    };
    let query = b"CCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTT";
    let target = b"GGTAATGACATTTATTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTTATTGATTTCTACCATCTTCTACTCCGGCTTTTTTAGCAGCGAAGCGTTTGATAAGCGAACCAATCGAGTCAGTACCGATGTAGCCGATAAACACGCTCGTTATATAAGCGAGATTGCTACTTAGTCCGGCGAAGTCGAGAAGGTCACGA";
    // let target = b"TTTGTTAGTCTACTTCGTTCAGTTACGTATTGCTAAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGCTTACGGTTCACTACTCACGACGATGTTTTTTTTGGTACCTTTTTTTCACCGGAGAGGACCCCGTAAAGTGATAATGATTATCATCTACATATCACAACGTGCGTGGAGGCCATCAAACCACGTCAAATAATCAATTATGACGCAGGTATCGTATTA";
    // let target = b"GTTTTGTTTAACCTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCGCAGCACCTGTGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATCAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTCAAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAACCTAAGGAGCGTCACACTTGTCCTGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACAGTCAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTGGTGCAGTTGGTTGAAGCCACCACGCCCTTCCGGCGTGGATAATGCAGCCTCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCTCAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAAT";
    // let target = b"TTTAGCCTGTGCTTCGTTTAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCAGCACCTAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACGTCAAGCGGTCTGTAGTCGTGCGTAAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCGACTCGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATGGTATCTGTCAAGAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCAT";
    let mut local_matrix = LocalAlignMatrix::new(query.len(), target.len());
    let mut align_result = AlignResult::new();
    smith_waterman(&query[..], &target[..], &align_config, &mut local_matrix, &mut align_result).unwrap();
    println!("{}", &align_result);
    println!("{:?}", align_result);
    let dur = start.elapsed();
    println!("{:.10?}", dur)
}
