use crate::alignment::{AlignConfig, AlignMatrix, AlignResult, Direction, GlobalAlignMatrix, LocalAlignMatrix};
use lazy_static::lazy_static;
use std::collections::{HashMap, HashSet};

// local
lazy_static! {
    static ref BASES: HashMap<u8, HashSet<u8>> = {
        let mut bases = HashMap::new();
        bases.insert('R' as u8, HashSet::from(['A' as u8, 'G' as u8]));
        bases.insert('Y' as u8, HashSet::from(['C' as u8, 'T' as u8]));
        bases.insert('M' as u8, HashSet::from(['C' as u8, 'A' as u8]));
        bases.insert('K' as u8, HashSet::from(['G' as u8, 'T' as u8]));
        bases.insert('S' as u8, HashSet::from(['C' as u8, 'G' as u8]));
        bases.insert('W' as u8, HashSet::from(['A' as u8, 'T' as u8]));
        bases.insert('H' as u8, HashSet::from(['A' as u8, 'T' as u8, 'C' as u8]));
        bases.insert('B' as u8, HashSet::from(['G' as u8, 'T' as u8, 'C' as u8]));
        bases.insert('V' as u8, HashSet::from(['G' as u8, 'A' as u8, 'C' as u8]));
        bases.insert('D' as u8, HashSet::from(['G' as u8, 'A' as u8, 'T' as u8]));
        bases.insert(
            'N' as u8,
            HashSet::from(['G' as u8, 'A' as u8, 'T' as u8, 'C' as u8]),
        );
        bases
    };
}
#[inline]
fn is_match(query_base: &u8, target_base: &u8) -> bool {
    query_base == target_base || {
        if let Some(target_bases) = BASES.get(query_base) {
            target_bases.contains(target_base)
        } else {
            false
        }
    }
}
fn smith_waterman(
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
            let bases_align_score = if is_match(&query[query_row - 1], &target[target_col - 1]) {
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
    // trace back depend on the direction_matrix
    let (mut trace_query_row, mut trace_target_col) = align_result.idx_of_max_score;
        // (query_row_of_max_score, target_col_of_max_score);
    // let mut align_query: Vec<u8> = vec![];
    // let mut align_target: Vec<u8> = vec![];
    // let mut line: Vec<u8> = vec![];
    while trace_query_row > 0
        && trace_target_col > 0
        && align_matrix.get_score(trace_query_row, trace_target_col) > 0
    {
        let this_direction = align_matrix.get_direction(trace_query_row, trace_target_col);
        if this_direction == Direction::Diag {
            align_result.align_query.push(query[trace_query_row - 1]);
            align_result.align_target.push(target[trace_target_col - 1]);
            if is_match(&query[trace_query_row - 1], &target[trace_target_col-1]) {
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
            align_result.align_line.push(b' ');
        } else if this_direction == Direction::Left {
            align_result.align_query.push(b'-');
            align_result.align_target.push(target[trace_target_col - 1]);
            trace_target_col -= 1;
            align_result.align_line.push(b' ');
        } else {
        }
    }
    align_result.idx_of_start = (trace_query_row, trace_target_col);
    // println!("{trace_query_row}-{query_row_of_max_score}");
    // println!("{trace_target_col}-{target_col_of_max_score}");
    align_result.reverse();
    println!("{}", &align_result);
    // align_target.reverse();
    // align_query.reverse();
    // line.reverse();
    // println!("target: {}", std::str::from_utf8(&align_target)?);
    // println!("line:   {}", std::str::from_utf8(&line)?);
    // println!("query:  {}", std::str::from_utf8(&align_query)?);
    Ok(())
}

// global
fn needleman_wunsch(query: &[u8], target: &[u8], align_config: AlignConfig) {}

// #[cfg(test)]
#[test]
pub fn test_local_align() {
    let align_config = AlignConfig {
        match_score: 3,
        mismatch_score: -3,
        gap_open_score: -7,
        gap_extend_score: -1,
    };
    let query = b"CCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTT";
    // let target = b"GGTAATGACATTTATTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTTATTGATTTCTACCATCTTCTACTCCGGCTTTTTTAGCAGCGAAGCGTTTGATAAGCGAACCAATCGAGTCAGTACCGATGTAGCCGATAAACACGCTCGTTATATAAGCGAGATTGCTACTTAGTCCGGCGAAGTCGAGAAGGTCACGA";
    // let target = b"TTTGTTAGTCTACTTCGTTCAGTTACGTATTGCTAAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGCTTACGGTTCACTACTCACGACGATGTTTTTTTTGGTACCTTTTTTTCACCGGAGAGGACCCCGTAAAGTGATAATGATTATCATCTACATATCACAACGTGCGTGGAGGCCATCAAACCACGTCAAATAATCAATTATGACGCAGGTATCGTATTA";
    // let target = b"GTTTTGTTTAACCTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCGCAGCACCTGTGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATCAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTCAAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAACCTAAGGAGCGTCACACTTGTCCTGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACAGTCAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTGGTGCAGTTGGTTGAAGCCACCACGCCCTTCCGGCGTGGATAATGCAGCCTCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCTCAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAAT";
    let target = b"TTTAGCCTGTGCTTCGTTTAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTCAGCACCTAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACGTCAAGCGGTCTGTAGTCGTGCGTAAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCGACTCGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATGGTATCTGTCAAGAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCAT";
    let mut local_matrix = LocalAlignMatrix::new(query.len(), target.len());
    let mut align_result = AlignResult::new();
    smith_waterman(&query[..], &target[..], &align_config, &mut local_matrix, &mut align_result).unwrap();
    println!("{}", &align_result);
}
