/*
use std::fmt::{Display, Formatter};
use anyhow::Error;
use crate::alignment::Direction::{Diag, Left, Up, Null};

pub mod pairwise_align;

#[derive(PartialEq, Clone, Copy)]
enum Direction {
    Diag,
    Left,
    Up,
    Null
}

#[derive(Clone)]
struct AlignConfig {
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
}

struct LocalAlignMatrix {
    query_len: usize,
    target_len: usize,
    score_matrix: Vec<Vec<i32>>,
    direction_matrix: Vec<Vec<Direction>>,
}

impl LocalAlignMatrix {
    fn new(query_len: usize, target_len: usize) -> Self {
        // target seq as col, query seq as row
        let score_row = vec![0i32; target_len + 1];
        let direction_row = vec![Null; target_len + 1];
        let mut score_matrix = vec![score_row; query_len + 1];
        let mut direction_matrix = vec![direction_row; query_len + 1];
        direction_matrix[0][0] = Diag;
        for col in 1..target_len + 1 {
            *(&mut direction_matrix[0][col]) = Left
        }
        for row in 1..query_len + 1 {
            *(&mut direction_matrix[row][0]) = Up
        }
        LocalAlignMatrix {
            query_len,
            target_len,
            score_matrix,
            direction_matrix,
        }
    }
}

struct GlobalAlignMatrix {
    
}

trait AlignMatrix {
    fn get_score(&self, row: usize, col: usize) -> i32;
    fn set_score(&mut self, row:usize, col:usize, score: i32) -> Result<(), anyhow::Error>;
    fn get_direction(&self, row:usize, col: usize) -> Direction;
    fn set_direction(&mut self, row: usize, col: usize, direction: Direction) -> Result<(),anyhow::Error>;
}

impl AlignMatrix  for LocalAlignMatrix{
    fn get_score(&self, row: usize, col: usize) -> i32 {
        self.score_matrix[row][col]
    }

    fn set_score(&mut self, row: usize, col: usize, score: i32) -> Result<(), Error> {
        self.score_matrix[row][col] =score;
        Ok(())
    }

    fn get_direction(&self, row: usize, col: usize) -> Direction {
        self.direction_matrix[row][col]
    }

    fn set_direction(&mut self, row: usize, col: usize, direction: Direction) -> Result<(), Error> {
        self.direction_matrix[row][col] = direction;
        Ok(())
    }
}

#[derive(Debug)]
struct AlignResult {
    align_target: Vec<u8>,
    align_line: Vec<u8>,
    align_query: Vec<u8>,
    max_score: i32,
    idx_of_max_score: (usize, usize), 
    idx_of_start: (usize, usize), // (query_row_idx, target_col_idx)
}

impl AlignResult {
    fn new() -> Self {
        AlignResult {
            align_target: Vec::<u8>::new(),
            align_line: Vec::<u8>::new(),
            align_query: Vec::<u8>::new(),
            max_score: 0i32,
            idx_of_max_score: (0usize, 0usize),
            idx_of_start: (0usize, 0usize),
        }
    }
    fn clear(&mut self) {
        self.align_target.clear();
        self.align_query.clear();
        self.align_line.clear();
        self.max_score = 0;
        self.idx_of_start = (0, 0);
        self.idx_of_max_score = (0, 0);
    }
    
    fn reverse(&mut self) {
        self.align_target.reverse();
        self.align_line.reverse();
        self.align_query.reverse();
    }
}

impl Display for AlignResult{
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let x = std::str::from_utf8(&self.align_target).unwrap();
        unsafe {
            write!(f, "target: {}\n        {}\nquery:  {}\n", std::str::from_utf8_unchecked(&self.align_target),
                   std::str::from_utf8_unchecked(&self.align_line),
                   std::str::from_utf8_unchecked(&self.align_query)
            )
        }
    }
}

 */