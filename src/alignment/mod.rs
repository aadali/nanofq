use crate::utils::IS_MATCHED;

static REF: &[u8; 0] = b"";
static READ: &[u8; 0] = b"";
#[derive(Debug, Clone, PartialEq, Copy)]
enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins,
    Null,
}

#[derive(PartialEq)]
pub enum ReadEnd {
    End5,
    End3,
}

#[derive(Debug, Clone, Default, Copy)]
pub struct Scores {
    pub(crate) match_: i32,
    pub(crate) mismatch: i32,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
}

#[derive(Debug)]
pub struct LocalAlignment<'a> {
    // all index in this struct used 1-based coordinate system
    reference: &'a [u8],
    read_end: &'a [u8],
    pub ref_range: (usize, usize),
    pub read_range: (usize, usize),
    max_score: i32,
    read_map_ref_operations: Vec<AlignmentOperation>,
}

impl<'a> LocalAlignment<'a> {
    #[inline]
    pub fn get_ident(&self) -> (usize, f64) {
        let ident = self
            .read_map_ref_operations
            .iter()
            .filter(|x| *x == &AlignmentOperation::Match)
            .count();
        (
            ident,
            ident as f64 / self.read_map_ref_operations.len() as f64,
        )
    }

    #[inline]
    pub fn get_percent(&self) -> f64 {
        (self.ref_range.1 - self.ref_range.0 + 1) as f64 / self.reference.len() as f64
    }

    pub fn pretty(&self, end: ReadEnd) -> String {
        let mut align_ref_vec = Vec::<u8>::new();
        let mut line_vec = Vec::<u8>::new();
        let mut align_read_vec = Vec::<u8>::new();
        let mut ref_idx = self.ref_range.0 - 1;
        let mut read_idx = self.read_range.0 - 1;
        for operation in &self.read_map_ref_operations {
            match operation {
                &AlignmentOperation::Match => {
                    align_ref_vec.push(self.reference[ref_idx]);
                    line_vec.push(b'|');
                    align_read_vec.push(self.read_end[read_idx]);
                    ref_idx += 1;
                    read_idx += 1;
                }
                &AlignmentOperation::Subst => {
                    align_ref_vec.push(self.reference[ref_idx]);
                    line_vec.push(b'\\');
                    align_read_vec.push(self.read_end[read_idx]);
                    ref_idx += 1;
                    read_idx += 1;
                }
                &AlignmentOperation::Del => {
                    align_read_vec.push(b'-');
                    line_vec.push(b'x');
                    align_ref_vec.push(self.reference[ref_idx]);
                    ref_idx += 1;
                }
                &AlignmentOperation::Ins => {
                    align_read_vec.push(self.read_end[read_idx]);
                    line_vec.push(b'x');
                    align_ref_vec.push(b'-');
                    read_idx += 1;
                }
                _ => {
                    panic!("you could not in here")
                }
            }
        }
        debug_assert_eq!(ref_idx, self.ref_range.1);
        debug_assert_eq!(read_idx, self.read_range.1);
        let mut read_range_start = self.read_range.0 as isize;
        let mut read_range_end = self.read_range.1 as isize;
        let align_read_str = match end {
            ReadEnd::End5 => {
                format!(
                    "{:>5} {} {:<5}",
                    self.read_range.0,
                    unsafe { std::str::from_utf8_unchecked(&align_read_vec) },
                    self.read_range.1
                )
            }
            ReadEnd::End3 => {
                read_range_start = 0 - (self.read_end.len() as isize - read_range_start + 1);
                read_range_end = 0 - (self.read_end.len() as isize - read_range_end + 1);
                format!(
                    "{:>5} {} {:<5}",
                    read_range_start,
                    unsafe { std::str::from_utf8_unchecked(&align_read_vec) },
                    read_range_end
                )
            }
        };
        let line = format!("      {}      ", unsafe {
            std::str::from_utf8_unchecked(&line_vec)
        });
        let align_ref_string = format!(
            "{:>5} {} {:<5}",
            self.ref_range.0,
            unsafe { std::str::from_utf8_unchecked(&align_ref_vec) },
            self.ref_range.1
        );
        let which_end = if end == ReadEnd::End5 {
            "Align5'End:\n"
        } else {
            "Align3'End:\n"
        };
        // format!(
        //     "{}AlignPercent: {}\nIdentPercent:{}\nReadAlignRange:{}, [{}, {}]\nRefAlignRange:{}, [{}, {}]\n{}\n{}\n{}\n",
        //     which_end,
        //     self.get_percent(),
        //     self.get_ident().1,
        //     self.read_end.len(),
        //     read_range_start,
        //     read_range_end,
        //     self.reference.len(),
        //     self.ref_range.0,
        //     self.ref_range.1,
        //     align_read_str,
        //     line,
        //     align_ref_string
        // )
        format!(
            "{}ReadAlignRange:{}, [{}, {}]\nRefAlignRange:{}, [{}, {}]\n{}\n{}\n{}\n",
            which_end,
            self.read_end.len(),
            read_range_start,
            read_range_end,
            self.reference.len(),
            self.ref_range.0,
            self.ref_range.1,
            align_read_str,
            line,
            align_ref_string
        )
    }
}

impl<'a> Default for LocalAlignment<'a> {
    fn default() -> Self {
        LocalAlignment {
            reference: &REF[..],
            read_end: &READ[..],
            ref_range: (1, 1),
            read_range: (1, 1),
            max_score: i32::MIN,
            read_map_ref_operations: vec![],
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct LocalAligner {
    row: usize,
    col: usize,
    scores: Scores,
    score_matrix: Vec<Vec<i32>>,
    align_matrix: Vec<Vec<AlignmentOperation>>,
}

impl LocalAligner {
    pub fn new(row: usize, col: usize, scores: Scores) -> Self {
        debug_assert!(scores.gap_extend < 0, "gap_extend must be negative int");
        debug_assert!(scores.gap_open < 0, "gap_open must be negative int");
        debug_assert!(scores.match_ > 0, "match must be positive int");
        debug_assert!(scores.mismatch < 0, "mismatch must be negative int");
        let score_matrix_row = vec![0i32; col + 1];
        let score_matrix = vec![score_matrix_row; row + 1];
        let mut align_matrix = vec![vec![AlignmentOperation::Null; col + 1]; row + 1];
        for i in 0..col + 1 {
            align_matrix[0][i] = AlignmentOperation::Ins;
        }
        for j in 0..row + 1 {
            align_matrix[j][0] = AlignmentOperation::Del;
        }
        LocalAligner {
            row,
            col,
            scores,
            score_matrix,
            align_matrix,
        }
    }

    pub fn update(&mut self, row: usize, col: usize, scores: Scores) {
        let tmp = Self::new(row, col, scores);
        self.row = tmp.row;
        self.col = tmp.col;
        self.scores = tmp.scores;
        self.score_matrix = tmp.score_matrix;
        self.align_matrix = tmp.align_matrix;
    }

    pub fn align<'a>(&mut self, reference: &'a [u8], read: &'a [u8]) -> LocalAlignment<'a> {
        debug_assert!(self.row >= reference.len());
        debug_assert!(self.col >= read.len());
        let mut local_alignment = LocalAlignment::default();
        local_alignment.reference = reference;
        local_alignment.read_end = read;
        for ref_row in 1..(reference.len() + 1) {
            for read_col in 1..(read.len() + 1) {
                let prev_score_from_diag = self.score_matrix[ref_row - 1][read_col - 1];
                let prev_score_from_up = self.score_matrix[ref_row - 1][read_col];
                let prev_score_from_left = self.score_matrix[ref_row][read_col - 1];
                let (this_score_from_diag, diag_bases_operation) =
                    if IS_MATCHED(&reference[ref_row - 1], &read[read_col - 1]) {
                        (
                            self.scores.match_ + prev_score_from_diag,
                            AlignmentOperation::Match,
                        )
                    } else {
                        (
                            self.scores.mismatch + prev_score_from_diag,
                            AlignmentOperation::Subst,
                        )
                    };
                let this_score_from_up =
                    if self.align_matrix[ref_row - 1][read_col] == AlignmentOperation::Del {
                        prev_score_from_up + self.scores.gap_extend
                    } else {
                        prev_score_from_up + self.scores.gap_open
                    };
                let this_score_from_left =
                    if self.align_matrix[ref_row][read_col - 1] == AlignmentOperation::Ins {
                        prev_score_from_left + self.scores.gap_extend
                    } else {
                        prev_score_from_left + self.scores.gap_open
                    };
                let this_cell_score = *[
                    this_score_from_diag,
                    this_score_from_left,
                    this_score_from_up,
                    0,
                ]
                .iter()
                .max()
                .unwrap();
                if this_cell_score > local_alignment.max_score {
                    local_alignment.max_score = this_cell_score;
                    local_alignment.ref_range.1 = ref_row;
                    local_alignment.read_range.1 = read_col;
                }
                self.score_matrix[ref_row][read_col] = this_cell_score;
                if this_cell_score == this_score_from_diag {
                    self.align_matrix[ref_row][read_col] = diag_bases_operation;
                } else if this_cell_score == this_score_from_left {
                    self.align_matrix[ref_row][read_col] = AlignmentOperation::Ins;
                } else {
                    self.align_matrix[ref_row][read_col] = AlignmentOperation::Del
                }
            }
        }

        let mut trace_ref_row = local_alignment.ref_range.1;
        let mut trace_read_col = local_alignment.read_range.1;
        // println!(
        //     "{}-{}",
        //     local_alignment.ref_range.1, local_alignment.read_range.1
        // );
        loop {
            let this_operation = self.align_matrix[trace_ref_row][trace_read_col];
            local_alignment.read_map_ref_operations.push(this_operation);
            if this_operation == AlignmentOperation::Match
                || this_operation == AlignmentOperation::Subst
            {
                trace_ref_row -= 1;
                trace_read_col -= 1;
                if !(trace_read_col > 0
                    && trace_read_col > 0
                    && self.score_matrix[trace_ref_row][trace_read_col] > 0)
                {
                    local_alignment.ref_range.0 = trace_ref_row + 1;
                    local_alignment.read_range.0 = trace_read_col + 1;
                    break;
                }
            } else if this_operation == AlignmentOperation::Del {
                trace_ref_row -= 1;
                if !(trace_read_col > 0
                    && trace_read_col > 0
                    && self.score_matrix[trace_ref_row][trace_read_col] > 0)
                {
                    local_alignment.ref_range.0 = trace_ref_row + 1;
                    local_alignment.read_range.0 = trace_read_col;
                    break;
                }
            } else {
                trace_read_col -= 1;
                if !(trace_read_col > 0
                    && trace_read_col > 0
                    && self.score_matrix[trace_ref_row][trace_read_col] > 0)
                {
                    local_alignment.ref_range.0 = trace_ref_row;
                    local_alignment.read_range.0 = trace_read_col + 1;
                    break;
                }
            }
        }
        local_alignment.read_map_ref_operations.reverse();
        local_alignment
    }
}

#[test]
fn test_max() {
    let v = [1, 23, 32];
    let a = 20;
    let b = 234;
    let c = 90;
    let x = *[a, b, c].iter().max().unwrap();
    println!("{x}");
    let x = *v.iter().max().unwrap();
    println!("{x}");
}

pub fn test_alignment() {
    let scores = Scores {
        match_: 3,
        mismatch: -3,
        gap_open: -5,
        gap_extend: -1,
    };
    let mut aligner = LocalAligner::new(80, 200, scores);
    let reference = b"AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT";
    // let read = b"TTTGTTAGTCTACTTCGTTCAGTTACGTATTGCTAAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTGCTTACGGTTCACTACTCACGACGAT";
    let read = b"GTTAACCTACTTCGTTCAGTTACGTATTGCTAAGGTTAACACAAAGACACCGACAACTTCTTCAGCACCTGCTTACGGTTCACTACTCACGAC";
    // let read = b"ATGGA";
    // let reference = b"ATGCGA";
    let alignment = aligner.align(reference, read);
    println!("{:?}", alignment);
    let x = alignment.pretty(ReadEnd::End5);
    println!("{x}");
}
