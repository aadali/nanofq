use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

static DEGE_BASES: OnceLock<HashMap<u8, HashSet<u8>>> = OnceLock::new();
static BASES: OnceLock<HashMap<u8, u8>> = OnceLock::new();
static Q2P_TABLE: OnceLock<[f64; 128]> = OnceLock::new();

pub fn get_dege_bases() -> &'static HashMap<u8, HashSet<u8>> {
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

pub fn get_bases() -> &'static HashMap<u8, u8> {
    BASES.get_or_init(|| {
        HashMap::from([
            (b'A', b'T'),
            (b'T', b'A'),
            (b'G', b'C'),
            (b'C', b'G'),
            (b'R', b'Y'),
            (b'Y', b'R'),
            (b'M', b'K'),
            (b'K', b'M'),
            (b'S', b'S'),
            (b'W', b'W'),
            (b'H', b'D'),
            (b'B', b'V'),
            (b'V', b'B'),
            (b'D', b'H'),
            (b'a', b'T'),
            (b't', b'A'),
            (b'g', b'C'),
            (b'c', b'G'),
            (b'N', b'N'),
            (b'n', b'N'),
        ])
    })
}

pub fn get_q2p_table() -> &'static [f64; 128] {
    Q2P_TABLE.get_or_init(|| {
        let mut arr = [f64::NAN; 128];
        for q in 33..127usize {
            arr[q] = 10.0f64.powf((q - 33) as f64 / -10.0)
        }
        arr
    })
}

pub fn rev_com(seq: &str) -> String {
    seq.as_bytes()
        .iter()
        .map(|x| *get_bases().get(x).unwrap() as char)
        .rev()
        .collect::<String>()
}
// ref_base from primer or reference can be degenerate base
pub static IS_MATCHED: fn(&u8, &u8) -> bool = |ref_base, read_base| {
    ref_base == read_base
        || get_dege_bases()
            .get(ref_base)
            .map_or(false, |x| x.contains(read_base))
};

pub const SEP_LINE: &str =
    "----------------------------------------------------------------------\n";

#[test]
fn test_dege_base() {
    assert!(IS_MATCHED(&b'V', &b'A'));
    assert!(IS_MATCHED(&b'A', &b'A'));
    assert!(!IS_MATCHED(&b'C', &b'A'));
    // assert!(IS_MATCHED(&b'G', &b'V'));
    assert!(IS_MATCHED(&b'B', &b'C'));
    assert!(IS_MATCHED(&b'B', &b'T'));
    assert!(IS_MATCHED(&b'B', &b'G'));
    assert!(IS_MATCHED(&b'W', &b'T'));
}
