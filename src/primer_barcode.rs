use crate::utils::quit_with_error;
use ahash::{HashMap, RandomState};
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::myers::Myers;
use std::fmt::{Display, Formatter};
use std::hash::Hash;

#[derive(PartialEq)]
pub enum PO {
    // Primer Orient
    F, // Forward
    R, // Reverse
}

#[derive(Hash, Clone)]
pub struct Primer {
    pub name: String,
    fwd: Vec<u8>,
    rev: Vec<u8>,
    fwd_rc: Vec<u8>,
    rev_rc: Vec<u8>,
}

impl Primer {
    pub fn new(name: &str, fwd: &[u8], rev: &[u8]) -> Self {
        Primer {
            name: name.to_string(),
            fwd: fwd.to_vec(),
            rev: rev.to_vec(),
            fwd_rc: revcomp(fwd),
            rev_rc: revcomp(rev),
        }
    }

    pub fn parse_primer_from_cli(primers_arg: &str) -> HashMap<String, Primer> {
        let mut primers = HashMap::with_hasher(RandomState::new());
        for each_field in primers_arg.split(";") {
            let mut primer_fields = each_field.split(",");
            let name = primer_fields.next();
            let fwd_primer = primer_fields.next();
            let rev_primer = primer_fields.next();
            if name.is_none() || fwd_primer.is_none() || rev_primer.is_none() {
                quit_with_error(
                    "Failed to parse primers from cli. Primer argument should be <PrimerName>,<Froward>,<Reverse>[;<PrimerName>,<Froward>,<Reverse>]...",
                )
            }
            let name = name.unwrap();
            let fwd_primer = fwd_primer.unwrap().as_bytes();
            let rev_primr = rev_primer.unwrap().as_bytes();
            let primer = Primer::new(name, fwd_primer, rev_primr);
            if primers.insert(name.to_string(), primer).is_some() {
                quit_with_error(&format!("Duplicate primer name found: {name}"))
            }
        }
        primers
    }

    pub fn parse_primer_from_file(primer_file: &str) -> HashMap<String, Primer> {
        let mut primers = HashMap::with_hasher(RandomState::new());
        let primers_content = std::fs::read_to_string(primer_file)
            .expect(&format!("Failed to read primer file: {}", primer_file));
        for (line_number, line) in primers_content.lines().enumerate() {
            if line.starts_with("#") {
                continue;
            }
            let mut primer_fields = line.split("\t");
            let name = primer_fields.next();
            let fwd_primer = primer_fields.next();
            let rev_primer = primer_fields.next();
            if name.is_none() || fwd_primer.is_none() || rev_primer.is_none() {
                quit_with_error(&format!(
                    "Failed to parse primer at line {}",
                    line_number + 1
                ))
            }
            let name = name.unwrap();
            let fwd_primer = fwd_primer.unwrap().as_bytes();
            let rev_primr = rev_primer.unwrap().as_bytes();
            let primer = Primer::new(name, fwd_primer, rev_primr);
            if primers.insert(name.to_string(), primer).is_some() {
                quit_with_error(&format!("Duplicate primer name found: {name}"))
            }
        }
        primers
    }

    pub fn primer_seq2_primer_name(
        primers: &HashMap<String, Primer>,
    ) -> HashMap<&[u8], (PO, &str)> {
        let mut primers_seq =
            HashMap::with_capacity_and_hasher(primers.len() * 2, RandomState::new());
        let primer_min_len = primers
            .iter()
            .map(|each| [each.1.fwd.len(), each.1.rev.len()])
            .flatten()
            .min()
            .expect("Failed to get primer with min length ");
        for (primer_name, primer) in primers {
            primers_seq.insert(&primer.fwd[..primer_min_len], (PO::F, primer_name.as_ref()));
            primers_seq.insert(&primer.rev[..primer_min_len], (PO::R, primer_name.as_ref()));
        }
        primers_seq
    }

    pub fn fwd_myers(&self) -> Myers {
        Myers::new(&self.fwd)
    }

    pub fn rev_myers(&self) -> Myers {
        Myers::new(&self.rev)
    }

    pub fn fwd_rc_myers(&self) -> Myers {
        Myers::new(&self.fwd_rc)
    }

    pub fn rev_rc_myers(&self) -> Myers {
        Myers::new(&self.rev_rc)
    }
}

pub fn get_myers_from_primers(primers: &HashMap<String, Primer>) -> HashMap<String, [Myers; 4]> {
    primers
        .iter()
        .map(|(primer_name, primer)| {
            (
                primer_name.clone(),
                [
                    primer.fwd_myers(),
                    primer.rev_rc_myers(),
                    primer.rev_myers(),
                    primer.fwd_rc_myers(),
                ],
            )
        })
        .collect::<HashMap<_, _>>()
}

impl Display for Primer {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Forward: {}\tForwardRc: {}\nReverse: {}\tReverseRc: {}",
            str::from_utf8(&self.fwd).unwrap(),
            str::from_utf8(&self.fwd_rc).unwrap(),
            str::from_utf8(&self.rev).unwrap(),
            str::from_utf8(&self.rev_rc).unwrap()
        )
    }
}

pub struct Barcode {
    barcode: Vec<u8>,
    barcode_rc: Vec<u8>,
}

impl Barcode {
    pub const NBL: &str = "AAGGTTAA";
    pub const NBR: &str = "CAGCACCT";

    pub fn new(barcode: &[u8]) -> Self {
        Barcode {
            barcode: barcode.to_owned(),
            barcode_rc: revcomp(barcode),
        }
    }

    fn patten(&self) -> Vec<u8> {
        self.barcode
            .iter()
            .chain(Self::NBR.as_bytes())
            .map(|x| *x)
            .collect::<Vec<_>>()
    }

    pub fn front_myers(&self) -> Myers {
        Myers::new(self.patten())
    }

    pub fn rear_myers(&self) -> Myers {
        Myers::new(revcomp(self.patten()))
    }
}

impl Display for Barcode {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Front: {} {}\nRear:  {} {}",
            str::from_utf8(&self.barcode).unwrap(),
            Self::NBR,
            str::from_utf8(&*revcomp(Self::NBR.as_bytes())).unwrap(),
            str::from_utf8(&self.barcode_rc).unwrap()
        )
    }
}

pub const BARCODES: [&str; 97] = [
    "TACTTCGTTCAGTTACGTATTGCT", // Ligation Adapter
    "CACAAAGACACCGACAACTTTCTT",
    "ACAGACGACTACAAACGGAATCGA",
    "CCTGGTAACTGGGACACAAGACTC",
    "TAGGGAAACACGATAGAATCCGAA", // barcode04
    "AAGGTTACACAAACCCTGGACAAG",
    "GACTACTTTCTGCCTTTGCGAGAA",
    "AAGGATTCATTCCCACGGTAACAC",
    "ACGTAACTTGGTTTGTTCCCTGAA",
    "AACCAAGACTCGCTGTGCCTAGTT",
    "GAGAGGACAAAGGTTTCAACGCTT",
    "TCCATTCCCTCCGATAGATGAAAC",
    "TCCGATTCTGCTTCTTTCTACCTG",
    "AGAACGACTTCCATACTCGTGTGA",
    "AACGAGTCTCTTGGGACCCATAGA",
    "AGGTCTACCTCGCTAACACCACTG",
    "CGTCAACTGACAGTGGTTCGTACT",
    "ACCCTCCAGGAAAGTACCTCTGAT",
    "CCAAACCCAACAACCTAGATAGGC",
    "GTTCCTCGTGCAGTGTCAAGAGAT",
    "TTGCGTCCTGTTACGAGAACTCAT",
    "GAGCCTCTCATTGTCCGTTCTCTA",
    "ACCACTGCCATGTATCAAAGTACG",
    "CTTACTACCCAGTGAACCTCCTCG",
    "GCATAGTTCTGCATGATGGGTTAG",
    "GTAAGTTGGGTATGCAACGCAATG",
    "CATACAGCGACTACGCATTCTCAT",
    "CGACGGTTAGATTCACCTCTTACA",
    "TGAAACCTAAGAAGGCACCGTATC",
    "CTAGACACCTTGGGTTGACAGACC",
    "TCAGTGAGGATCTACTTCGACCCA",
    "TGCGTACAGCAATCAGTTACATTG",
    "CCAGTAGAAGTCCGACAACGTCAT",
    "CAGACTTGGTACGGTTGGGTAACT",
    "GGACGAAGAACTCAAGTCAAAGGC",
    "CTACTTACGAAGCTGAGGGACTGC",
    "ATGTCCCAGTTAGAGGAGGAAACA",
    "GCTTGCGATTGATGCTTAGTATCA",
    "ACCACAGGAGGACGATACAGAGAA",
    "CCACAGTGTCAACTAGAGCCTCTC",
    "TAGTTTGGATGACCAAGGATAGCC",
    "GGAGTTCGTCCAGAGAAGTACACG",
    "CTACGTGTAAGGCATACCTGCCAG",
    "CTTTCGTTGTTGACTCGACGGTAG",
    "AGTAGAAAGGGTTCCTTCCCACTC",
    "GATCCAACAGAGATGCCTTCAGTG",
    "GCTGTGTTCCACTTCATTCTCCTG",
    "GTGCAACTTTCCCACAGGTAGTTC",
    "CATCTGGAACGTGGTACACCTGTA",
    "ACTGGTGCAGCTTTGAACATCTAG",
    "ATGGACTTTGGTAACTTCCTGCGT",
    "GTTGAATGAGCCTACTGGGTCCTC",
    "TGAGAGACAAGATTGTTCGTGGAC",
    "AGATTCAGACCGTCTCATGCAAAG",
    "CAAGAGCTTTGACTAAGGAGCATG",
    "TGGAAGATGAGACCCTGATCTACG",
    "TCACTACTCAACAGGTGGCATGAA",
    "GCTAGGTCAATCTCCTTCGGAAGT",
    "CAGGTTACTCCTCCGTGAGTCTGA",
    "TCAATCAAGAAGGGAAAGCAAGGT",
    "CATGTTCAACCAAGGCTTCTATGG",
    "AGAGGGTACTATGTGCCTCAGCAC",
    "CACCCACACTTACTTCAGGACGTA",
    "TTCTGAAGTTCCTGGGTCTTGAAC",
    "GACAGACACCGTTCATCGACTTTC",
    "TTCTCAGTCTTCCTCCAGACAAGG",
    "CCGATCCTTGTGGCTTCTAACTTC",
    "GTTTGTCATACTCGTGTGCTCACC",
    "GAATCTAAGCAAACACGAAGGTGG",
    "TACAGTCCGAGCCTCATGTGATCT",
    "ACCGAGATCCTACGAATGGAGTGT",
    "CCTGGGAGCATCAGGTAGTAACAG",
    "TAGCTGACTGTCTTCCATACCGAC",
    "AAGAAACAGGATGACAGAACCCTC",
    "TACAAGCATCCCAACACTTCCACT",
    "GACCATTGTGATGAACCCTGTTGT",
    "ATGCTTGTTACATCAACCCTGGAC",
    "CGACCTGTTTCTCAGGGATACAAC",
    "AACAACCGAACCTTTGAATCAGAA",
    "TCTCGGAGATAGTTCTCACTGCTG",
    "CGGATGAACATAGGATAGCGATTC",
    "CCTCATCTTGTGAAGTTGTTTCGG",
    "ACGGTATGTCGAGTTCCAGGACTA",
    "TGGCTTGATCTAGGTAAGGTCGAA",
    "GTAGTGGACCTAGAACCTGTGCCA",
    "AACGGAGGAGTTAGTTGGATGATC",
    "AGGTGATCCCAACAAGCGTAAGTA",
    "TACATGCTCCTGTTGTTAGGGAGG",
    "TCTTCTACTACCGATCCGAAGCAG",
    "ACAGCATCAATGTTTGGCTAGTTG",
    "GATGTAGAGGGTACGGTTTGAGGC",
    "GGCTCCATAGGAACTCACGCTACT",
    "TTGTGAGTGGAAAGATACAGGACC",
    "AGTTTCCATCACTTCAGACTTGGG",
    "GATTGTCCTCAAACTGCCACCTAC",
    "CCTGTCTGGAAGAAGAATGGACTT",
    "CTGAACGGTCATAGAGTCCACCAT",
];

#[cfg(test)]
mod barcode_primer_test {
    use super::*;

    #[test]
    #[ignore]
    fn bar() {
        let barcode = Barcode::new(BARCODES[1].as_bytes());
        println!("{}", str::from_utf8(&barcode.patten()).unwrap());
        println!("{barcode}")
    }

    #[test]
    fn primer() {
        let primer = Primer::new(
            "test",
            "GCAACAACAACCTTTCATCCT".as_bytes(),
            revcomp("TATTTGACAGGATTTATGTGTA".as_bytes()).as_slice(),
        );
        println!("{primer}")
    }
}
