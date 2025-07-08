use clap::parser::ValueSource;
use std::collections::HashMap;
use std::sync::OnceLock;

static SEQ_INFO: OnceLock<HashMap<&str, &SequenceInfo>> = OnceLock::new();

pub fn get_seq_info() -> &'static HashMap<&'static str, &'static SequenceInfo> {
    SEQ_INFO.get_or_init(|| {
        let mut seq_info = HashMap::from([
            ("LSK", &LSK),
            ("RAD", &RAD),
            ("ULK", &ULK),
            ("NBD_1", &NBD_1),
            ("NBD_2", &NBD_2),
            ("NBD_3", &NBD_3),
            ("NBD_4", &NBD_4),
            ("NBD_5", &NBD_5),
            ("NBD_6", &NBD_6),
            ("NBD_7", &NBD_7),
            ("NBD_8", &NBD_8),
            ("NBD_9", &NBD_9),
            ("NBD_10", &NBD_10),
            ("NBD_11", &NBD_11),
            ("NBD_12", &NBD_12),
            ("NBD_13", &NBD_13),
            ("NBD_14", &NBD_14),
            ("NBD_15", &NBD_15),
            ("NBD_16", &NBD_16),
            ("NBD_17", &NBD_17),
            ("NBD_18", &NBD_18),
            ("NBD_19", &NBD_19),
            ("NBD_20", &NBD_20),
            ("NBD_21", &NBD_21),
            ("NBD_22", &NBD_22),
            ("NBD_23", &NBD_23),
            ("NBD_24", &NBD_24),
            ("NBD_25", &NBD_25),
            ("NBD_26", &NBD_26),
            ("NBD_27", &NBD_27),
            ("NBD_28", &NBD_28),
            ("NBD_29", &NBD_29),
            ("NBD_30", &NBD_30),
            ("NBD_31", &NBD_31),
            ("NBD_32", &NBD_32),
            ("NBD_33", &NBD_33),
            ("NBD_34", &NBD_34),
            ("NBD_35", &NBD_35),
            ("NBD_36", &NBD_36),
            ("NBD_37", &NBD_37),
            ("NBD_38", &NBD_38),
            ("NBD_39", &NBD_39),
            ("NBD_40", &NBD_40),
            ("NBD_41", &NBD_41),
            ("NBD_42", &NBD_42),
            ("NBD_43", &NBD_43),
            ("NBD_44", &NBD_44),
            ("NBD_45", &NBD_45),
            ("NBD_46", &NBD_46),
            ("NBD_47", &NBD_47),
            ("NBD_48", &NBD_48),
            ("NBD_49", &NBD_49),
            ("NBD_50", &NBD_50),
            ("NBD_51", &NBD_51),
            ("NBD_52", &NBD_52),
            ("NBD_53", &NBD_53),
            ("NBD_54", &NBD_54),
            ("NBD_55", &NBD_55),
            ("NBD_56", &NBD_56),
            ("NBD_57", &NBD_57),
            ("NBD_58", &NBD_58),
            ("NBD_59", &NBD_59),
            ("NBD_60", &NBD_60),
            ("NBD_61", &NBD_61),
            ("NBD_62", &NBD_62),
            ("NBD_63", &NBD_63),
            ("NBD_64", &NBD_64),
            ("NBD_65", &NBD_65),
            ("NBD_66", &NBD_66),
            ("NBD_67", &NBD_67),
            ("NBD_68", &NBD_68),
            ("NBD_69", &NBD_69),
            ("NBD_70", &NBD_70),
            ("NBD_71", &NBD_71),
            ("NBD_72", &NBD_72),
            ("NBD_73", &NBD_73),
            ("NBD_74", &NBD_74),
            ("NBD_75", &NBD_75),
            ("NBD_76", &NBD_76),
            ("NBD_77", &NBD_77),
            ("NBD_78", &NBD_78),
            ("NBD_79", &NBD_79),
            ("NBD_80", &NBD_80),
            ("NBD_81", &NBD_81),
            ("NBD_82", &NBD_82),
            ("NBD_83", &NBD_83),
            ("NBD_84", &NBD_84),
            ("NBD_85", &NBD_85),
            ("NBD_86", &NBD_86),
            ("NBD_87", &NBD_87),
            ("NBD_88", &NBD_88),
            ("NBD_89", &NBD_89),
            ("NBD_90", &NBD_90),
            ("NBD_91", &NBD_91),
            ("NBD_92", &NBD_92),
            ("NBD_93", &NBD_93),
            ("NBD_94", &NBD_94),
            ("NBD_95", &NBD_95),
            ("NBD_96", &NBD_96),
            ("RBK_1", &RBK_1),
            ("RBK_2", &RBK_2),
            ("RBK_3", &RBK_3),
            ("RBK_4", &RBK_4),
            ("RBK_5", &RBK_5),
            ("RBK_6", &RBK_6),
            ("RBK_7", &RBK_7),
            ("RBK_8", &RBK_8),
            ("RBK_9", &RBK_9),
            ("RBK_10", &RBK_10),
            ("RBK_11", &RBK_11),
            ("RBK_12", &RBK_12),
            ("RBK_13", &RBK_13),
            ("RBK_14", &RBK_14),
            ("RBK_15", &RBK_15),
            ("RBK_16", &RBK_16),
            ("RBK_17", &RBK_17),
            ("RBK_18", &RBK_18),
            ("RBK_19", &RBK_19),
            ("RBK_20", &RBK_20),
            ("RBK_21", &RBK_21),
            ("RBK_22", &RBK_22),
            ("RBK_23", &RBK_23),
            ("RBK_24", &RBK_24),
            ("RBK_25", &RBK_25),
            ("RBK_26", &RBK_26),
            ("RBK_27", &RBK_27),
            ("RBK_28", &RBK_28),
            ("RBK_29", &RBK_29),
            ("RBK_30", &RBK_30),
            ("RBK_31", &RBK_31),
            ("RBK_32", &RBK_32),
            ("RBK_33", &RBK_33),
            ("RBK_34", &RBK_34),
            ("RBK_35", &RBK_35),
            ("RBK_36", &RBK_36),
            ("RBK_37", &RBK_37),
            ("RBK_38", &RBK_38),
            ("RBK_39", &RBK_39),
            ("RBK_40", &RBK_40),
            ("RBK_41", &RBK_41),
            ("RBK_42", &RBK_42),
            ("RBK_43", &RBK_43),
            ("RBK_44", &RBK_44),
            ("RBK_45", &RBK_45),
            ("RBK_46", &RBK_46),
            ("RBK_47", &RBK_47),
            ("RBK_48", &RBK_48),
            ("RBK_49", &RBK_49),
            ("RBK_50", &RBK_50),
            ("RBK_51", &RBK_51),
            ("RBK_52", &RBK_52),
            ("RBK_53", &RBK_53),
            ("RBK_54", &RBK_54),
            ("RBK_55", &RBK_55),
            ("RBK_56", &RBK_56),
            ("RBK_57", &RBK_57),
            ("RBK_58", &RBK_58),
            ("RBK_59", &RBK_59),
            ("RBK_60", &RBK_60),
            ("RBK_61", &RBK_61),
            ("RBK_62", &RBK_62),
            ("RBK_63", &RBK_63),
            ("RBK_64", &RBK_64),
            ("RBK_65", &RBK_65),
            ("RBK_66", &RBK_66),
            ("RBK_67", &RBK_67),
            ("RBK_68", &RBK_68),
            ("RBK_69", &RBK_69),
            ("RBK_70", &RBK_70),
            ("RBK_71", &RBK_71),
            ("RBK_72", &RBK_72),
            ("RBK_73", &RBK_73),
            ("RBK_74", &RBK_74),
            ("RBK_75", &RBK_75),
            ("RBK_76", &RBK_76),
            ("RBK_77", &RBK_77),
            ("RBK_78", &RBK_78),
            ("RBK_79", &RBK_79),
            ("RBK_80", &RBK_80),
            ("RBK_81", &RBK_81),
            ("RBK_82", &RBK_82),
            ("RBK_83", &RBK_83),
            ("RBK_84", &RBK_84),
            ("RBK_85", &RBK_85),
            ("RBK_86", &RBK_86),
            ("RBK_87", &RBK_87),
            ("RBK_88", &RBK_88),
            ("RBK_89", &RBK_89),
            ("RBK_90", &RBK_90),
            ("RBK_91", &RBK_91),
            ("RBK_92", &RBK_92),
            ("RBK_93", &RBK_93),
            ("RBK_94", &RBK_94),
            ("RBK_95", &RBK_95),
            ("RBK_96", &RBK_96),
            ("PCS", &PCS),
            ("PCB_1", &PCB_1),
            ("PCB_2", &PCB_2),
            ("PCB_3", &PCB_3),
            ("PCB_4", &PCB_4),
            ("PCB_5", &PCB_5),
            ("PCB_6", &PCB_6),
            ("PCB_7", &PCB_7),
            ("PCB_8", &PCB_8),
            ("PCB_9", &PCB_9),
            ("PCB_10", &PCB_10),
            ("PCB_11", &PCB_11),
            ("PCB_12", &PCB_12),
            ("PCB_13", &PCB_13),
            ("PCB_14", &PCB_14),
            ("PCB_15", &PCB_15),
            ("PCB_16", &PCB_16),
            ("PCB_17", &PCB_17),
            ("PCB_18", &PCB_18),
            ("PCB_19", &PCB_19),
            ("PCB_20", &PCB_20),
            ("PCB_21", &PCB_21),
            ("PCB_22", &PCB_22),
            ("PCB_23", &PCB_23),
            ("PCB_24", &PCB_24),
        ]);
        seq_info
    })
}

pub type End = (usize, f64, f64);
#[derive(Clone, Debug)]
pub struct SequenceInfo {
    pub name: &'static str,
    pub end5: (&'static str, End),
    pub end3: Option<(&'static str, End)>,
    pub rev_com_end5: Option<(&'static str, End)>,
    pub rev_com_end3: Option<(&'static str, End)>,
}

impl SequenceInfo {
    fn single_update(
        end: &mut End,
        len: (ValueSource, usize),
        pct: (ValueSource, f64),
        ident: (ValueSource, f64),
    ) {
        if len.0 == ValueSource::CommandLine {
            end.0 = len.1
        }
        if pct.0 == ValueSource::CommandLine {
            end.1 = pct.1;
        }

        if ident.0 == ValueSource::CommandLine {
            end.2 = ident.1;
        }
    }
    pub fn update(
        &mut self,
        end5_len: (ValueSource, usize),
        end5_pct: (ValueSource, f64),
        end5_ident: (ValueSource, f64),
        end3_len: (ValueSource, usize),
        end3_pct: (ValueSource, f64),
        end3_ident: (ValueSource, f64),
        rev_com_end5_len: (ValueSource, usize),
        rev_com_end5_pct: (ValueSource, f64),
        rev_com_end5_ident: (ValueSource, f64),
        rev_com_end3_len: (ValueSource, usize),
        rev_com_end3_pct: (ValueSource, f64),
        rev_com_end3_ident: (ValueSource, f64),
    ) {
        if end5_len.0 == ValueSource::CommandLine {
            self.end5.1.0 = end5_len.1
        }
        if end5_pct.0 == ValueSource::CommandLine {
            self.end5.1.1 = end5_pct.1
        }
        if end5_ident.0 == ValueSource::CommandLine {
            self.end5.1.2 = end5_ident.1
        }

        if self.end3.is_some() {
            let mut this_end3 = self.end3.unwrap().1.clone();
            Self::single_update(&mut this_end3, end3_len, end3_pct, end3_ident);
            self.end3 = Some((self.end3.unwrap().0, this_end3));
        }
        if self.rev_com_end5.is_some() {
            let mut this_rev_com_end5 = self.rev_com_end5.unwrap().1.clone();
            Self::single_update(
                &mut this_rev_com_end5,
                rev_com_end5_len,
                rev_com_end5_pct,
                rev_com_end5_pct,
            );
            self.rev_com_end5 = Some((self.rev_com_end5.unwrap().0, this_rev_com_end5))
        }

        if self.rev_com_end3.is_some() {
            let mut this_rev_com_end3 = self.rev_com_end3.unwrap().1.clone();
            Self::single_update(
                &mut this_rev_com_end3,
                rev_com_end3_len,
                rev_com_end3_pct,
                rev_com_end3_pct,
            );
            self.rev_com_end3 = Some((self.rev_com_end3.unwrap().0, this_rev_com_end3))
        }
    }
}

const LSK_END5: End = (100, 0.5, 0.75);
const LSK_END3: End = (80, 0.5, 0.75);
const RAD_END5: End = (180, 0.5, 0.75);
const NBD_END5: End = (150, 0.6, 0.75);
const NBD_END3: End = (120, 0.6, 0.75);
const RBK_END5: End = (180, 0.5, 0.75);
const PCS_END5: End = (150, 0.6, 0.75);
const PCS_END3: End = (150, 0.4, 0.75);
const PCS_REV_COM_END5: End = (150, 0.6, 0.75);
const PCS_REV_COM_END3: End = (150, 0.4, 0.75);
const PCB_END5: End = (180, 0.4, 0.75);
const PCB_END3: End = (180, 0.3, 0.75);
const PCB_REV_COM_END5: End = (180, 0.5, 0.75);
const PCB_REV_COM_END3: End = (180, 0.3, 0.75);

/*
SQK-LSK114
LSK114 library reads structure
          |--->  LA_ADAPTER_5   <----| | insert Seq | |--->   LA_ADAPTER_3   <---|
5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-..............-AGCAATACGTAACTGAACGAAGTACAGG-3

3' end always is truncated
 */
const LSK: SequenceInfo = SequenceInfo {
    name: "LSK",
    end5: ("CCTGTACTTCGTTCAGTTACGTATTGCT", LSK_END5),
    end3: Some(("AGCAATACGTAACTGAACGAAGTACAGG", LSK_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};

/*
      SQK-RAD114; SQK-ULK114
      the rapid adapter(RA) from ont document is 5’-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-3', but RA_ADAPTER will be used when trimming reads with Rapid Adapter(R    A)
      Always consider only the adapter at 5' end for Rapid library

      structure of reads with RA, but no barcode
        |RA_ADAPTER we want to trim from reads                             | insert Seq
      5-GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA-...................-3

*/
const RAD: SequenceInfo = SequenceInfo {
    name: "RAD",
    end5: (
        "GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RAD_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const ULK: SequenceInfo = SequenceInfo {
    name: "ULK",
    end5: (
        "GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RAD_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};

/*
SQK-NBD114-24; SQK-NBD114-96
NBD114-24/96 library reads structure
Example for Native Barcode01
         |NA_ADAPTER_5                |L_F_5   |Barcode01 rev com       |R_F_5   |insert Seq         |L_F_3   |Barcode01               |R_F_3         |NA_ADAPTER_3
5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT AAGGTTAA CACAAAGACACCGACAACTTTCTT CAGCACCT ................... AGGTGCTG AAGAAAGTTGTCGGTGTCTTTGTG TTAACCTTAGCAAT ACGTAACTGAACGAAGTACAGG-3
we use barcode_left_flanking + barcode + barcode_right_flanking as query to trim nbd reads
*/
const NBD_1: SequenceInfo = SequenceInfo {
    name: "NBD_1",
    end5: ("AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_2: SequenceInfo = SequenceInfo {
    name: "NBD_2",
    end5: ("AAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_3: SequenceInfo = SequenceInfo {
    name: "NBD_3",
    end5: ("AAGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_4: SequenceInfo = SequenceInfo {
    name: "NBD_4",
    end5: ("AAGGTTAATAGGGAAACACGATAGAATCCGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_5: SequenceInfo = SequenceInfo {
    name: "NBD_5",
    end5: ("AAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_6: SequenceInfo = SequenceInfo {
    name: "NBD_6",
    end5: ("AAGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_7: SequenceInfo = SequenceInfo {
    name: "NBD_7",
    end5: ("AAGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_8: SequenceInfo = SequenceInfo {
    name: "NBD_8",
    end5: ("AAGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_9: SequenceInfo = SequenceInfo {
    name: "NBD_9",
    end5: ("AAGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_10: SequenceInfo = SequenceInfo {
    name: "NBD_10",
    end5: ("AAGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_11: SequenceInfo = SequenceInfo {
    name: "NBD_11",
    end5: ("AAGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_12: SequenceInfo = SequenceInfo {
    name: "NBD_12",
    end5: ("AAGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_13: SequenceInfo = SequenceInfo {
    name: "NBD_13",
    end5: ("AAGGTTAAAGAACGACTTCCATACTCGTGTGACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTCACACGAGTATGGAAGTCGTTCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_14: SequenceInfo = SequenceInfo {
    name: "NBD_14",
    end5: ("AAGGTTAAAACGAGTCTCTTGGGACCCATAGACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTCTATGGGTCCCAAGAGACTCGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_15: SequenceInfo = SequenceInfo {
    name: "NBD_15",
    end5: ("AAGGTTAAAGGTCTACCTCGCTAACACCACTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAGTGGTGTTAGCGAGGTAGACCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_16: SequenceInfo = SequenceInfo {
    name: "NBD_16",
    end5: ("AAGGTTAACGTCAACTGACAGTGGTTCGTACTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAGTACGAACCACTGTCAGTTGACGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_17: SequenceInfo = SequenceInfo {
    name: "NBD_17",
    end5: ("AAGGTTAAACCCTCCAGGAAAGTACCTCTGATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATCAGAGGTACTTTCCTGGAGGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_18: SequenceInfo = SequenceInfo {
    name: "NBD_18",
    end5: ("AAGGTTAACCAAACCCAACAACCTAGATAGGCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGCCTATCTAGGTTGTTGGGTTTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_19: SequenceInfo = SequenceInfo {
    name: "NBD_19",
    end5: ("AAGGTTAAGTTCCTCGTGCAGTGTCAAGAGATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATCTCTTGACACTGCACGAGGAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_20: SequenceInfo = SequenceInfo {
    name: "NBD_20",
    end5: ("AAGGTTAATTGCGTCCTGTTACGAGAACTCATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATGAGTTCTCGTAACAGGACGCAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_21: SequenceInfo = SequenceInfo {
    name: "NBD_21",
    end5: ("AAGGTTAAGAGCCTCTCATTGTCCGTTCTCTACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTAGAGAACGGACAATGAGAGGCTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_22: SequenceInfo = SequenceInfo {
    name: "NBD_22",
    end5: ("AAGGTTAAACCACTGCCATGTATCAAAGTACGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCGTACTTTGATACATGGCAGTGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_23: SequenceInfo = SequenceInfo {
    name: "NBD_23",
    end5: ("AAGGTTAACTTACTACCCAGTGAACCTCCTCGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCGAGGAGGTTCACTGGGTAGTAAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_24: SequenceInfo = SequenceInfo {
    name: "NBD_24",
    end5: ("AAGGTTAAGCATAGTTCTGCATGATGGGTTAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTAACCCATCATGCAGAACTATGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_25: SequenceInfo = SequenceInfo {
    name: "NBD_25",
    end5: ("AAGGTTAAGTAAGTTGGGTATGCAACGCAATGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCATTGCGTTGCATACCCAACTTACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_26: SequenceInfo = SequenceInfo {
    name: "NBD_26",
    end5: ("AAGGTTAACATACAGCGACTACGCATTCTCATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATGAGAATGCGTAGTCGCTGTATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_27: SequenceInfo = SequenceInfo {
    name: "NBD_27",
    end5: ("AAGGTTAACGACGGTTAGATTCACCTCTTACACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTGTAAGAGGTGAATCTAACCGTCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_28: SequenceInfo = SequenceInfo {
    name: "NBD_28",
    end5: ("AAGGTTAATGAAACCTAAGAAGGCACCGTATCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGATACGGTGCCTTCTTAGGTTTCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_29: SequenceInfo = SequenceInfo {
    name: "NBD_29",
    end5: ("AAGGTTAACTAGACACCTTGGGTTGACAGACCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGGTCTGTCAACCCAAGGTGTCTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_30: SequenceInfo = SequenceInfo {
    name: "NBD_30",
    end5: ("AAGGTTAATCAGTGAGGATCTACTTCGACCCACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTGGGTCGAAGTAGATCCTCACTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_31: SequenceInfo = SequenceInfo {
    name: "NBD_31",
    end5: ("AAGGTTAATGCGTACAGCAATCAGTTACATTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAATGTAACTGATTGCTGTACGCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_32: SequenceInfo = SequenceInfo {
    name: "NBD_32",
    end5: ("AAGGTTAACCAGTAGAAGTCCGACAACGTCATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATGACGTTGTCGGACTTCTACTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_33: SequenceInfo = SequenceInfo {
    name: "NBD_33",
    end5: ("AAGGTTAACAGACTTGGTACGGTTGGGTAACTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAGTTACCCAACCGTACCAAGTCTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_34: SequenceInfo = SequenceInfo {
    name: "NBD_34",
    end5: ("AAGGTTAAGGACGAAGAACTCAAGTCAAAGGCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGCCTTTGACTTGAGTTCTTCGTCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_35: SequenceInfo = SequenceInfo {
    name: "NBD_35",
    end5: ("AAGGTTAACTACTTACGAAGCTGAGGGACTGCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGCAGTCCCTCAGCTTCGTAAGTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_36: SequenceInfo = SequenceInfo {
    name: "NBD_36",
    end5: ("AAGGTTAAATGTCCCAGTTAGAGGAGGAAACACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTGTTTCCTCCTCTAACTGGGACATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_37: SequenceInfo = SequenceInfo {
    name: "NBD_37",
    end5: ("AAGGTTAAGCTTGCGATTGATGCTTAGTATCACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTGATACTAAGCATCAATCGCAAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_38: SequenceInfo = SequenceInfo {
    name: "NBD_38",
    end5: ("AAGGTTAAACCACAGGAGGACGATACAGAGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCTCTGTATCGTCCTCCTGTGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_39: SequenceInfo = SequenceInfo {
    name: "NBD_39",
    end5: ("AAGGTTAACCACAGTGTCAACTAGAGCCTCTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAGAGGCTCTAGTTGACACTGTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_40: SequenceInfo = SequenceInfo {
    name: "NBD_40",
    end5: ("AAGGTTAATAGTTTGGATGACCAAGGATAGCCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGGCTATCCTTGGTCATCCAAACTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_41: SequenceInfo = SequenceInfo {
    name: "NBD_41",
    end5: ("AAGGTTAAGGAGTTCGTCCAGAGAAGTACACGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCGTGTACTTCTCTGGACGAACTCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_42: SequenceInfo = SequenceInfo {
    name: "NBD_42",
    end5: ("AAGGTTAACTACGTGTAAGGCATACCTGCCAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTGGCAGGTATGCCTTACACGTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_43: SequenceInfo = SequenceInfo {
    name: "NBD_43",
    end5: ("AAGGTTAACTTTCGTTGTTGACTCGACGGTAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTACCGTCGAGTCAACAACGAAAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_44: SequenceInfo = SequenceInfo {
    name: "NBD_44",
    end5: ("AAGGTTAAAGTAGAAAGGGTTCCTTCCCACTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAGTGGGAAGGAACCCTTTCTACTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_45: SequenceInfo = SequenceInfo {
    name: "NBD_45",
    end5: ("AAGGTTAAGATCCAACAGAGATGCCTTCAGTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCACTGAAGGCATCTCTGTTGGATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_46: SequenceInfo = SequenceInfo {
    name: "NBD_46",
    end5: ("AAGGTTAAGCTGTGTTCCACTTCATTCTCCTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAGGAGAATGAAGTGGAACACAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_47: SequenceInfo = SequenceInfo {
    name: "NBD_47",
    end5: ("AAGGTTAAGTGCAACTTTCCCACAGGTAGTTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAACTACCTGTGGGAAAGTTGCACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_48: SequenceInfo = SequenceInfo {
    name: "NBD_48",
    end5: ("AAGGTTAACATCTGGAACGTGGTACACCTGTACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTACAGGTGTACCACGTTCCAGATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_49: SequenceInfo = SequenceInfo {
    name: "NBD_49",
    end5: ("AAGGTTAAACTGGTGCAGCTTTGAACATCTAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTAGATGTTCAAAGCTGCACCAGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_50: SequenceInfo = SequenceInfo {
    name: "NBD_50",
    end5: ("AAGGTTAAATGGACTTTGGTAACTTCCTGCGTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGACGCAGGAAGTTACCAAAGTCCATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_51: SequenceInfo = SequenceInfo {
    name: "NBD_51",
    end5: ("AAGGTTAAGTTGAATGAGCCTACTGGGTCCTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAGGACCCAGTAGGCTCATTCAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_52: SequenceInfo = SequenceInfo {
    name: "NBD_52",
    end5: ("AAGGTTAATGAGAGACAAGATTGTTCGTGGACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTCCACGAACAATCTTGTCTCTCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_53: SequenceInfo = SequenceInfo {
    name: "NBD_53",
    end5: ("AAGGTTAAAGATTCAGACCGTCTCATGCAAAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTTTGCATGAGACGGTCTGAATCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_54: SequenceInfo = SequenceInfo {
    name: "NBD_54",
    end5: ("AAGGTTAACAAGAGCTTTGACTAAGGAGCATGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCATGCTCCTTAGTCAAAGCTCTTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_55: SequenceInfo = SequenceInfo {
    name: "NBD_55",
    end5: ("AAGGTTAATGGAAGATGAGACCCTGATCTACGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCGTAGATCAGGGTCTCATCTTCCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_56: SequenceInfo = SequenceInfo {
    name: "NBD_56",
    end5: ("AAGGTTAATCACTACTCAACAGGTGGCATGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCATGCCACCTGTTGAGTAGTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_57: SequenceInfo = SequenceInfo {
    name: "NBD_57",
    end5: ("AAGGTTAAGCTAGGTCAATCTCCTTCGGAAGTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGACTTCCGAAGGAGATTGACCTAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_58: SequenceInfo = SequenceInfo {
    name: "NBD_58",
    end5: ("AAGGTTAACAGGTTACTCCTCCGTGAGTCTGACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTCAGACTCACGGAGGAGTAACCTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_59: SequenceInfo = SequenceInfo {
    name: "NBD_59",
    end5: ("AAGGTTAATCAATCAAGAAGGGAAAGCAAGGTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGACCTTGCTTTCCCTTCTTGATTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_60: SequenceInfo = SequenceInfo {
    name: "NBD_60",
    end5: ("AAGGTTAACATGTTCAACCAAGGCTTCTATGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCATAGAAGCCTTGGTTGAACATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_61: SequenceInfo = SequenceInfo {
    name: "NBD_61",
    end5: ("AAGGTTAAAGAGGGTACTATGTGCCTCAGCACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTGCTGAGGCACATAGTACCCTCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_62: SequenceInfo = SequenceInfo {
    name: "NBD_62",
    end5: ("AAGGTTAACACCCACACTTACTTCAGGACGTACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTACGTCCTGAAGTAAGTGTGGGTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_63: SequenceInfo = SequenceInfo {
    name: "NBD_63",
    end5: ("AAGGTTAATTCTGAAGTTCCTGGGTCTTGAACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTTCAAGACCCAGGAACTTCAGAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_64: SequenceInfo = SequenceInfo {
    name: "NBD_64",
    end5: ("AAGGTTAAGACAGACACCGTTCATCGACTTTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAAAGTCGATGAACGGTGTCTGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_65: SequenceInfo = SequenceInfo {
    name: "NBD_65",
    end5: ("AAGGTTAATTCTCAGTCTTCCTCCAGACAAGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCTTGTCTGGAGGAAGACTGAGAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_66: SequenceInfo = SequenceInfo {
    name: "NBD_66",
    end5: ("AAGGTTAACCGATCCTTGTGGCTTCTAACTTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAAGTTAGAAGCCACAAGGATCGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_67: SequenceInfo = SequenceInfo {
    name: "NBD_67",
    end5: ("AAGGTTAAGTTTGTCATACTCGTGTGCTCACCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGGTGAGCACACGAGTATGACAAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_68: SequenceInfo = SequenceInfo {
    name: "NBD_68",
    end5: ("AAGGTTAAGAATCTAAGCAAACACGAAGGTGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCACCTTCGTGTTTGCTTAGATTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_69: SequenceInfo = SequenceInfo {
    name: "NBD_69",
    end5: ("AAGGTTAATACAGTCCGAGCCTCATGTGATCTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAGATCACATGAGGCTCGGACTGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_70: SequenceInfo = SequenceInfo {
    name: "NBD_70",
    end5: ("AAGGTTAAACCGAGATCCTACGAATGGAGTGTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGACACTCCATTCGTAGGATCTCGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_71: SequenceInfo = SequenceInfo {
    name: "NBD_71",
    end5: ("AAGGTTAACCTGGGAGCATCAGGTAGTAACAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTGTTACTACCTGATGCTCCCAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_72: SequenceInfo = SequenceInfo {
    name: "NBD_72",
    end5: ("AAGGTTAATAGCTGACTGTCTTCCATACCGACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTCGGTATGGAAGACAGTCAGCTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_73: SequenceInfo = SequenceInfo {
    name: "NBD_73",
    end5: ("AAGGTTAAAAGAAACAGGATGACAGAACCCTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAGGGTTCTGTCATCCTGTTTCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_74: SequenceInfo = SequenceInfo {
    name: "NBD_74",
    end5: ("AAGGTTAATACAAGCATCCCAACACTTCCACTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAGTGGAAGTGTTGGGATGCTTGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_75: SequenceInfo = SequenceInfo {
    name: "NBD_75",
    end5: ("AAGGTTAAGACCATTGTGATGAACCCTGTTGTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGACAACAGGGTTCATCACAATGGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_76: SequenceInfo = SequenceInfo {
    name: "NBD_76",
    end5: ("AAGGTTAAATGCTTGTTACATCAACCCTGGACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTCCAGGGTTGATGTAACAAGCATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_77: SequenceInfo = SequenceInfo {
    name: "NBD_77",
    end5: ("AAGGTTAACGACCTGTTTCTCAGGGATACAACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTTGTATCCCTGAGAAACAGGTCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_78: SequenceInfo = SequenceInfo {
    name: "NBD_78",
    end5: ("AAGGTTAAAACAACCGAACCTTTGAATCAGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCTGATTCAAAGGTTCGGTTGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_79: SequenceInfo = SequenceInfo {
    name: "NBD_79",
    end5: ("AAGGTTAATCTCGGAGATAGTTCTCACTGCTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAGCAGTGAGAACTATCTCCGAGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_80: SequenceInfo = SequenceInfo {
    name: "NBD_80",
    end5: ("AAGGTTAACGGATGAACATAGGATAGCGATTCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGAATCGCTATCCTATGTTCATCCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_81: SequenceInfo = SequenceInfo {
    name: "NBD_81",
    end5: ("AAGGTTAACCTCATCTTGTGAAGTTGTTTCGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCGAAACAACTTCACAAGATGAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_82: SequenceInfo = SequenceInfo {
    name: "NBD_82",
    end5: ("AAGGTTAAACGGTATGTCGAGTTCCAGGACTACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTAGTCCTGGAACTCGACATACCGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_83: SequenceInfo = SequenceInfo {
    name: "NBD_83",
    end5: ("AAGGTTAATGGCTTGATCTAGGTAAGGTCGAACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTTCGACCTTACCTAGATCAAGCCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_84: SequenceInfo = SequenceInfo {
    name: "NBD_84",
    end5: ("AAGGTTAAGTAGTGGACCTAGAACCTGTGCCACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTGGCACAGGTTCTAGGTCCACTACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_85: SequenceInfo = SequenceInfo {
    name: "NBD_85",
    end5: ("AAGGTTAAAACGGAGGAGTTAGTTGGATGATCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGATCATCCAACTAACTCCTCCGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_86: SequenceInfo = SequenceInfo {
    name: "NBD_86",
    end5: ("AAGGTTAAAGGTGATCCCAACAAGCGTAAGTACAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGTACTTACGCTTGTTGGGATCACCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_87: SequenceInfo = SequenceInfo {
    name: "NBD_87",
    end5: ("AAGGTTAATACATGCTCCTGTTGTTAGGGAGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCTCCCTAACAACAGGAGCATGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_88: SequenceInfo = SequenceInfo {
    name: "NBD_88",
    end5: ("AAGGTTAATCTTCTACTACCGATCCGAAGCAGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCTGCTTCGGATCGGTAGTAGAAGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_89: SequenceInfo = SequenceInfo {
    name: "NBD_89",
    end5: ("AAGGTTAAACAGCATCAATGTTTGGCTAGTTGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCAACTAGCCAAACATTGATGCTGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_90: SequenceInfo = SequenceInfo {
    name: "NBD_90",
    end5: ("AAGGTTAAGATGTAGAGGGTACGGTTTGAGGCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGCCTCAAACCGTACCCTCTACATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_91: SequenceInfo = SequenceInfo {
    name: "NBD_91",
    end5: ("AAGGTTAAGGCTCCATAGGAACTCACGCTACTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAGTAGCGTGAGTTCCTATGGAGCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_92: SequenceInfo = SequenceInfo {
    name: "NBD_92",
    end5: ("AAGGTTAATTGTGAGTGGAAAGATACAGGACCCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGGTCCTGTATCTTTCCACTCACAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_93: SequenceInfo = SequenceInfo {
    name: "NBD_93",
    end5: ("AAGGTTAAAGTTTCCATCACTTCAGACTTGGGCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGCCCAAGTCTGAAGTGATGGAAACTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_94: SequenceInfo = SequenceInfo {
    name: "NBD_94",
    end5: ("AAGGTTAAGATTGTCCTCAAACTGCCACCTACCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGGTAGGTGGCAGTTTGAGGACAATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_95: SequenceInfo = SequenceInfo {
    name: "NBD_95",
    end5: ("AAGGTTAACCTGTCTGGAAGAAGAATGGACTTCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGAAGTCCATTCTTCTTCCAGACAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_96: SequenceInfo = SequenceInfo {
    name: "NBD_96",
    end5: ("AAGGTTAACTGAACGGTCATAGAGTCCACCATCAGCACCT", NBD_END5),
    end3: Some(("AGGTGCTGATGGTGGACTCTATGACCGTTCAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_1: SequenceInfo = SequenceInfo {
    name: "RBK_1",
    end5: (
        "AAGAAAGTTGTCGGTGTCTTTGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_2: SequenceInfo = SequenceInfo {
    name: "RBK_2",
    end5: (
        "TCGATTCCGTTTGTAGTCGTCTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_3: SequenceInfo = SequenceInfo {
    name: "RBK_3",
    end5: (
        "GAGTCTTGTGTCCCAGTTACCAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_4: SequenceInfo = SequenceInfo {
    name: "RBK_4",
    end5: (
        "TTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_5: SequenceInfo = SequenceInfo {
    name: "RBK_5",
    end5: (
        "CTTGTCCAGGGTTTGTGTAACCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_6: SequenceInfo = SequenceInfo {
    name: "RBK_6",
    end5: (
        "TTCTCGCAAAGGCAGAAAGTAGTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_7: SequenceInfo = SequenceInfo {
    name: "RBK_7",
    end5: (
        "GTGTTACCGTGGGAATGAATCCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_8: SequenceInfo = SequenceInfo {
    name: "RBK_8",
    end5: (
        "TTCAGGGAACAAACCAAGTTACGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_9: SequenceInfo = SequenceInfo {
    name: "RBK_9",
    end5: (
        "AACTAGGCACAGCGAGTCTTGGTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_10: SequenceInfo = SequenceInfo {
    name: "RBK_10",
    end5: (
        "AAGCGTTGAAACCTTTGTCCTCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_11: SequenceInfo = SequenceInfo {
    name: "RBK_11",
    end5: (
        "GTTTCATCTATCGGAGGGAATGGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_12: SequenceInfo = SequenceInfo {
    name: "RBK_12",
    end5: (
        "CAGGTAGAAAGAAGCAGAATCGGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_13: SequenceInfo = SequenceInfo {
    name: "RBK_13",
    end5: (
        "AGAACGACTTCCATACTCGTGTGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_14: SequenceInfo = SequenceInfo {
    name: "RBK_14",
    end5: (
        "AACGAGTCTCTTGGGACCCATAGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_15: SequenceInfo = SequenceInfo {
    name: "RBK_15",
    end5: (
        "AGGTCTACCTCGCTAACACCACTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_16: SequenceInfo = SequenceInfo {
    name: "RBK_16",
    end5: (
        "CGTCAACTGACAGTGGTTCGTACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_17: SequenceInfo = SequenceInfo {
    name: "RBK_17",
    end5: (
        "ACCCTCCAGGAAAGTACCTCTGATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_18: SequenceInfo = SequenceInfo {
    name: "RBK_18",
    end5: (
        "CCAAACCCAACAACCTAGATAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_19: SequenceInfo = SequenceInfo {
    name: "RBK_19",
    end5: (
        "GTTCCTCGTGCAGTGTCAAGAGATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_20: SequenceInfo = SequenceInfo {
    name: "RBK_20",
    end5: (
        "TTGCGTCCTGTTACGAGAACTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_21: SequenceInfo = SequenceInfo {
    name: "RBK_21",
    end5: (
        "GAGCCTCTCATTGTCCGTTCTCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_22: SequenceInfo = SequenceInfo {
    name: "RBK_22",
    end5: (
        "ACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_23: SequenceInfo = SequenceInfo {
    name: "RBK_23",
    end5: (
        "CTTACTACCCAGTGAACCTCCTCGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_24: SequenceInfo = SequenceInfo {
    name: "RBK_24",
    end5: (
        "GCATAGTTCTGCATGATGGGTTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_25: SequenceInfo = SequenceInfo {
    name: "RBK_25",
    end5: (
        "GTAAGTTGGGTATGCAACGCAATGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_26: SequenceInfo = SequenceInfo {
    name: "RBK_26",
    end5: (
        "CATACAGCGACTACGCATTCTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_27: SequenceInfo = SequenceInfo {
    name: "RBK_27",
    end5: (
        "CGACGGTTAGATTCACCTCTTACAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_28: SequenceInfo = SequenceInfo {
    name: "RBK_28",
    end5: (
        "TGAAACCTAAGAAGGCACCGTATCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_29: SequenceInfo = SequenceInfo {
    name: "RBK_29",
    end5: (
        "CTAGACACCTTGGGTTGACAGACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_30: SequenceInfo = SequenceInfo {
    name: "RBK_30",
    end5: (
        "TCAGTGAGGATCTACTTCGACCCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_31: SequenceInfo = SequenceInfo {
    name: "RBK_31",
    end5: (
        "TGCGTACAGCAATCAGTTACATTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_32: SequenceInfo = SequenceInfo {
    name: "RBK_32",
    end5: (
        "CCAGTAGAAGTCCGACAACGTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_33: SequenceInfo = SequenceInfo {
    name: "RBK_33",
    end5: (
        "CAGACTTGGTACGGTTGGGTAACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_34: SequenceInfo = SequenceInfo {
    name: "RBK_34",
    end5: (
        "GGACGAAGAACTCAAGTCAAAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_35: SequenceInfo = SequenceInfo {
    name: "RBK_35",
    end5: (
        "CTACTTACGAAGCTGAGGGACTGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_36: SequenceInfo = SequenceInfo {
    name: "RBK_36",
    end5: (
        "ATGTCCCAGTTAGAGGAGGAAACAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_37: SequenceInfo = SequenceInfo {
    name: "RBK_37",
    end5: (
        "GCTTGCGATTGATGCTTAGTATCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_38: SequenceInfo = SequenceInfo {
    name: "RBK_38",
    end5: (
        "ACCACAGGAGGACGATACAGAGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_39: SequenceInfo = SequenceInfo {
    name: "RBK_39",
    end5: (
        "CCACAGTGTCAACTAGAGCCTCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_40: SequenceInfo = SequenceInfo {
    name: "RBK_40",
    end5: (
        "TAGTTTGGATGACCAAGGATAGCCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_41: SequenceInfo = SequenceInfo {
    name: "RBK_41",
    end5: (
        "GGAGTTCGTCCAGAGAAGTACACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_42: SequenceInfo = SequenceInfo {
    name: "RBK_42",
    end5: (
        "CTACGTGTAAGGCATACCTGCCAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_43: SequenceInfo = SequenceInfo {
    name: "RBK_43",
    end5: (
        "CTTTCGTTGTTGACTCGACGGTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_44: SequenceInfo = SequenceInfo {
    name: "RBK_44",
    end5: (
        "AGTAGAAAGGGTTCCTTCCCACTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_45: SequenceInfo = SequenceInfo {
    name: "RBK_45",
    end5: (
        "GATCCAACAGAGATGCCTTCAGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_46: SequenceInfo = SequenceInfo {
    name: "RBK_46",
    end5: (
        "GCTGTGTTCCACTTCATTCTCCTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_47: SequenceInfo = SequenceInfo {
    name: "RBK_47",
    end5: (
        "GTGCAACTTTCCCACAGGTAGTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_48: SequenceInfo = SequenceInfo {
    name: "RBK_48",
    end5: (
        "CATCTGGAACGTGGTACACCTGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_49: SequenceInfo = SequenceInfo {
    name: "RBK_49",
    end5: (
        "ACTGGTGCAGCTTTGAACATCTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_50: SequenceInfo = SequenceInfo {
    name: "RBK_50",
    end5: (
        "ATGGACTTTGGTAACTTCCTGCGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_51: SequenceInfo = SequenceInfo {
    name: "RBK_51",
    end5: (
        "GTTGAATGAGCCTACTGGGTCCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_52: SequenceInfo = SequenceInfo {
    name: "RBK_52",
    end5: (
        "TGAGAGACAAGATTGTTCGTGGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_53: SequenceInfo = SequenceInfo {
    name: "RBK_53",
    end5: (
        "AGATTCAGACCGTCTCATGCAAAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_54: SequenceInfo = SequenceInfo {
    name: "RBK_54",
    end5: (
        "CAAGAGCTTTGACTAAGGAGCATGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_55: SequenceInfo = SequenceInfo {
    name: "RBK_55",
    end5: (
        "TGGAAGATGAGACCCTGATCTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_56: SequenceInfo = SequenceInfo {
    name: "RBK_56",
    end5: (
        "TCACTACTCAACAGGTGGCATGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_57: SequenceInfo = SequenceInfo {
    name: "RBK_57",
    end5: (
        "GCTAGGTCAATCTCCTTCGGAAGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_58: SequenceInfo = SequenceInfo {
    name: "RBK_58",
    end5: (
        "CAGGTTACTCCTCCGTGAGTCTGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_59: SequenceInfo = SequenceInfo {
    name: "RBK_59",
    end5: (
        "TCAATCAAGAAGGGAAAGCAAGGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_60: SequenceInfo = SequenceInfo {
    name: "RBK_60",
    end5: (
        "CATGTTCAACCAAGGCTTCTATGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_61: SequenceInfo = SequenceInfo {
    name: "RBK_61",
    end5: (
        "AGAGGGTACTATGTGCCTCAGCACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_62: SequenceInfo = SequenceInfo {
    name: "RBK_62",
    end5: (
        "CACCCACACTTACTTCAGGACGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_63: SequenceInfo = SequenceInfo {
    name: "RBK_63",
    end5: (
        "TTCTGAAGTTCCTGGGTCTTGAACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_64: SequenceInfo = SequenceInfo {
    name: "RBK_64",
    end5: (
        "GACAGACACCGTTCATCGACTTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_65: SequenceInfo = SequenceInfo {
    name: "RBK_65",
    end5: (
        "TTCTCAGTCTTCCTCCAGACAAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_66: SequenceInfo = SequenceInfo {
    name: "RBK_66",
    end5: (
        "CCGATCCTTGTGGCTTCTAACTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_67: SequenceInfo = SequenceInfo {
    name: "RBK_67",
    end5: (
        "GTTTGTCATACTCGTGTGCTCACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_68: SequenceInfo = SequenceInfo {
    name: "RBK_68",
    end5: (
        "GAATCTAAGCAAACACGAAGGTGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_69: SequenceInfo = SequenceInfo {
    name: "RBK_69",
    end5: (
        "TACAGTCCGAGCCTCATGTGATCTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_70: SequenceInfo = SequenceInfo {
    name: "RBK_70",
    end5: (
        "ACCGAGATCCTACGAATGGAGTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_71: SequenceInfo = SequenceInfo {
    name: "RBK_71",
    end5: (
        "CCTGGGAGCATCAGGTAGTAACAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_72: SequenceInfo = SequenceInfo {
    name: "RBK_72",
    end5: (
        "TAGCTGACTGTCTTCCATACCGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_73: SequenceInfo = SequenceInfo {
    name: "RBK_73",
    end5: (
        "AAGAAACAGGATGACAGAACCCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_74: SequenceInfo = SequenceInfo {
    name: "RBK_74",
    end5: (
        "TACAAGCATCCCAACACTTCCACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_75: SequenceInfo = SequenceInfo {
    name: "RBK_75",
    end5: (
        "GACCATTGTGATGAACCCTGTTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_76: SequenceInfo = SequenceInfo {
    name: "RBK_76",
    end5: (
        "ATGCTTGTTACATCAACCCTGGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_77: SequenceInfo = SequenceInfo {
    name: "RBK_77",
    end5: (
        "CGACCTGTTTCTCAGGGATACAACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_78: SequenceInfo = SequenceInfo {
    name: "RBK_78",
    end5: (
        "AACAACCGAACCTTTGAATCAGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_79: SequenceInfo = SequenceInfo {
    name: "RBK_79",
    end5: (
        "TCTCGGAGATAGTTCTCACTGCTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_80: SequenceInfo = SequenceInfo {
    name: "RBK_80",
    end5: (
        "CGGATGAACATAGGATAGCGATTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_81: SequenceInfo = SequenceInfo {
    name: "RBK_81",
    end5: (
        "CCTCATCTTGTGAAGTTGTTTCGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_82: SequenceInfo = SequenceInfo {
    name: "RBK_82",
    end5: (
        "ACGGTATGTCGAGTTCCAGGACTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_83: SequenceInfo = SequenceInfo {
    name: "RBK_83",
    end5: (
        "TGGCTTGATCTAGGTAAGGTCGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_84: SequenceInfo = SequenceInfo {
    name: "RBK_84",
    end5: (
        "GTAGTGGACCTAGAACCTGTGCCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_85: SequenceInfo = SequenceInfo {
    name: "RBK_85",
    end5: (
        "AACGGAGGAGTTAGTTGGATGATCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_86: SequenceInfo = SequenceInfo {
    name: "RBK_86",
    end5: (
        "AGGTGATCCCAACAAGCGTAAGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_87: SequenceInfo = SequenceInfo {
    name: "RBK_87",
    end5: (
        "TACATGCTCCTGTTGTTAGGGAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_88: SequenceInfo = SequenceInfo {
    name: "RBK_88",
    end5: (
        "TCTTCTACTACCGATCCGAAGCAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_89: SequenceInfo = SequenceInfo {
    name: "RBK_89",
    end5: (
        "ACAGCATCAATGTTTGGCTAGTTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_90: SequenceInfo = SequenceInfo {
    name: "RBK_90",
    end5: (
        "GATGTAGAGGGTACGGTTTGAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_91: SequenceInfo = SequenceInfo {
    name: "RBK_91",
    end5: (
        "GGCTCCATAGGAACTCACGCTACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_92: SequenceInfo = SequenceInfo {
    name: "RBK_92",
    end5: (
        "TTGTGAGTGGAAAGATACAGGACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_93: SequenceInfo = SequenceInfo {
    name: "RBK_93",
    end5: (
        "AGTTTCCATCACTTCAGACTTGGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_94: SequenceInfo = SequenceInfo {
    name: "RBK_94",
    end5: (
        "GATTGTCCTCAAACTGCCACCTACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_95: SequenceInfo = SequenceInfo {
    name: "RBK_95",
    end5: (
        "CCTGTCTGGAAGAAGAATGGACTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const RBK_96: SequenceInfo = SequenceInfo {
    name: "RBK_96",
    end5: (
        "CTGAACGGTCATAGAGTCCACCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    ),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};

/*
cDNA-PCR Sequencing Kit: SQK-PCS114
cDNA-PCR Barcoding Kit V14: SQK-PCB114.24
Strand Switching Primer II (SSPII)
cDNA RT Adapter (CRTA)

SQK-PCS114 structure
     |SSPII                                                | insert Seq with polyA  | CRTA
5-...TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG .........AAAAAAAAAAAAAAA CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG...-3
3-...CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA .........TTTTTTTTTTTTTTT GAACGCCCGCCGCCTGAGAGGAGACTTCTATCTCGCTGTCCGTTC...-5


SQK-PCB114.24 structure
     | BP01                   | SSPII                                               | insert Seq with polyA | CRTA                                             | BP01 reverse com
5-...AAGAAAGTTGTCGGTGTCTTTGTG TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG .........AAAAAAAAAAAAAAA CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT     CACAAAGACACCGACAACTTTCTT...-3
3-...TTCTTTCAACAGCCACAGAAACAC AAAGACAACCACGACTATAACGAAABBBBAABBBBAABBBBAABBBBAAACCC .........TTTTTTTTTTTTTTT GAACGCCCGCCGCCTGAGAGGAGACTTCTATCTCGCTGTCCGTTCA     GTGTTTCTGTGGCTGTTGAAAGAA...-5

The cDNA RT Adapter (CRTA) is a double stranded adapter with a poly(T) overhang
which anneals to the very end of the poly(A) tail of the RNA strand.
This ensures that the full length of the RNA is reverse transcribed and
that the poly(A) length can be estimated accurately.
Annealing Buffer (AB) has been included to improve CRTA ligation.
The full structure of CRTA is like below:
CRTA:                      5'-CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT-3'
CRTA_REV_COM:   3'-TTTTTTTTTTTGAACGCCCGCCGCCTGAGAGGAGACTTCTATCTCGCTGTCCGTTCA-5'
*/
const PCS: SequenceInfo = SequenceInfo {
    name: "PCS",
    end5: (
        "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCS_END5,
    ),
    end3: Some(("CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT", PCS_END3)),
    rev_com_end5: Some((
        "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCS_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA",
        PCS_REV_COM_END3,
    )),
};
const PCB_1: SequenceInfo = SequenceInfo {
    name: "PCB_1",
    end5: (
        "AAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCACAAAGACACCGACAACTTTCTT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACACAAAGACACCGACAACTTTCTT",
        PCB_REV_COM_END3,
    )),
};
const PCB_2: SequenceInfo = SequenceInfo {
    name: "PCB_2",
    end5: (
        "TCGATTCCGTTTGTAGTCGTCTGTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTACAGACGACTACAAACGGAATCGA",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "TCGATTCCGTTTGTAGTCGTCTGTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAACAGACGACTACAAACGGAATCGA",
        PCB_REV_COM_END3,
    )),
};
const PCB_3: SequenceInfo = SequenceInfo {
    name: "PCB_3",
    end5: (
        "GAGTCTTGTGTCCCAGTTACCAGGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCCTGGTAACTGGGACACAAGACTC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GAGTCTTGTGTCCCAGTTACCAGGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACCTGGTAACTGGGACACAAGACTC",
        PCB_REV_COM_END3,
    )),
};
const PCB_4: SequenceInfo = SequenceInfo {
    name: "PCB_4",
    end5: (
        "TTCGGATTCTATCGTGTTTCCCTATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTAGGGAAACACGATAGAATCCGAA",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "TTCGGATTCTATCGTGTTTCCCTAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATAGGGAAACACGATAGAATCCGAA",
        PCB_REV_COM_END3,
    )),
};
const PCB_5: SequenceInfo = SequenceInfo {
    name: "PCB_5",
    end5: (
        "CTTGTCCAGGGTTTGTGTAACCTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAAGGTTACACAAACCCTGGACAAG",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "CTTGTCCAGGGTTTGTGTAACCTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAAGGTTACACAAACCCTGGACAAG",
        PCB_REV_COM_END3,
    )),
};
const PCB_6: SequenceInfo = SequenceInfo {
    name: "PCB_6",
    end5: (
        "TTCTCGCAAAGGCAGAAAGTAGTCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGACTACTTTCTGCCTTTGCGAGAA",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "TTCTCGCAAAGGCAGAAAGTAGTCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGACTACTTTCTGCCTTTGCGAGAA",
        PCB_REV_COM_END3,
    )),
};
const PCB_7: SequenceInfo = SequenceInfo {
    name: "PCB_7",
    end5: (
        "GTGTTACCGTGGGAATGAATCCTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAAGGATTCATTCCCACGGTAACAC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GTGTTACCGTGGGAATGAATCCTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAAGGATTCATTCCCACGGTAACAC",
        PCB_REV_COM_END3,
    )),
};
const PCB_8: SequenceInfo = SequenceInfo {
    name: "PCB_8",
    end5: (
        "TTCAGGGAACAAACCAAGTTACGTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTACGTAACTTGGTTTGTTCCCTGAA",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "TTCAGGGAACAAACCAAGTTACGTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAACGTAACTTGGTTTGTTCCCTGAA",
        PCB_REV_COM_END3,
    )),
};
const PCB_9: SequenceInfo = SequenceInfo {
    name: "PCB_9",
    end5: (
        "AACTAGGCACAGCGAGTCTTGGTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAACCAAGACTCGCTGTGCCTAGTT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AACTAGGCACAGCGAGTCTTGGTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAACCAAGACTCGCTGTGCCTAGTT",
        PCB_REV_COM_END3,
    )),
};
const PCB_10: SequenceInfo = SequenceInfo {
    name: "PCB_10",
    end5: (
        "AAGCGTTGAAACCTTTGTCCTCTCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGAGAGGACAAAGGTTTCAACGCTT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AAGCGTTGAAACCTTTGTCCTCTCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGAGAGGACAAAGGTTTCAACGCTT",
        PCB_REV_COM_END3,
    )),
};
const PCB_11: SequenceInfo = SequenceInfo {
    name: "PCB_11",
    end5: (
        "GTTTCATCTATCGGAGGGAATGGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCCATTCCCTCCGATAGATGAAAC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GTTTCATCTATCGGAGGGAATGGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCCATTCCCTCCGATAGATGAAAC",
        PCB_REV_COM_END3,
    )),
};
const PCB_12: SequenceInfo = SequenceInfo {
    name: "PCB_12",
    end5: (
        "CAGGTAGAAAGAAGCAGAATCGGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCCGATTCTGCTTCTTTCTACCTG",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "CAGGTAGAAAGAAGCAGAATCGGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCCGATTCTGCTTCTTTCTACCTG",
        PCB_REV_COM_END3,
    )),
};
const PCB_13: SequenceInfo = SequenceInfo {
    name: "PCB_13",
    end5: (
        "AGAACGACTTCCATACTCGTGTGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCACACGAGTATGGAAGTCGTTCT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AGAACGACTTCCATACTCGTGTGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCACACGAGTATGGAAGTCGTTCT",
        PCB_REV_COM_END3,
    )),
};
const PCB_14: SequenceInfo = SequenceInfo {
    name: "PCB_14",
    end5: (
        "AACGAGTCTCTTGGGACCCATAGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCTATGGGTCCCAAGAGACTCGTT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AACGAGTCTCTTGGGACCCATAGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCTATGGGTCCCAAGAGACTCGTT",
        PCB_REV_COM_END3,
    )),
};
const PCB_15: SequenceInfo = SequenceInfo {
    name: "PCB_15",
    end5: (
        "AGGTCTACCTCGCTAACACCACTGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCAGTGGTGTTAGCGAGGTAGACCT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "AGGTCTACCTCGCTAACACCACTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACAGTGGTGTTAGCGAGGTAGACCT",
        PCB_REV_COM_END3,
    )),
};
const PCB_16: SequenceInfo = SequenceInfo {
    name: "PCB_16",
    end5: (
        "CGTCAACTGACAGTGGTTCGTACTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAGTACGAACCACTGTCAGTTGACG",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "CGTCAACTGACAGTGGTTCGTACTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAGTACGAACCACTGTCAGTTGACG",
        PCB_REV_COM_END3,
    )),
};
const PCB_17: SequenceInfo = SequenceInfo {
    name: "PCB_17",
    end5: (
        "ACCCTCCAGGAAAGTACCTCTGATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATCAGAGGTACTTTCCTGGAGGGT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "ACCCTCCAGGAAAGTACCTCTGATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATCAGAGGTACTTTCCTGGAGGGT",
        PCB_REV_COM_END3,
    )),
};
const PCB_18: SequenceInfo = SequenceInfo {
    name: "PCB_18",
    end5: (
        "CCAAACCCAACAACCTAGATAGGCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGCCTATCTAGGTTGTTGGGTTTGG",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "CCAAACCCAACAACCTAGATAGGCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGCCTATCTAGGTTGTTGGGTTTGG",
        PCB_REV_COM_END3,
    )),
};
const PCB_19: SequenceInfo = SequenceInfo {
    name: "PCB_19",
    end5: (
        "GTTCCTCGTGCAGTGTCAAGAGATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATCTCTTGACACTGCACGAGGAAC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GTTCCTCGTGCAGTGTCAAGAGATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATCTCTTGACACTGCACGAGGAAC",
        PCB_REV_COM_END3,
    )),
};
const PCB_20: SequenceInfo = SequenceInfo {
    name: "PCB_20",
    end5: (
        "TTGCGTCCTGTTACGAGAACTCATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATGAGTTCTCGTAACAGGACGCAA",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "TTGCGTCCTGTTACGAGAACTCATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATGAGTTCTCGTAACAGGACGCAA",
        PCB_REV_COM_END3,
    )),
};
const PCB_21: SequenceInfo = SequenceInfo {
    name: "PCB_21",
    end5: (
        "GAGCCTCTCATTGTCCGTTCTCTATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTAGAGAACGGACAATGAGAGGCTC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GAGCCTCTCATTGTCCGTTCTCTAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATAGAGAACGGACAATGAGAGGCTC",
        PCB_REV_COM_END3,
    )),
};
const PCB_22: SequenceInfo = SequenceInfo {
    name: "PCB_22",
    end5: (
        "ACCACTGCCATGTATCAAAGTACGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCGTACTTTGATACATGGCAGTGGT",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "ACCACTGCCATGTATCAAAGTACGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACGTACTTTGATACATGGCAGTGGT",
        PCB_REV_COM_END3,
    )),
};
const PCB_23: SequenceInfo = SequenceInfo {
    name: "PCB_23",
    end5: (
        "CTTACTACCCAGTGAACCTCCTCGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCGAGGAGGTTCACTGGGTAGTAAG",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "CTTACTACCCAGTGAACCTCCTCGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACGAGGAGGTTCACTGGGTAGTAAG",
        PCB_REV_COM_END3,
    )),
};
const PCB_24: SequenceInfo = SequenceInfo {
    name: "PCB_24",
    end5: (
        "GCATAGTTCTGCATGATGGGTTAGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    ),
    end3: Some((
        "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCTAACCCATCATGCAGAACTATGC",
        PCB_END3,
    )),
    rev_com_end5: Some((
        "GCATAGTTCTGCATGATGGGTTAGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACTAACCCATCATGCAGAACTATGC",
        PCB_REV_COM_END3,
    )),
};
