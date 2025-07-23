use clap::parser::ValueSource;
use std::collections::HashMap;
use std::iter::repeat_n;

pub fn get_trim_cfg<'a>() -> HashMap<&'a str, TrimConfig<'a>> {
    HashMap::from([
        ("LSK", LSK),
        ("RAD", RAD),
        ("ULK", ULK),
        ("NBD_1", NBD_1),
        ("NBD_2", NBD_2),
        ("NBD_3", NBD_3),
        ("NBD_4", NBD_4),
        ("NBD_5", NBD_5),
        ("NBD_6", NBD_6),
        ("NBD_7", NBD_7),
        ("NBD_8", NBD_8),
        ("NBD_9", NBD_9),
        ("NBD_10", NBD_10),
        ("NBD_11", NBD_11),
        ("NBD_12", NBD_12),
        ("NBD_13", NBD_13),
        ("NBD_14", NBD_14),
        ("NBD_15", NBD_15),
        ("NBD_16", NBD_16),
        ("NBD_17", NBD_17),
        ("NBD_18", NBD_18),
        ("NBD_19", NBD_19),
        ("NBD_20", NBD_20),
        ("NBD_21", NBD_21),
        ("NBD_22", NBD_22),
        ("NBD_23", NBD_23),
        ("NBD_24", NBD_24),
        ("NBD_25", NBD_25),
        ("NBD_26", NBD_26),
        ("NBD_27", NBD_27),
        ("NBD_28", NBD_28),
        ("NBD_29", NBD_29),
        ("NBD_30", NBD_30),
        ("NBD_31", NBD_31),
        ("NBD_32", NBD_32),
        ("NBD_33", NBD_33),
        ("NBD_34", NBD_34),
        ("NBD_35", NBD_35),
        ("NBD_36", NBD_36),
        ("NBD_37", NBD_37),
        ("NBD_38", NBD_38),
        ("NBD_39", NBD_39),
        ("NBD_40", NBD_40),
        ("NBD_41", NBD_41),
        ("NBD_42", NBD_42),
        ("NBD_43", NBD_43),
        ("NBD_44", NBD_44),
        ("NBD_45", NBD_45),
        ("NBD_46", NBD_46),
        ("NBD_47", NBD_47),
        ("NBD_48", NBD_48),
        ("NBD_49", NBD_49),
        ("NBD_50", NBD_50),
        ("NBD_51", NBD_51),
        ("NBD_52", NBD_52),
        ("NBD_53", NBD_53),
        ("NBD_54", NBD_54),
        ("NBD_55", NBD_55),
        ("NBD_56", NBD_56),
        ("NBD_57", NBD_57),
        ("NBD_58", NBD_58),
        ("NBD_59", NBD_59),
        ("NBD_60", NBD_60),
        ("NBD_61", NBD_61),
        ("NBD_62", NBD_62),
        ("NBD_63", NBD_63),
        ("NBD_64", NBD_64),
        ("NBD_65", NBD_65),
        ("NBD_66", NBD_66),
        ("NBD_67", NBD_67),
        ("NBD_68", NBD_68),
        ("NBD_69", NBD_69),
        ("NBD_70", NBD_70),
        ("NBD_71", NBD_71),
        ("NBD_72", NBD_72),
        ("NBD_73", NBD_73),
        ("NBD_74", NBD_74),
        ("NBD_75", NBD_75),
        ("NBD_76", NBD_76),
        ("NBD_77", NBD_77),
        ("NBD_78", NBD_78),
        ("NBD_79", NBD_79),
        ("NBD_80", NBD_80),
        ("NBD_81", NBD_81),
        ("NBD_82", NBD_82),
        ("NBD_83", NBD_83),
        ("NBD_84", NBD_84),
        ("NBD_85", NBD_85),
        ("NBD_86", NBD_86),
        ("NBD_87", NBD_87),
        ("NBD_88", NBD_88),
        ("NBD_89", NBD_89),
        ("NBD_90", NBD_90),
        ("NBD_91", NBD_91),
        ("NBD_92", NBD_92),
        ("NBD_93", NBD_93),
        ("NBD_94", NBD_94),
        ("NBD_95", NBD_95),
        ("NBD_96", NBD_96),
        ("RBK", RBK_1),
        // ("RBK_2", &RBK_2),
        // ("RBK_3", &RBK_3),
        // ("RBK_4", &RBK_4),
        // ("RBK_5", &RBK_5),
        // ("RBK_6", &RBK_6),
        // ("RBK_7", &RBK_7),
        // ("RBK_8", &RBK_8),
        // ("RBK_9", &RBK_9),
        // ("RBK_10", &RBK_10),
        // ("RBK_11", &RBK_11),
        // ("RBK_12", &RBK_12),
        // ("RBK_13", &RBK_13),
        // ("RBK_14", &RBK_14),
        // ("RBK_15", &RBK_15),
        // ("RBK_16", &RBK_16),
        // ("RBK_17", &RBK_17),
        // ("RBK_18", &RBK_18),
        // ("RBK_19", &RBK_19),
        // ("RBK_20", &RBK_20),
        // ("RBK_21", &RBK_21),
        // ("RBK_22", &RBK_22),
        // ("RBK_23", &RBK_23),
        // ("RBK_24", &RBK_24),
        // ("RBK_25", &RBK_25),
        // ("RBK_26", &RBK_26),
        // ("RBK_27", &RBK_27),
        // ("RBK_28", &RBK_28),
        // ("RBK_29", &RBK_29),
        // ("RBK_30", &RBK_30),
        // ("RBK_31", &RBK_31),
        // ("RBK_32", &RBK_32),
        // ("RBK_33", &RBK_33),
        // ("RBK_34", &RBK_34),
        // ("RBK_35", &RBK_35),
        // ("RBK_36", &RBK_36),
        // ("RBK_37", &RBK_37),
        // ("RBK_38", &RBK_38),
        // ("RBK_39", &RBK_39),
        // ("RBK_40", &RBK_40),
        // ("RBK_41", &RBK_41),
        // ("RBK_42", &RBK_42),
        // ("RBK_43", &RBK_43),
        // ("RBK_44", &RBK_44),
        // ("RBK_45", &RBK_45),
        // ("RBK_46", &RBK_46),
        // ("RBK_47", &RBK_47),
        // ("RBK_48", &RBK_48),
        // ("RBK_49", &RBK_49),
        // ("RBK_50", &RBK_50),
        // ("RBK_51", &RBK_51),
        // ("RBK_52", &RBK_52),
        // ("RBK_53", &RBK_53),
        // ("RBK_54", &RBK_54),
        // ("RBK_55", &RBK_55),
        // ("RBK_56", &RBK_56),
        // ("RBK_57", &RBK_57),
        // ("RBK_58", &RBK_58),
        // ("RBK_59", &RBK_59),
        // ("RBK_60", &RBK_60),
        // ("RBK_61", &RBK_61),
        // ("RBK_62", &RBK_62),
        // ("RBK_63", &RBK_63),
        // ("RBK_64", &RBK_64),
        // ("RBK_65", &RBK_65),
        // ("RBK_66", &RBK_66),
        // ("RBK_67", &RBK_67),
        // ("RBK_68", &RBK_68),
        // ("RBK_69", &RBK_69),
        // ("RBK_70", &RBK_70),
        // ("RBK_71", &RBK_71),
        // ("RBK_72", &RBK_72),
        // ("RBK_73", &RBK_73),
        // ("RBK_74", &RBK_74),
        // ("RBK_75", &RBK_75),
        // ("RBK_76", &RBK_76),
        // ("RBK_77", &RBK_77),
        // ("RBK_78", &RBK_78),
        // ("RBK_79", &RBK_79),
        // ("RBK_80", &RBK_80),
        // ("RBK_81", &RBK_81),
        // ("RBK_82", &RBK_82),
        // ("RBK_83", &RBK_83),
        // ("RBK_84", &RBK_84),
        // ("RBK_85", &RBK_85),
        // ("RBK_86", &RBK_86),
        // ("RBK_87", &RBK_87),
        // ("RBK_88", &RBK_88),
        // ("RBK_89", &RBK_89),
        // ("RBK_90", &RBK_90),
        // ("RBK_91", &RBK_91),
        // ("RBK_92", &RBK_92),
        // ("RBK_93", &RBK_93),
        // ("RBK_94", &RBK_94),
        // ("RBK_95", &RBK_95),
        // ("RBK_96", &RBK_96),
        ("PCS", PCS),
        ("PCB", PCB),
        // ("PCB_2", &PCB_2),
        // ("PCB_3", &PCB_3),
        // ("PCB_4", &PCB_4),
        // ("PCB_5", &PCB_5),
        // ("PCB_6", &PCB_6),
        // ("PCB_7", &PCB_7),
        // ("PCB_8", &PCB_8),
        // ("PCB_9", &PCB_9),
        // ("PCB_10", &PCB_10),
        // ("PCB_11", &PCB_11),
        // ("PCB_12", &PCB_12),
        // ("PCB_13", &PCB_13),
        // ("PCB_14", &PCB_14),
        // ("PCB_15", &PCB_15),
        // ("PCB_16", &PCB_16),
        // ("PCB_17", &PCB_17),
        // ("PCB_18", &PCB_18),
        // ("PCB_19", &PCB_19),
        // ("PCB_20", &PCB_20),
        // ("PCB_21", &PCB_21),
        // ("PCB_22", &PCB_22),
        // ("PCB_23", &PCB_23),
        // ("PCB_24", &PCB_24),
    ])
}

pub type EndAlignPara = (usize, f64, f64);
pub type EndConfig<'a> = Option<(&'a str, EndAlignPara)>;
#[derive(Clone, Debug)]
pub struct TrimConfig<'a> {
    pub kit_name: &'static str,
    pub end5: EndConfig<'a>,
    pub end3: EndConfig<'a>,
    pub rev_com_end5: EndConfig<'a>,
    pub rev_com_end3: EndConfig<'a>,
}

impl<'a> Default for TrimConfig<'a> {
    fn default() -> Self {
        TrimConfig {
            kit_name: "customer",
            end5: None,
            end3: None,
            rev_com_end5: None,
            rev_com_end3: None,
        }
    }
}

impl<'a> TrimConfig<'a> {
    fn single_update(
        end: &mut EndAlignPara,
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

    pub fn get_info(&self) -> String {
        let mut info = String::new();
        info.push_str(self.kit_name);
        info.push('\n');
        if self.end5.is_some() {
            let end5 = self.end5.as_ref().unwrap();
            info.push_str(&format!(
                "Expect sequence1 in 5'end: {}, length: {}, config: {:?}\n",
                end5.0,
                end5.0.len(),
                end5.1
            ));
        }
        if self.end3.is_some() {
            let end3 = self.end3.as_ref().unwrap();
            info.push_str(&format!(
                "Expect sequence1 in 3'end: {}, length: {}, config: {:?}\n",
                end3.0,
                end3.0.len(),
                end3.1
            ));
        }

        if self.rev_com_end5.is_some() {
            let rev_com_end5 = self.rev_com_end5.as_ref().unwrap();
            info.push_str(&format!(
                "Expect sequence2 in 5'end: {}, length: {}, config: {:?}\n",
                rev_com_end5.0,
                rev_com_end5.0.len(),
                rev_com_end5.1
            ));
        }
        if self.rev_com_end3.is_some() {
            let rev_com_end3 = self.rev_com_end3.as_ref().unwrap();
            info.push_str(&format!(
                "Expect sequence2 in 3'end: {}, length: {}, config: {:?}\n",
                rev_com_end3.0,
                rev_com_end3.0.len(),
                rev_com_end3.1
            ));
        }
        let x: String = repeat_n('=', 100).collect::<String>() + "\n";
        info.push_str(&x);
        info
    }

    pub fn get_dim(&self) -> (usize, usize) {
        let mut ref_lengths = vec![];
        let mut reads_lengths = vec![];
        if let Some(end5) = self.end5.as_ref() {
            ref_lengths.push(end5.0.len());
            reads_lengths.push(end5.1.0);
        }

        if let Some(end3) = self.end3.as_ref() {
            ref_lengths.push(end3.0.len());
            reads_lengths.push(end3.1.0);
        }

        if let Some(rev_com_end5) = self.rev_com_end5.as_ref() {
            ref_lengths.push(rev_com_end5.0.len());
            reads_lengths.push(rev_com_end5.1.0);
        }

        if let Some(rev_com_end3) = self.rev_com_end3.as_ref() {
            ref_lengths.push(rev_com_end3.0.len());
            reads_lengths.push(rev_com_end3.1.0);
        }
        (
            *ref_lengths.iter().max().unwrap(),
            *reads_lengths.iter().max().unwrap(),
        )
    }

    #[inline]
    pub fn may_trim_end5(&self) -> bool {
        self.end5.is_some()
    }

    #[inline]
    pub fn may_trim_end3(&self) -> bool {
        self.end3.is_some()
    }

    #[inline]
    pub fn may_trim_rev_com_end5(&self) -> bool {
        self.rev_com_end5.is_some()
    }

    #[inline]
    pub fn may_trim_rev_com_end3(&self) -> bool {
        self.rev_com_end3.is_some()
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
        if self.end5.is_some() {
            let mut this_end5 = self.end5.as_ref().unwrap().1.clone();
            Self::single_update(&mut this_end5, end5_len, end5_pct, end5_ident);
            self.end5.as_mut().map(|x| x.1 = this_end5);
        }
        if self.end3.is_some() {
            let mut this_end3 = self.end3.as_ref().unwrap().1.clone();
            Self::single_update(&mut this_end3, end3_len, end3_pct, end3_ident);
            self.end3.as_mut().map(|x| x.1 = this_end3);
            // self.end3 = Some((self.end3.unwrap().0, this_end3));
        }
        if self.rev_com_end5.is_some() {
            let mut this_rev_com_end5 = self.rev_com_end5.as_ref().unwrap().1.clone();
            Self::single_update(
                &mut this_rev_com_end5,
                rev_com_end5_len,
                rev_com_end5_pct,
                rev_com_end5_ident,
            );
            self.rev_com_end5.as_mut().map(|x| x.1 = this_rev_com_end5);
            // self.rev_com_end5 = Some((self.rev_com_end5.unwrap().0, this_rev_com_end5))
        }

        if self.rev_com_end3.is_some() {
            let mut this_rev_com_end3 = self.rev_com_end3.as_ref().unwrap().1.clone();
            Self::single_update(
                &mut this_rev_com_end3,
                rev_com_end3_len,
                rev_com_end3_pct,
                rev_com_end3_ident,
            );
            self.rev_com_end3.as_mut().map(|x| x.1 = this_rev_com_end3);
            // self.rev_com_end3 = Some((self.rev_com_end3.unwrap().0, this_rev_com_end3))
        }
    }
}

const LSK_END5: EndAlignPara = (100, 0.75, 0.8);
const LSK_END3: EndAlignPara = (80, 0.6, 0.8);
const RAD_END5: EndAlignPara = (180, 0.75, 0.8);
const NBD_END5: EndAlignPara = (150, 0.75, 0.8);
const NBD_END3: EndAlignPara = (120, 0.6, 0.8);
const RBK_END5: EndAlignPara = (180, 0.75, 0.8);
const PCS_END5: EndAlignPara = (150, 0.6, 0.75);
const PCS_END3: EndAlignPara = (120, 0.6, 0.75);
const PCS_REV_COM_END5: EndAlignPara = (150, 0.6, 0.75);
const PCS_REV_COM_END3: EndAlignPara = (120, 0.4, 0.75);
const PCB_END5: EndAlignPara = (180, 0.6, 0.75);
const PCB_END3: EndAlignPara = (120, 0.6, 0.75);
const PCB_REV_COM_END5: EndAlignPara = (180, 0.6, 0.75);
const PCB_REV_COM_END3: EndAlignPara = (120, 0.6, 0.75);

/*
SQK-LSK114
LSK114 library reads structure
          |--->  LA_ADAPTER_5   <----| | insert Seq | |--->   LA_ADAPTER_3   <---|
5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-..............-AGCAATACGTAACTGAACGAAGTACAGG-3

3' end always is truncated
 */
const LSK: TrimConfig = TrimConfig {
    kit_name: "LSK",
    end5: Some(("CCTGTACTTCGTTCAGTTACGTATTGCT", LSK_END5)),
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
const RAD: TrimConfig = TrimConfig {
    kit_name: "RAD",
    end5: Some((
        "GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RAD_END5,
    )),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
const ULK: TrimConfig = TrimConfig {
    kit_name: "ULK",
    end5: Some((
        "GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RAD_END5,
    )),
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
const NBD_1: TrimConfig = TrimConfig {
    kit_name: "NBD_1",
    end5: Some(("AAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_2: TrimConfig = TrimConfig {
    kit_name: "NBD_2",
    end5: Some(("AAGGTTAAACAGACGACTACAAACGGAATCGACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTCGATTCCGTTTGTAGTCGTCTGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_3: TrimConfig = TrimConfig {
    kit_name: "NBD_3",
    end5: Some(("AAGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAGTCTTGTGTCCCAGTTACCAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_4: TrimConfig = TrimConfig {
    kit_name: "NBD_4",
    end5: Some(("AAGGTTAATAGGGAAACACGATAGAATCCGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCGGATTCTATCGTGTTTCCCTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_5: TrimConfig = TrimConfig {
    kit_name: "NBD_5",
    end5: Some(("AAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTTGTCCAGGGTTTGTGTAACCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_6: TrimConfig = TrimConfig {
    kit_name: "NBD_6",
    end5: Some(("AAGGTTAAGACTACTTTCTGCCTTTGCGAGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCTCGCAAAGGCAGAAAGTAGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_7: TrimConfig = TrimConfig {
    kit_name: "NBD_7",
    end5: Some(("AAGGTTAAAAGGATTCATTCCCACGGTAACACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_8: TrimConfig = TrimConfig {
    kit_name: "NBD_8",
    end5: Some(("AAGGTTAAACGTAACTTGGTTTGTTCCCTGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCAGGGAACAAACCAAGTTACGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_9: TrimConfig = TrimConfig {
    kit_name: "NBD_9",
    end5: Some(("AAGGTTAAAACCAAGACTCGCTGTGCCTAGTTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAACTAGGCACAGCGAGTCTTGGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_10: TrimConfig = TrimConfig {
    kit_name: "NBD_10",
    end5: Some(("AAGGTTAAGAGAGGACAAAGGTTTCAACGCTTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAAGCGTTGAAACCTTTGTCCTCTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_11: TrimConfig = TrimConfig {
    kit_name: "NBD_11",
    end5: Some(("AAGGTTAATCCATTCCCTCCGATAGATGAAACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTTTCATCTATCGGAGGGAATGGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_12: TrimConfig = TrimConfig {
    kit_name: "NBD_12",
    end5: Some(("AAGGTTAATCCGATTCTGCTTCTTTCTACCTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAGGTAGAAAGAAGCAGAATCGGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_13: TrimConfig = TrimConfig {
    kit_name: "NBD_13",
    end5: Some(("AAGGTTAAAGAACGACTTCCATACTCGTGTGACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTCACACGAGTATGGAAGTCGTTCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_14: TrimConfig = TrimConfig {
    kit_name: "NBD_14",
    end5: Some(("AAGGTTAAAACGAGTCTCTTGGGACCCATAGACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTCTATGGGTCCCAAGAGACTCGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_15: TrimConfig = TrimConfig {
    kit_name: "NBD_15",
    end5: Some(("AAGGTTAAAGGTCTACCTCGCTAACACCACTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAGTGGTGTTAGCGAGGTAGACCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_16: TrimConfig = TrimConfig {
    kit_name: "NBD_16",
    end5: Some(("AAGGTTAACGTCAACTGACAGTGGTTCGTACTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAGTACGAACCACTGTCAGTTGACGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_17: TrimConfig = TrimConfig {
    kit_name: "NBD_17",
    end5: Some(("AAGGTTAAACCCTCCAGGAAAGTACCTCTGATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATCAGAGGTACTTTCCTGGAGGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_18: TrimConfig = TrimConfig {
    kit_name: "NBD_18",
    end5: Some(("AAGGTTAACCAAACCCAACAACCTAGATAGGCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGCCTATCTAGGTTGTTGGGTTTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_19: TrimConfig = TrimConfig {
    kit_name: "NBD_19",
    end5: Some(("AAGGTTAAGTTCCTCGTGCAGTGTCAAGAGATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATCTCTTGACACTGCACGAGGAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_20: TrimConfig = TrimConfig {
    kit_name: "NBD_20",
    end5: Some(("AAGGTTAATTGCGTCCTGTTACGAGAACTCATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATGAGTTCTCGTAACAGGACGCAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_21: TrimConfig = TrimConfig {
    kit_name: "NBD_21",
    end5: Some(("AAGGTTAAGAGCCTCTCATTGTCCGTTCTCTACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTAGAGAACGGACAATGAGAGGCTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_22: TrimConfig = TrimConfig {
    kit_name: "NBD_22",
    end5: Some(("AAGGTTAAACCACTGCCATGTATCAAAGTACGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCGTACTTTGATACATGGCAGTGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_23: TrimConfig = TrimConfig {
    kit_name: "NBD_23",
    end5: Some(("AAGGTTAACTTACTACCCAGTGAACCTCCTCGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCGAGGAGGTTCACTGGGTAGTAAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_24: TrimConfig = TrimConfig {
    kit_name: "NBD_24",
    end5: Some(("AAGGTTAAGCATAGTTCTGCATGATGGGTTAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTAACCCATCATGCAGAACTATGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_25: TrimConfig = TrimConfig {
    kit_name: "NBD_25",
    end5: Some(("AAGGTTAAGTAAGTTGGGTATGCAACGCAATGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCATTGCGTTGCATACCCAACTTACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_26: TrimConfig = TrimConfig {
    kit_name: "NBD_26",
    end5: Some(("AAGGTTAACATACAGCGACTACGCATTCTCATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATGAGAATGCGTAGTCGCTGTATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_27: TrimConfig = TrimConfig {
    kit_name: "NBD_27",
    end5: Some(("AAGGTTAACGACGGTTAGATTCACCTCTTACACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTGTAAGAGGTGAATCTAACCGTCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_28: TrimConfig = TrimConfig {
    kit_name: "NBD_28",
    end5: Some(("AAGGTTAATGAAACCTAAGAAGGCACCGTATCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGATACGGTGCCTTCTTAGGTTTCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_29: TrimConfig = TrimConfig {
    kit_name: "NBD_29",
    end5: Some(("AAGGTTAACTAGACACCTTGGGTTGACAGACCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGGTCTGTCAACCCAAGGTGTCTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_30: TrimConfig = TrimConfig {
    kit_name: "NBD_30",
    end5: Some(("AAGGTTAATCAGTGAGGATCTACTTCGACCCACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTGGGTCGAAGTAGATCCTCACTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_31: TrimConfig = TrimConfig {
    kit_name: "NBD_31",
    end5: Some(("AAGGTTAATGCGTACAGCAATCAGTTACATTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAATGTAACTGATTGCTGTACGCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_32: TrimConfig = TrimConfig {
    kit_name: "NBD_32",
    end5: Some(("AAGGTTAACCAGTAGAAGTCCGACAACGTCATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATGACGTTGTCGGACTTCTACTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_33: TrimConfig = TrimConfig {
    kit_name: "NBD_33",
    end5: Some(("AAGGTTAACAGACTTGGTACGGTTGGGTAACTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAGTTACCCAACCGTACCAAGTCTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_34: TrimConfig = TrimConfig {
    kit_name: "NBD_34",
    end5: Some(("AAGGTTAAGGACGAAGAACTCAAGTCAAAGGCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGCCTTTGACTTGAGTTCTTCGTCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_35: TrimConfig = TrimConfig {
    kit_name: "NBD_35",
    end5: Some(("AAGGTTAACTACTTACGAAGCTGAGGGACTGCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGCAGTCCCTCAGCTTCGTAAGTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_36: TrimConfig = TrimConfig {
    kit_name: "NBD_36",
    end5: Some(("AAGGTTAAATGTCCCAGTTAGAGGAGGAAACACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTGTTTCCTCCTCTAACTGGGACATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_37: TrimConfig = TrimConfig {
    kit_name: "NBD_37",
    end5: Some(("AAGGTTAAGCTTGCGATTGATGCTTAGTATCACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTGATACTAAGCATCAATCGCAAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_38: TrimConfig = TrimConfig {
    kit_name: "NBD_38",
    end5: Some(("AAGGTTAAACCACAGGAGGACGATACAGAGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCTCTGTATCGTCCTCCTGTGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_39: TrimConfig = TrimConfig {
    kit_name: "NBD_39",
    end5: Some(("AAGGTTAACCACAGTGTCAACTAGAGCCTCTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAGAGGCTCTAGTTGACACTGTGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_40: TrimConfig = TrimConfig {
    kit_name: "NBD_40",
    end5: Some(("AAGGTTAATAGTTTGGATGACCAAGGATAGCCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGGCTATCCTTGGTCATCCAAACTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_41: TrimConfig = TrimConfig {
    kit_name: "NBD_41",
    end5: Some(("AAGGTTAAGGAGTTCGTCCAGAGAAGTACACGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCGTGTACTTCTCTGGACGAACTCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_42: TrimConfig = TrimConfig {
    kit_name: "NBD_42",
    end5: Some(("AAGGTTAACTACGTGTAAGGCATACCTGCCAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTGGCAGGTATGCCTTACACGTAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_43: TrimConfig = TrimConfig {
    kit_name: "NBD_43",
    end5: Some(("AAGGTTAACTTTCGTTGTTGACTCGACGGTAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTACCGTCGAGTCAACAACGAAAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_44: TrimConfig = TrimConfig {
    kit_name: "NBD_44",
    end5: Some(("AAGGTTAAAGTAGAAAGGGTTCCTTCCCACTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAGTGGGAAGGAACCCTTTCTACTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_45: TrimConfig = TrimConfig {
    kit_name: "NBD_45",
    end5: Some(("AAGGTTAAGATCCAACAGAGATGCCTTCAGTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCACTGAAGGCATCTCTGTTGGATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_46: TrimConfig = TrimConfig {
    kit_name: "NBD_46",
    end5: Some(("AAGGTTAAGCTGTGTTCCACTTCATTCTCCTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAGGAGAATGAAGTGGAACACAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_47: TrimConfig = TrimConfig {
    kit_name: "NBD_47",
    end5: Some(("AAGGTTAAGTGCAACTTTCCCACAGGTAGTTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAACTACCTGTGGGAAAGTTGCACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_48: TrimConfig = TrimConfig {
    kit_name: "NBD_48",
    end5: Some(("AAGGTTAACATCTGGAACGTGGTACACCTGTACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTACAGGTGTACCACGTTCCAGATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_49: TrimConfig = TrimConfig {
    kit_name: "NBD_49",
    end5: Some(("AAGGTTAAACTGGTGCAGCTTTGAACATCTAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTAGATGTTCAAAGCTGCACCAGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_50: TrimConfig = TrimConfig {
    kit_name: "NBD_50",
    end5: Some(("AAGGTTAAATGGACTTTGGTAACTTCCTGCGTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGACGCAGGAAGTTACCAAAGTCCATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_51: TrimConfig = TrimConfig {
    kit_name: "NBD_51",
    end5: Some(("AAGGTTAAGTTGAATGAGCCTACTGGGTCCTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAGGACCCAGTAGGCTCATTCAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_52: TrimConfig = TrimConfig {
    kit_name: "NBD_52",
    end5: Some(("AAGGTTAATGAGAGACAAGATTGTTCGTGGACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTCCACGAACAATCTTGTCTCTCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_53: TrimConfig = TrimConfig {
    kit_name: "NBD_53",
    end5: Some(("AAGGTTAAAGATTCAGACCGTCTCATGCAAAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTTTGCATGAGACGGTCTGAATCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_54: TrimConfig = TrimConfig {
    kit_name: "NBD_54",
    end5: Some(("AAGGTTAACAAGAGCTTTGACTAAGGAGCATGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCATGCTCCTTAGTCAAAGCTCTTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_55: TrimConfig = TrimConfig {
    kit_name: "NBD_55",
    end5: Some(("AAGGTTAATGGAAGATGAGACCCTGATCTACGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCGTAGATCAGGGTCTCATCTTCCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_56: TrimConfig = TrimConfig {
    kit_name: "NBD_56",
    end5: Some(("AAGGTTAATCACTACTCAACAGGTGGCATGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCATGCCACCTGTTGAGTAGTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_57: TrimConfig = TrimConfig {
    kit_name: "NBD_57",
    end5: Some(("AAGGTTAAGCTAGGTCAATCTCCTTCGGAAGTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGACTTCCGAAGGAGATTGACCTAGCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_58: TrimConfig = TrimConfig {
    kit_name: "NBD_58",
    end5: Some(("AAGGTTAACAGGTTACTCCTCCGTGAGTCTGACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTCAGACTCACGGAGGAGTAACCTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_59: TrimConfig = TrimConfig {
    kit_name: "NBD_59",
    end5: Some(("AAGGTTAATCAATCAAGAAGGGAAAGCAAGGTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGACCTTGCTTTCCCTTCTTGATTGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_60: TrimConfig = TrimConfig {
    kit_name: "NBD_60",
    end5: Some(("AAGGTTAACATGTTCAACCAAGGCTTCTATGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCATAGAAGCCTTGGTTGAACATGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_61: TrimConfig = TrimConfig {
    kit_name: "NBD_61",
    end5: Some(("AAGGTTAAAGAGGGTACTATGTGCCTCAGCACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTGCTGAGGCACATAGTACCCTCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_62: TrimConfig = TrimConfig {
    kit_name: "NBD_62",
    end5: Some(("AAGGTTAACACCCACACTTACTTCAGGACGTACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTACGTCCTGAAGTAAGTGTGGGTGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_63: TrimConfig = TrimConfig {
    kit_name: "NBD_63",
    end5: Some(("AAGGTTAATTCTGAAGTTCCTGGGTCTTGAACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTTCAAGACCCAGGAACTTCAGAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_64: TrimConfig = TrimConfig {
    kit_name: "NBD_64",
    end5: Some(("AAGGTTAAGACAGACACCGTTCATCGACTTTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAAAGTCGATGAACGGTGTCTGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_65: TrimConfig = TrimConfig {
    kit_name: "NBD_65",
    end5: Some(("AAGGTTAATTCTCAGTCTTCCTCCAGACAAGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCTTGTCTGGAGGAAGACTGAGAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_66: TrimConfig = TrimConfig {
    kit_name: "NBD_66",
    end5: Some(("AAGGTTAACCGATCCTTGTGGCTTCTAACTTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAAGTTAGAAGCCACAAGGATCGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_67: TrimConfig = TrimConfig {
    kit_name: "NBD_67",
    end5: Some(("AAGGTTAAGTTTGTCATACTCGTGTGCTCACCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGGTGAGCACACGAGTATGACAAACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_68: TrimConfig = TrimConfig {
    kit_name: "NBD_68",
    end5: Some(("AAGGTTAAGAATCTAAGCAAACACGAAGGTGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCACCTTCGTGTTTGCTTAGATTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_69: TrimConfig = TrimConfig {
    kit_name: "NBD_69",
    end5: Some(("AAGGTTAATACAGTCCGAGCCTCATGTGATCTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAGATCACATGAGGCTCGGACTGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_70: TrimConfig = TrimConfig {
    kit_name: "NBD_70",
    end5: Some(("AAGGTTAAACCGAGATCCTACGAATGGAGTGTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGACACTCCATTCGTAGGATCTCGGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_71: TrimConfig = TrimConfig {
    kit_name: "NBD_71",
    end5: Some(("AAGGTTAACCTGGGAGCATCAGGTAGTAACAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTGTTACTACCTGATGCTCCCAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_72: TrimConfig = TrimConfig {
    kit_name: "NBD_72",
    end5: Some(("AAGGTTAATAGCTGACTGTCTTCCATACCGACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTCGGTATGGAAGACAGTCAGCTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_73: TrimConfig = TrimConfig {
    kit_name: "NBD_73",
    end5: Some(("AAGGTTAAAAGAAACAGGATGACAGAACCCTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAGGGTTCTGTCATCCTGTTTCTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_74: TrimConfig = TrimConfig {
    kit_name: "NBD_74",
    end5: Some(("AAGGTTAATACAAGCATCCCAACACTTCCACTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAGTGGAAGTGTTGGGATGCTTGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_75: TrimConfig = TrimConfig {
    kit_name: "NBD_75",
    end5: Some(("AAGGTTAAGACCATTGTGATGAACCCTGTTGTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGACAACAGGGTTCATCACAATGGTCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_76: TrimConfig = TrimConfig {
    kit_name: "NBD_76",
    end5: Some(("AAGGTTAAATGCTTGTTACATCAACCCTGGACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTCCAGGGTTGATGTAACAAGCATTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_77: TrimConfig = TrimConfig {
    kit_name: "NBD_77",
    end5: Some(("AAGGTTAACGACCTGTTTCTCAGGGATACAACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTTGTATCCCTGAGAAACAGGTCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_78: TrimConfig = TrimConfig {
    kit_name: "NBD_78",
    end5: Some(("AAGGTTAAAACAACCGAACCTTTGAATCAGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCTGATTCAAAGGTTCGGTTGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_79: TrimConfig = TrimConfig {
    kit_name: "NBD_79",
    end5: Some(("AAGGTTAATCTCGGAGATAGTTCTCACTGCTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAGCAGTGAGAACTATCTCCGAGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_80: TrimConfig = TrimConfig {
    kit_name: "NBD_80",
    end5: Some(("AAGGTTAACGGATGAACATAGGATAGCGATTCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGAATCGCTATCCTATGTTCATCCGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_81: TrimConfig = TrimConfig {
    kit_name: "NBD_81",
    end5: Some(("AAGGTTAACCTCATCTTGTGAAGTTGTTTCGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCGAAACAACTTCACAAGATGAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_82: TrimConfig = TrimConfig {
    kit_name: "NBD_82",
    end5: Some(("AAGGTTAAACGGTATGTCGAGTTCCAGGACTACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTAGTCCTGGAACTCGACATACCGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_83: TrimConfig = TrimConfig {
    kit_name: "NBD_83",
    end5: Some(("AAGGTTAATGGCTTGATCTAGGTAAGGTCGAACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTTCGACCTTACCTAGATCAAGCCATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_84: TrimConfig = TrimConfig {
    kit_name: "NBD_84",
    end5: Some(("AAGGTTAAGTAGTGGACCTAGAACCTGTGCCACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTGGCACAGGTTCTAGGTCCACTACTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_85: TrimConfig = TrimConfig {
    kit_name: "NBD_85",
    end5: Some(("AAGGTTAAAACGGAGGAGTTAGTTGGATGATCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGATCATCCAACTAACTCCTCCGTTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_86: TrimConfig = TrimConfig {
    kit_name: "NBD_86",
    end5: Some(("AAGGTTAAAGGTGATCCCAACAAGCGTAAGTACAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGTACTTACGCTTGTTGGGATCACCTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_87: TrimConfig = TrimConfig {
    kit_name: "NBD_87",
    end5: Some(("AAGGTTAATACATGCTCCTGTTGTTAGGGAGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCTCCCTAACAACAGGAGCATGTATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_88: TrimConfig = TrimConfig {
    kit_name: "NBD_88",
    end5: Some(("AAGGTTAATCTTCTACTACCGATCCGAAGCAGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCTGCTTCGGATCGGTAGTAGAAGATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_89: TrimConfig = TrimConfig {
    kit_name: "NBD_89",
    end5: Some(("AAGGTTAAACAGCATCAATGTTTGGCTAGTTGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCAACTAGCCAAACATTGATGCTGTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_90: TrimConfig = TrimConfig {
    kit_name: "NBD_90",
    end5: Some(("AAGGTTAAGATGTAGAGGGTACGGTTTGAGGCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGCCTCAAACCGTACCCTCTACATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_91: TrimConfig = TrimConfig {
    kit_name: "NBD_91",
    end5: Some(("AAGGTTAAGGCTCCATAGGAACTCACGCTACTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAGTAGCGTGAGTTCCTATGGAGCCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_92: TrimConfig = TrimConfig {
    kit_name: "NBD_92",
    end5: Some(("AAGGTTAATTGTGAGTGGAAAGATACAGGACCCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGGTCCTGTATCTTTCCACTCACAATTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_93: TrimConfig = TrimConfig {
    kit_name: "NBD_93",
    end5: Some(("AAGGTTAAAGTTTCCATCACTTCAGACTTGGGCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGCCCAAGTCTGAAGTGATGGAAACTTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_94: TrimConfig = TrimConfig {
    kit_name: "NBD_94",
    end5: Some(("AAGGTTAAGATTGTCCTCAAACTGCCACCTACCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGGTAGGTGGCAGTTTGAGGACAATCTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_95: TrimConfig = TrimConfig {
    kit_name: "NBD_95",
    end5: Some(("AAGGTTAACCTGTCTGGAAGAAGAATGGACTTCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGAAGTCCATTCTTCTTCCAGACAGGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
const NBD_96: TrimConfig = TrimConfig {
    kit_name: "NBD_96",
    end5: Some(("AAGGTTAACTGAACGGTCATAGAGTCCACCATCAGCACCT", NBD_END5)),
    end3: Some(("AGGTGCTGATGGTGGACTCTATGACCGTTCAGTTAACCTTAGCAAT", NBD_END3)),
    rev_com_end5: None,
    rev_com_end3: None,
};
/*
SQK-RBK114.24; SQK-RBK114.96
structure of reads with RA, with rapid barcode
Example for Rapid Barcode01
  |L_F             |Rapid Barcode01         |R_F                                               | insert Seq
5-GCTTGGGTGTTTAACC AAGAAAGTTGTCGGTGTCTTTGTG GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA .................-3
*/
const RBK_1: TrimConfig = TrimConfig {
    kit_name: "RBK_1",
    end5: Some((
        "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        // "AAGAAAGTTGTCGGTGTCTTTGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
        RBK_END5,
    )),
    end3: None,
    rev_com_end5: None,
    rev_com_end3: None,
};
// const RBK_2: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_2",
//     end5: Some((
//         "TCGATTCCGTTTGTAGTCGTCTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_3: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_3",
//     end5: Some((
//         "GAGTCTTGTGTCCCAGTTACCAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_4: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_4",
//     end5: Some((
//         "TTCGGATTCTATCGTGTTTCCCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_5: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_5",
//     end5: Some((
//         "CTTGTCCAGGGTTTGTGTAACCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_6: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_6",
//     end5: Some((
//         "TTCTCGCAAAGGCAGAAAGTAGTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_7: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_7",
//     end5: Some((
//         "GTGTTACCGTGGGAATGAATCCTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_8: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_8",
//     end5: Some((
//         "TTCAGGGAACAAACCAAGTTACGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_9: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_9",
//     end5: Some((
//         "AACTAGGCACAGCGAGTCTTGGTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_10: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_10",
//     end5: Some((
//         "AAGCGTTGAAACCTTTGTCCTCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_11: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_11",
//     end5: Some((
//         "GTTTCATCTATCGGAGGGAATGGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_12: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_12",
//     end5: Some((
//         "CAGGTAGAAAGAAGCAGAATCGGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_13: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_13",
//     end5: Some((
//         "AGAACGACTTCCATACTCGTGTGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_14: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_14",
//     end5: Some((
//         "AACGAGTCTCTTGGGACCCATAGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_15: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_15",
//     end5: Some((
//         "AGGTCTACCTCGCTAACACCACTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_16: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_16",
//     end5: Some((
//         "CGTCAACTGACAGTGGTTCGTACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_17: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_17",
//     end5: Some((
//         "ACCCTCCAGGAAAGTACCTCTGATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_18: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_18",
//     end5: Some((
//         "CCAAACCCAACAACCTAGATAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_19: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_19",
//     end5: Some((
//         "GTTCCTCGTGCAGTGTCAAGAGATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_20: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_20",
//     end5: Some((
//         "TTGCGTCCTGTTACGAGAACTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_21: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_21",
//     end5: Some((
//         "GAGCCTCTCATTGTCCGTTCTCTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_22: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_22",
//     end5: Some((
//         "ACCACTGCCATGTATCAAAGTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_23: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_23",
//     end5: Some((
//         "CTTACTACCCAGTGAACCTCCTCGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_24: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_24",
//     end5: Some((
//         "GCATAGTTCTGCATGATGGGTTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_25: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_25",
//     end5: Some((
//         "GTAAGTTGGGTATGCAACGCAATGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_26: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_26",
//     end5: Some((
//         "CATACAGCGACTACGCATTCTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_27: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_27",
//     end5: Some((
//         "CGACGGTTAGATTCACCTCTTACAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_28: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_28",
//     end5: Some((
//         "TGAAACCTAAGAAGGCACCGTATCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_29: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_29",
//     end5: Some((
//         "CTAGACACCTTGGGTTGACAGACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_30: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_30",
//     end5: Some((
//         "TCAGTGAGGATCTACTTCGACCCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_31: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_31",
//     end5: Some((
//         "TGCGTACAGCAATCAGTTACATTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_32: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_32",
//     end5: Some((
//         "CCAGTAGAAGTCCGACAACGTCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_33: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_33",
//     end5: Some((
//         "CAGACTTGGTACGGTTGGGTAACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_34: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_34",
//     end5: Some((
//         "GGACGAAGAACTCAAGTCAAAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_35: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_35",
//     end5: Some((
//         "CTACTTACGAAGCTGAGGGACTGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_36: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_36",
//     end5: Some((
//         "ATGTCCCAGTTAGAGGAGGAAACAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_37: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_37",
//     end5: Some((
//         "GCTTGCGATTGATGCTTAGTATCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_38: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_38",
//     end5: Some((
//         "ACCACAGGAGGACGATACAGAGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_39: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_39",
//     end5: Some((
//         "CCACAGTGTCAACTAGAGCCTCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_40: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_40",
//     end5: Some((
//         "TAGTTTGGATGACCAAGGATAGCCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_41: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_41",
//     end5: Some((
//         "GGAGTTCGTCCAGAGAAGTACACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_42: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_42",
//     end5: Some((
//         "CTACGTGTAAGGCATACCTGCCAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_43: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_43",
//     end5: Some((
//         "CTTTCGTTGTTGACTCGACGGTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_44: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_44",
//     end5: Some((
//         "AGTAGAAAGGGTTCCTTCCCACTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_45: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_45",
//     end5: Some((
//         "GATCCAACAGAGATGCCTTCAGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_46: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_46",
//     end5: Some((
//         "GCTGTGTTCCACTTCATTCTCCTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_47: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_47",
//     end5: Some((
//         "GTGCAACTTTCCCACAGGTAGTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_48: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_48",
//     end5: Some((
//         "CATCTGGAACGTGGTACACCTGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_49: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_49",
//     end5: Some((
//         "ACTGGTGCAGCTTTGAACATCTAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_50: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_50",
//     end5: Some((
//         "ATGGACTTTGGTAACTTCCTGCGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_51: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_51",
//     end5: Some((
//         "GTTGAATGAGCCTACTGGGTCCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_52: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_52",
//     end5: Some((
//         "TGAGAGACAAGATTGTTCGTGGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_53: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_53",
//     end5: Some((
//         "AGATTCAGACCGTCTCATGCAAAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_54: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_54",
//     end5: Some((
//         "CAAGAGCTTTGACTAAGGAGCATGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_55: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_55",
//     end5: Some((
//         "TGGAAGATGAGACCCTGATCTACGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_56: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_56",
//     end5: Some((
//         "TCACTACTCAACAGGTGGCATGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_57: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_57",
//     end5: Some((
//         "GCTAGGTCAATCTCCTTCGGAAGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_58: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_58",
//     end5: Some((
//         "CAGGTTACTCCTCCGTGAGTCTGAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_59: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_59",
//     end5: Some((
//         "TCAATCAAGAAGGGAAAGCAAGGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_60: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_60",
//     end5: Some((
//         "CATGTTCAACCAAGGCTTCTATGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_61: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_61",
//     end5: Some((
//         "AGAGGGTACTATGTGCCTCAGCACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_62: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_62",
//     end5: Some((
//         "CACCCACACTTACTTCAGGACGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_63: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_63",
//     end5: Some((
//         "TTCTGAAGTTCCTGGGTCTTGAACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_64: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_64",
//     end5: Some((
//         "GACAGACACCGTTCATCGACTTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_65: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_65",
//     end5: Some((
//         "TTCTCAGTCTTCCTCCAGACAAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_66: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_66",
//     end5: Some((
//         "CCGATCCTTGTGGCTTCTAACTTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_67: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_67",
//     end5: Some((
//         "GTTTGTCATACTCGTGTGCTCACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_68: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_68",
//     end5: Some((
//         "GAATCTAAGCAAACACGAAGGTGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_69: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_69",
//     end5: Some((
//         "TACAGTCCGAGCCTCATGTGATCTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_70: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_70",
//     end5: Some((
//         "ACCGAGATCCTACGAATGGAGTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_71: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_71",
//     end5: Some((
//         "CCTGGGAGCATCAGGTAGTAACAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_72: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_72",
//     end5: Some((
//         "TAGCTGACTGTCTTCCATACCGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_73: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_73",
//     end5: Some((
//         "AAGAAACAGGATGACAGAACCCTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_74: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_74",
//     end5: Some((
//         "TACAAGCATCCCAACACTTCCACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_75: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_75",
//     end5: Some((
//         "GACCATTGTGATGAACCCTGTTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_76: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_76",
//     end5: Some((
//         "ATGCTTGTTACATCAACCCTGGACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_77: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_77",
//     end5: Some((
//         "CGACCTGTTTCTCAGGGATACAACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_78: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_78",
//     end5: Some((
//         "AACAACCGAACCTTTGAATCAGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_79: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_79",
//     end5: Some((
//         "TCTCGGAGATAGTTCTCACTGCTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_80: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_80",
//     end5: Some((
//         "CGGATGAACATAGGATAGCGATTCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_81: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_81",
//     end5: Some((
//         "CCTCATCTTGTGAAGTTGTTTCGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_82: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_82",
//     end5: Some((
//         "ACGGTATGTCGAGTTCCAGGACTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_83: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_83",
//     end5: Some((
//         "TGGCTTGATCTAGGTAAGGTCGAAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_84: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_84",
//     end5: Some((
//         "GTAGTGGACCTAGAACCTGTGCCAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_85: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_85",
//     end5: Some((
//         "AACGGAGGAGTTAGTTGGATGATCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_86: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_86",
//     end5: Some((
//         "AGGTGATCCCAACAAGCGTAAGTAGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_87: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_87",
//     end5: Some((
//         "TACATGCTCCTGTTGTTAGGGAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_88: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_88",
//     end5: Some((
//         "TCTTCTACTACCGATCCGAAGCAGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_89: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_89",
//     end5: Some((
//         "ACAGCATCAATGTTTGGCTAGTTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_90: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_90",
//     end5: Some((
//         "GATGTAGAGGGTACGGTTTGAGGCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_91: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_91",
//     end5: Some((
//         "GGCTCCATAGGAACTCACGCTACTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_92: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_92",
//     end5: Some((
//         "TTGTGAGTGGAAAGATACAGGACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_93: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_93",
//     end5: Some((
//         "AGTTTCCATCACTTCAGACTTGGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_94: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_94",
//     end5: Some((
//         "GATTGTCCTCAAACTGCCACCTACGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_95: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_95",
//     end5: Some((
//         "CCTGTCTGGAAGAAGAATGGACTTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };
// const RBK_96: SequenceInfo = SequenceInfo {
//     kit_name: "RBK_96",
//     end5: Some((
//         "CTGAACGGTCATAGAGTCCACCATGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
//         RBK_END5,
//     )),
//     end3: None,
//     rev_com_end5: None,
//     rev_com_end3: None,
// };

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
const PCS: TrimConfig = TrimConfig {
    kit_name: "PCS",
    end5: Some((
        "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCS_END5,
    )),
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
const PCB: TrimConfig = TrimConfig {
    kit_name: "PCB_1",
    end5: Some((
        "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
        PCB_END5,
    )),
    end3: Some(("CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT", PCB_END3)),
    rev_com_end5: Some((
        "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
        PCB_REV_COM_END5,
    )),
    rev_com_end3: Some((
        "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA",
        PCB_REV_COM_END3,
    )),
};
// const PCB_2: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_2",
//     end5: Some((
//         "TCGATTCCGTTTGTAGTCGTCTGTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTACAGACGACTACAAACGGAATCGA",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "TCGATTCCGTTTGTAGTCGTCTGTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAACAGACGACTACAAACGGAATCGA",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_3: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_3",
//     end5: Some((
//         "GAGTCTTGTGTCCCAGTTACCAGGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCCTGGTAACTGGGACACAAGACTC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GAGTCTTGTGTCCCAGTTACCAGGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACCTGGTAACTGGGACACAAGACTC",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_4: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_4",
//     end5: Some((
//         "TTCGGATTCTATCGTGTTTCCCTATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTAGGGAAACACGATAGAATCCGAA",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "TTCGGATTCTATCGTGTTTCCCTAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATAGGGAAACACGATAGAATCCGAA",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_5: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_5",
//     end5: Some((
//         "CTTGTCCAGGGTTTGTGTAACCTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAAGGTTACACAAACCCTGGACAAG",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "CTTGTCCAGGGTTTGTGTAACCTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAAGGTTACACAAACCCTGGACAAG",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_6: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_6",
//     end5: Some((
//         "TTCTCGCAAAGGCAGAAAGTAGTCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGACTACTTTCTGCCTTTGCGAGAA",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "TTCTCGCAAAGGCAGAAAGTAGTCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGACTACTTTCTGCCTTTGCGAGAA",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_7: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_7",
//     end5: Some((
//         "GTGTTACCGTGGGAATGAATCCTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAAGGATTCATTCCCACGGTAACAC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GTGTTACCGTGGGAATGAATCCTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAAGGATTCATTCCCACGGTAACAC",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_8: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_8",
//     end5: Some((
//         "TTCAGGGAACAAACCAAGTTACGTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTACGTAACTTGGTTTGTTCCCTGAA",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "TTCAGGGAACAAACCAAGTTACGTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAACGTAACTTGGTTTGTTCCCTGAA",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_9: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_9",
//     end5: Some((
//         "AACTAGGCACAGCGAGTCTTGGTTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAACCAAGACTCGCTGTGCCTAGTT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "AACTAGGCACAGCGAGTCTTGGTTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAACCAAGACTCGCTGTGCCTAGTT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_10: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_10",
//     end5: Some((
//         "AAGCGTTGAAACCTTTGTCCTCTCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGAGAGGACAAAGGTTTCAACGCTT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "AAGCGTTGAAACCTTTGTCCTCTCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGAGAGGACAAAGGTTTCAACGCTT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_11: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_11",
//     end5: Some((
//         "GTTTCATCTATCGGAGGGAATGGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCCATTCCCTCCGATAGATGAAAC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GTTTCATCTATCGGAGGGAATGGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCCATTCCCTCCGATAGATGAAAC",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_12: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_12",
//     end5: Some((
//         "CAGGTAGAAAGAAGCAGAATCGGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCCGATTCTGCTTCTTTCTACCTG",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "CAGGTAGAAAGAAGCAGAATCGGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCCGATTCTGCTTCTTTCTACCTG",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_13: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_13",
//     end5: Some((
//         "AGAACGACTTCCATACTCGTGTGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCACACGAGTATGGAAGTCGTTCT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "AGAACGACTTCCATACTCGTGTGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCACACGAGTATGGAAGTCGTTCT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_14: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_14",
//     end5: Some((
//         "AACGAGTCTCTTGGGACCCATAGATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTCTATGGGTCCCAAGAGACTCGTT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "AACGAGTCTCTTGGGACCCATAGAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATCTATGGGTCCCAAGAGACTCGTT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_15: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_15",
//     end5: Some((
//         "AGGTCTACCTCGCTAACACCACTGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCAGTGGTGTTAGCGAGGTAGACCT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "AGGTCTACCTCGCTAACACCACTGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACAGTGGTGTTAGCGAGGTAGACCT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_16: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_16",
//     end5: Some((
//         "CGTCAACTGACAGTGGTTCGTACTTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTAGTACGAACCACTGTCAGTTGACG",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "CGTCAACTGACAGTGGTTCGTACTACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAAGTACGAACCACTGTCAGTTGACG",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_17: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_17",
//     end5: Some((
//         "ACCCTCCAGGAAAGTACCTCTGATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATCAGAGGTACTTTCCTGGAGGGT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "ACCCTCCAGGAAAGTACCTCTGATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATCAGAGGTACTTTCCTGGAGGGT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_18: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_18",
//     end5: Some((
//         "CCAAACCCAACAACCTAGATAGGCTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTGCCTATCTAGGTTGTTGGGTTTGG",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "CCAAACCCAACAACCTAGATAGGCACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAGCCTATCTAGGTTGTTGGGTTTGG",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_19: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_19",
//     end5: Some((
//         "GTTCCTCGTGCAGTGTCAAGAGATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATCTCTTGACACTGCACGAGGAAC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GTTCCTCGTGCAGTGTCAAGAGATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATCTCTTGACACTGCACGAGGAAC",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_20: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_20",
//     end5: Some((
//         "TTGCGTCCTGTTACGAGAACTCATTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTATGAGTTCTCGTAACAGGACGCAA",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "TTGCGTCCTGTTACGAGAACTCATACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAAATGAGTTCTCGTAACAGGACGCAA",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_21: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_21",
//     end5: Some((
//         "GAGCCTCTCATTGTCCGTTCTCTATTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTTAGAGAACGGACAATGAGAGGCTC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GAGCCTCTCATTGTCCGTTCTCTAACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAATAGAGAACGGACAATGAGAGGCTC",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_22: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_22",
//     end5: Some((
//         "ACCACTGCCATGTATCAAAGTACGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCGTACTTTGATACATGGCAGTGGT",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "ACCACTGCCATGTATCAAAGTACGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACGTACTTTGATACATGGCAGTGGT",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_23: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_23",
//     end5: Some((
//         "CTTACTACCCAGTGAACCTCCTCGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCGAGGAGGTTCACTGGGTAGTAAG",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "CTTACTACCCAGTGAACCTCCTCGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACGAGGAGGTTCACTGGGTAGTAAG",
//         PCB_REV_COM_END3,
//     )),
// };
// const PCB_24: SequenceInfo = SequenceInfo {
//     kit_name: "PCB_24",
//     end5: Some((
//         "GCATAGTTCTGCATGATGGGTTAGTTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG",
//         PCB_END5,
//     )),
//     end3: Some((
//         "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGTCTAACCCATCATGCAGAACTATGC",
//         PCB_END3,
//     )),
//     rev_com_end5: Some((
//         "GCATAGTTCTGCATGATGGGTTAGACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG",
//         PCB_REV_COM_END5,
//     )),
//     rev_com_end3: Some((
//         "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAACTAACCCATCATGCAGAACTATGC",
//         PCB_REV_COM_END3,
//     )),
// };
