#ifndef NANOFQ_ADAPTER_H
#define NANOFQ_ADAPTER_H

#include <string>
#include <unordered_map>
#include <fmt/core.h>
#include "myUtility.h"
#include "SequenceInfo.h"

// LSK114
constexpr trim_end LSK_TOP5END{100, 0.75, 0.75}; // 28
constexpr trim_end LSK_TOP3END{60, 0.5, 0.75}; // 28

// NBD114
constexpr trim_end NBD_TOP5END{150, 0.4, 0.75}; // 68
constexpr trim_end NBD_TOP3END{150, 0.3, 0.75}; // 68

// RAD114
constexpr trim_end RAD_TOP5END{150, 0.5, 0.75}; // 66

// RBK114
constexpr trim_end RBK_TOP5END{200, 0.4, 0.75}; // 90

// PCS114 TOP
constexpr trim_end PCS_TOP5END{150, 0.6, 0.75}; // 53
constexpr trim_end PCS_TOP3END{150, 0.4, 0.75}; // 45
constexpr trim_end PCS_BOT5END{150, 0.6, 0.75}; // 45
constexpr trim_end PCS_BOT3END{150, 0.4, 0.75}; // 53

// PCB114 TOP
// constexpr trim_end PCB_TOP5END{180, 0.5, 0.75}; // 77
// constexpr trim_end PCB_TOP3END{180, 0.4, 0.75}; // 70
// constexpr trim_end PCB_BOT5END{180, 0.5, 0.75}; // 70
// constexpr trim_end PCB_BOT3END{180, 0.4, 0.75}; // 77
constexpr trim_end PCB_TOP5END{180, 0.3, 0.75}; // 77
constexpr trim_end PCB_TOP3END{180, 0.3, 0.75}; // 70
constexpr trim_end PCB_BOT5END{180, 0.3, 0.75}; // 70
constexpr trim_end PCB_BOT3END{180, 0.3, 0.75}; // 77

// https://nanoporetech.com/document/chemistry-technical-document#adapter-sequences
namespace barcode_info
{
    /*
    * SQK-LSK114
    * LSK114 library reads structure
    *           |--->  LA_ADAPTER_5   <----| | insert Seq | |--->   LA_ADAPTER_3   <---|
    * 5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-..............-AGCAATACGTAACTGAACGAAGTACAGG-3
    *
    * 3' end always is truncated
    */
    const std::string LA_ADAPTER_5 = "CCTGTACTTCGTTCAGTTACGTATTGCT";
    const std::string LA_ADAPTER_3 = "AGCAATACGTAACTGAACGAAGTACAGG";
    // the first 15 or even less bases is enough, because 3'end always is truncated


    /*
     * SQK-NBD114-24; SQK-NBD114-96
     * NBD114-24/96 library reads structure
     * Example for Native Barcode01
     *           |NA_ADAPTER_5                |L_F_5   |Barcode01 rev com       |R_F_5   |insert Seq         |L_F_3   |Barcode01               |R_F_3         |NA_ADAPTER_3
     * 5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT AAGGTTAA CACAAAGACACCGACAACTTTCTT CAGCACCT ................... AGGTGCTG AAGAAAGTTGTCGGTGTCTTTGTG TTAACCTTAGCAAT ACGTAACTGAACGAAGTACAGG-3
    */
    const std::string NA_ADAPTER_5 = "CCTGTACTTCGTTCAGTTACGTATTGCT";
    const std::string NA_ADAPTER_3 = "ACGTAACTGAACGAAGTACAGG";
    const std::string NB_LEFT_FLANKING_5 = "AAGGTTAA"; // the left flanking seq of barcode at 5' end
    const std::string NB_RIGHT_FLANKING_5 = "CAGCACCT"; // the right flanking seq of barcode at 3' end
    const std::string NB_LEFT_FLANKING_3 = "AGGTGCTG";
    const std::string NB_RIGHT_FLANKING_3 = "TTAACCTTAGCAAT";

    /*
     * SQK-RAD114; SQK-ULK114
     * the rapid adapter(RA) from ont document is 5’-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-3', but RA_ADAPTER will be used when trimming reads with Rapid Adapter(RA)
     * Always consider only the adapter at 5' end for Rapid library
     *
     * structure of reads with RA, but no barcode
     *   |RA_ADAPTER we want to trim from reads                             | insert Seq
     * 5-GCTTGGGTGTTTAACCGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA-...................-3
     *
     *
     * SQK-RBK114.24; SQK-RBK114.96
     * structure of reads with RA, with rapid barcode
     * Example for Rapid Barcode01
     *   |L_F             |Rapid Barcode01         |R_F                                               | insert Seq
     * 5-GCTTGGGTGTTTAACC AAGAAAGTTGTCGGTGTCTTTGTG GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA .................-3
    */
    const std::string RB_LEFT_FLANKING = "GCTTGGGTGTTTAACC";
    const std::string RB_RIGTH_FLANKING = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";
    const std::string RA_ADAPTER = RB_LEFT_FLANKING + RB_RIGTH_FLANKING;

    /*
     * cDNA-PCR Sequencing Kit: SQK-PCS114
     * cDNA-PCR Barcoding Kit V14: SQK-PCB114.24
     * Strand Switching Primer II (SSPII)
     * cDNA RT Adapter (CRTA)
     *
     * SQK-PCS114 structure
     *      |SSPII                                                | insert Seq with polyA  | CRTA
     * 5-...TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG .........AAAAAAAAAAAAAAA CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG...-3
     * 3-...CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA .........TTTTTTTTTTTTTTT GAACGCCCGCCGCCTGAGAGGAGACTTCTATCTCGCTGTCCGTTC...-5
     *
     *
     * SQK-PCB114.24 structure
     *      | BP01                   | SSPII                                               | insert Seqw with polyA | CRTA                                         | BP01 reverse com
     * 5-...AAGAAAGTTGTCGGTGTCTTTGTG TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG .........AAAAAAAAAAAAAAA CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT CACAAAGACACCGACAACTTTCTT...-3
     * 3-...TTCTTTCAACAGCCACAGAAACAC AAAGACAACCACGACTATAACGAAABBBBAABBBBAABBBBAABBBBAAACCC .........TTTTTTTTTTTTTTT GAACGCCCGCCGCCTGAGAGGAGACTTCTATCTCGCTGTCCGTTCA GTGTTTCTGTGGCTGTTGAAAGAA...-5
*/
    const std::string SSPII = "TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG"; //  V: [GCA]
    const std::string CRTA = "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAGT"; // more one T than document at 3' end
    const std::string CRTA_REV_COM = "ACTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG";
    const std::string SSPII_REV_COM = "CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA";

    /*BARCODE*/
    /*Native Barcode*/
    const std::string NB01 = "AAGAAAGTTGTCGGTGTCTTTGTG";
    const std::string NB02 = "TCGATTCCGTTTGTAGTCGTCTGT";
    const std::string NB03 = "GAGTCTTGTGTCCCAGTTACCAGG";
    const std::string NB04 = "TTCGGATTCTATCGTGTTTCCCTA";
    const std::string NB05 = "CTTGTCCAGGGTTTGTGTAACCTT";
    const std::string NB06 = "TTCTCGCAAAGGCAGAAAGTAGTC";
    const std::string NB07 = "GTGTTACCGTGGGAATGAATCCTT";
    const std::string NB08 = "TTCAGGGAACAAACCAAGTTACGT";
    const std::string NB09 = "AACTAGGCACAGCGAGTCTTGGTT";
    const std::string NB10 = "AAGCGTTGAAACCTTTGTCCTCTC";
    const std::string NB11 = "GTTTCATCTATCGGAGGGAATGGA";
    const std::string NB12 = "CAGGTAGAAAGAAGCAGAATCGGA";
    const std::string NB13 = "TCACACGAGTATGGAAGTCGTTCT";
    const std::string NB14 = "TCTATGGGTCCCAAGAGACTCGTT";
    const std::string NB15 = "CAGTGGTGTTAGCGAGGTAGACCT";
    const std::string NB16 = "AGTACGAACCACTGTCAGTTGACG";
    const std::string NB17 = "ATCAGAGGTACTTTCCTGGAGGGT";
    const std::string NB18 = "GCCTATCTAGGTTGTTGGGTTTGG";
    const std::string NB19 = "ATCTCTTGACACTGCACGAGGAAC";
    const std::string NB20 = "ATGAGTTCTCGTAACAGGACGCAA";
    const std::string NB21 = "TAGAGAACGGACAATGAGAGGCTC";
    const std::string NB22 = "CGTACTTTGATACATGGCAGTGGT";
    const std::string NB23 = "CGAGGAGGTTCACTGGGTAGTAAG";
    const std::string NB24 = "CTAACCCATCATGCAGAACTATGC";
    const std::string NB25 = "CATTGCGTTGCATACCCAACTTAC";
    const std::string NB26 = "ATGAGAATGCGTAGTCGCTGTATG";
    const std::string NB27 = "TGTAAGAGGTGAATCTAACCGTCG";
    const std::string NB28 = "GATACGGTGCCTTCTTAGGTTTCA";
    const std::string NB29 = "GGTCTGTCAACCCAAGGTGTCTAG";
    const std::string NB30 = "TGGGTCGAAGTAGATCCTCACTGA";
    const std::string NB31 = "CAATGTAACTGATTGCTGTACGCA";
    const std::string NB32 = "ATGACGTTGTCGGACTTCTACTGG";
    const std::string NB33 = "AGTTACCCAACCGTACCAAGTCTG";
    const std::string NB34 = "GCCTTTGACTTGAGTTCTTCGTCC";
    const std::string NB35 = "GCAGTCCCTCAGCTTCGTAAGTAG";
    const std::string NB36 = "TGTTTCCTCCTCTAACTGGGACAT";
    const std::string NB37 = "TGATACTAAGCATCAATCGCAAGC";
    const std::string NB38 = "TTCTCTGTATCGTCCTCCTGTGGT";
    const std::string NB39 = "GAGAGGCTCTAGTTGACACTGTGG";
    const std::string NB40 = "GGCTATCCTTGGTCATCCAAACTA";
    const std::string NB41 = "CGTGTACTTCTCTGGACGAACTCC";
    const std::string NB42 = "CTGGCAGGTATGCCTTACACGTAG";
    const std::string NB43 = "CTACCGTCGAGTCAACAACGAAAG";
    const std::string NB44 = "GAGTGGGAAGGAACCCTTTCTACT";
    const std::string NB45 = "CACTGAAGGCATCTCTGTTGGATC";
    const std::string NB46 = "CAGGAGAATGAAGTGGAACACAGC";
    const std::string NB47 = "GAACTACCTGTGGGAAAGTTGCAC";
    const std::string NB48 = "TACAGGTGTACCACGTTCCAGATG";
    const std::string NB49 = "CTAGATGTTCAAAGCTGCACCAGT";
    const std::string NB50 = "ACGCAGGAAGTTACCAAAGTCCAT";
    const std::string NB51 = "GAGGACCCAGTAGGCTCATTCAAC";
    const std::string NB52 = "GTCCACGAACAATCTTGTCTCTCA";
    const std::string NB53 = "CTTTGCATGAGACGGTCTGAATCT";
    const std::string NB54 = "CATGCTCCTTAGTCAAAGCTCTTG";
    const std::string NB55 = "CGTAGATCAGGGTCTCATCTTCCA";
    const std::string NB56 = "TTCATGCCACCTGTTGAGTAGTGA";
    const std::string NB57 = "ACTTCCGAAGGAGATTGACCTAGC";
    const std::string NB58 = "TCAGACTCACGGAGGAGTAACCTG";
    const std::string NB59 = "ACCTTGCTTTCCCTTCTTGATTGA";
    const std::string NB60 = "CCATAGAAGCCTTGGTTGAACATG";
    const std::string NB61 = "GTGCTGAGGCACATAGTACCCTCT";
    const std::string NB62 = "TACGTCCTGAAGTAAGTGTGGGTG";
    const std::string NB63 = "GTTCAAGACCCAGGAACTTCAGAA";
    const std::string NB64 = "GAAAGTCGATGAACGGTGTCTGTC";
    const std::string NB65 = "CCTTGTCTGGAGGAAGACTGAGAA";
    const std::string NB66 = "GAAGTTAGAAGCCACAAGGATCGG";
    const std::string NB67 = "GGTGAGCACACGAGTATGACAAAC";
    const std::string NB68 = "CCACCTTCGTGTTTGCTTAGATTC";
    const std::string NB69 = "AGATCACATGAGGCTCGGACTGTA";
    const std::string NB70 = "ACACTCCATTCGTAGGATCTCGGT";
    const std::string NB71 = "CTGTTACTACCTGATGCTCCCAGG";
    const std::string NB72 = "GTCGGTATGGAAGACAGTCAGCTA";
    const std::string NB73 = "GAGGGTTCTGTCATCCTGTTTCTT";
    const std::string NB74 = "AGTGGAAGTGTTGGGATGCTTGTA";
    const std::string NB75 = "ACAACAGGGTTCATCACAATGGTC";
    const std::string NB76 = "GTCCAGGGTTGATGTAACAAGCAT";
    const std::string NB77 = "GTTGTATCCCTGAGAAACAGGTCG";
    const std::string NB78 = "TTCTGATTCAAAGGTTCGGTTGTT";
    const std::string NB79 = "CAGCAGTGAGAACTATCTCCGAGA";
    const std::string NB80 = "GAATCGCTATCCTATGTTCATCCG";
    const std::string NB81 = "CCGAAACAACTTCACAAGATGAGG";
    const std::string NB82 = "TAGTCCTGGAACTCGACATACCGT";
    const std::string NB83 = "TTCGACCTTACCTAGATCAAGCCA";
    const std::string NB84 = "TGGCACAGGTTCTAGGTCCACTAC";
    const std::string NB85 = "GATCATCCAACTAACTCCTCCGTT";
    const std::string NB86 = "TACTTACGCTTGTTGGGATCACCT";
    const std::string NB87 = "CCTCCCTAACAACAGGAGCATGTA";
    const std::string NB88 = "CTGCTTCGGATCGGTAGTAGAAGA";
    const std::string NB89 = "CAACTAGCCAAACATTGATGCTGT";
    const std::string NB90 = "GCCTCAAACCGTACCCTCTACATC";
    const std::string NB91 = "AGTAGCGTGAGTTCCTATGGAGCC";
    const std::string NB92 = "GGTCCTGTATCTTTCCACTCACAA";
    const std::string NB93 = "CCCAAGTCTGAAGTGATGGAAACT";
    const std::string NB94 = "GTAGGTGGCAGTTTGAGGACAATC";
    const std::string NB95 = "AAGTCCATTCTTCTTCCAGACAGG";
    const std::string NB96 = "ATGGTGGACTCTATGACCGTTCAG";

    /*Rapid Barcode*/
    const std::string RB01 = NB01;
    const std::string RB02 = NB02;
    const std::string RB03 = NB03;
    const std::string RB04 = NB04;
    const std::string RB05 = NB05;
    const std::string RB06 = NB06;
    const std::string RB07 = NB07;
    const std::string RB08 = NB08;
    const std::string RB09 = NB09;
    const std::string RB10 = NB10;
    const std::string RB11 = NB11;
    const std::string RB12 = NB12;
    const std::string RB13 = "AGAACGACTTCCATACTCGTGTGA";
    const std::string RB14 = "AACGAGTCTCTTGGGACCCATAGA";
    const std::string RB15 = "AGGTCTACCTCGCTAACACCACTG";
    const std::string RB16 = "CGTCAACTGACAGTGGTTCGTACT";
    const std::string RB17 = "ACCCTCCAGGAAAGTACCTCTGAT";
    const std::string RB18 = "CCAAACCCAACAACCTAGATAGGC";
    const std::string RB19 = "GTTCCTCGTGCAGTGTCAAGAGAT";
    const std::string RB20 = "TTGCGTCCTGTTACGAGAACTCAT";
    const std::string RB21 = "GAGCCTCTCATTGTCCGTTCTCTA";
    const std::string RB22 = "ACCACTGCCATGTATCAAAGTACG";
    const std::string RB23 = "CTTACTACCCAGTGAACCTCCTCG";
    const std::string RB24 = "GCATAGTTCTGCATGATGGGTTAG";
    const std::string RB25 = "GTAAGTTGGGTATGCAACGCAATG";
    const std::string RB26 = "CATACAGCGACTACGCATTCTCAT";
    const std::string RB27 = "CGACGGTTAGATTCACCTCTTACA";
    const std::string RB28 = "TGAAACCTAAGAAGGCACCGTATC";
    const std::string RB29 = "CTAGACACCTTGGGTTGACAGACC";
    const std::string RB30 = "TCAGTGAGGATCTACTTCGACCCA";
    const std::string RB31 = "TGCGTACAGCAATCAGTTACATTG";
    const std::string RB32 = "CCAGTAGAAGTCCGACAACGTCAT";
    const std::string RB33 = "CAGACTTGGTACGGTTGGGTAACT";
    const std::string RB34 = "GGACGAAGAACTCAAGTCAAAGGC";
    const std::string RB35 = "CTACTTACGAAGCTGAGGGACTGC";
    const std::string RB36 = "ATGTCCCAGTTAGAGGAGGAAACA";
    const std::string RB37 = "GCTTGCGATTGATGCTTAGTATCA";
    const std::string RB38 = "ACCACAGGAGGACGATACAGAGAA";
    const std::string RB39 = "CCACAGTGTCAACTAGAGCCTCTC";
    const std::string RB40 = "TAGTTTGGATGACCAAGGATAGCC";
    const std::string RB41 = "GGAGTTCGTCCAGAGAAGTACACG";
    const std::string RB42 = "CTACGTGTAAGGCATACCTGCCAG";
    const std::string RB43 = "CTTTCGTTGTTGACTCGACGGTAG";
    const std::string RB44 = "AGTAGAAAGGGTTCCTTCCCACTC";
    const std::string RB45 = "GATCCAACAGAGATGCCTTCAGTG";
    const std::string RB46 = "GCTGTGTTCCACTTCATTCTCCTG";
    const std::string RB47 = "GTGCAACTTTCCCACAGGTAGTTC";
    const std::string RB48 = "CATCTGGAACGTGGTACACCTGTA";
    const std::string RB49 = "ACTGGTGCAGCTTTGAACATCTAG";
    const std::string RB50 = "ATGGACTTTGGTAACTTCCTGCGT";
    const std::string RB51 = "GTTGAATGAGCCTACTGGGTCCTC";
    const std::string RB52 = "TGAGAGACAAGATTGTTCGTGGAC";
    const std::string RB53 = "AGATTCAGACCGTCTCATGCAAAG";
    const std::string RB54 = "CAAGAGCTTTGACTAAGGAGCATG";
    const std::string RB55 = "TGGAAGATGAGACCCTGATCTACG";
    const std::string RB56 = "TCACTACTCAACAGGTGGCATGAA";
    const std::string RB57 = "GCTAGGTCAATCTCCTTCGGAAGT";
    const std::string RB58 = "CAGGTTACTCCTCCGTGAGTCTGA";
    const std::string RB59 = "TCAATCAAGAAGGGAAAGCAAGGT";
    const std::string RB60 = "CATGTTCAACCAAGGCTTCTATGG";
    const std::string RB61 = "AGAGGGTACTATGTGCCTCAGCAC";
    const std::string RB62 = "CACCCACACTTACTTCAGGACGTA";
    const std::string RB63 = "TTCTGAAGTTCCTGGGTCTTGAAC";
    const std::string RB64 = "GACAGACACCGTTCATCGACTTTC";
    const std::string RB65 = "TTCTCAGTCTTCCTCCAGACAAGG";
    const std::string RB66 = "CCGATCCTTGTGGCTTCTAACTTC";
    const std::string RB67 = "GTTTGTCATACTCGTGTGCTCACC";
    const std::string RB68 = "GAATCTAAGCAAACACGAAGGTGG";
    const std::string RB69 = "TACAGTCCGAGCCTCATGTGATCT";
    const std::string RB70 = "ACCGAGATCCTACGAATGGAGTGT";
    const std::string RB71 = "CCTGGGAGCATCAGGTAGTAACAG";
    const std::string RB72 = "TAGCTGACTGTCTTCCATACCGAC";
    const std::string RB73 = "AAGAAACAGGATGACAGAACCCTC";
    const std::string RB74 = "TACAAGCATCCCAACACTTCCACT";
    const std::string RB75 = "GACCATTGTGATGAACCCTGTTGT";
    const std::string RB76 = "ATGCTTGTTACATCAACCCTGGAC";
    const std::string RB77 = "CGACCTGTTTCTCAGGGATACAAC";
    const std::string RB78 = "AACAACCGAACCTTTGAATCAGAA";
    const std::string RB79 = "TCTCGGAGATAGTTCTCACTGCTG";
    const std::string RB80 = "CGGATGAACATAGGATAGCGATTC";
    const std::string RB81 = "CCTCATCTTGTGAAGTTGTTTCGG";
    const std::string RB82 = "ACGGTATGTCGAGTTCCAGGACTA";
    const std::string RB83 = "TGGCTTGATCTAGGTAAGGTCGAA";
    const std::string RB84 = "GTAGTGGACCTAGAACCTGTGCCA";
    const std::string RB85 = "AACGGAGGAGTTAGTTGGATGATC";
    const std::string RB86 = "AGGTGATCCCAACAAGCGTAAGTA";
    const std::string RB87 = "TACATGCTCCTGTTGTTAGGGAGG";
    const std::string RB88 = "TCTTCTACTACCGATCCGAAGCAG";
    const std::string RB89 = "ACAGCATCAATGTTTGGCTAGTTG";
    const std::string RB90 = "GATGTAGAGGGTACGGTTTGAGGC";
    const std::string RB91 = "GGCTCCATAGGAACTCACGCTACT";
    const std::string RB92 = "TTGTGAGTGGAAAGATACAGGACC";
    const std::string RB93 = "AGTTTCCATCACTTCAGACTTGGG";
    const std::string RB94 = "GATTGTCCTCAAACTGCCACCTAC";
    const std::string RB95 = "CCTGTCTGGAAGAAGAATGGACTT";
    const std::string RB96 = "CTGAACGGTCATAGAGTCCACCAT";
    const std::vector<std::string> NB_VEC{
        NB01, NB02, NB03, NB04, NB05, NB06, NB07, NB08, NB09, NB10, NB11, NB12,
        NB13, NB14, NB15, NB16, NB17, NB18, NB19, NB20, NB21, NB22, NB23, NB24,
        NB25, NB26, NB27, NB28, NB29, NB30, NB31, NB32, NB33, NB34, NB35, NB36,
        NB37, NB38, NB39, NB40, NB41, NB42, NB43, NB44, NB45, NB46, NB47, NB48,
        NB49, NB50, NB51, NB52, NB53, NB54, NB55, NB56, NB57, NB58, NB59, NB60,
        NB61, NB62, NB63, NB64, NB65, NB66, NB67, NB68, NB69, NB70, NB71, NB72,
        NB73, NB74, NB75, NB76, NB77, NB78, NB79, NB80, NB81, NB82, NB83, NB84,
        NB85, NB86, NB87, NB88, NB89, NB90, NB91, NB92, NB93, NB94, NB95, NB96
    };

    const std::vector<std::string> RB_VEC{
        RB01, RB02, RB03, RB04, RB05, RB06, RB07, RB08, RB09, RB10, RB11, RB12,
        RB13, RB14, RB15, RB16, RB17, RB18, RB19, RB20, RB21, RB22, RB23, RB24,
        RB25, RB26, RB27, RB28, RB29, RB30, RB31, RB32, RB33, RB34, RB35, RB36,
        RB37, RB38, RB39, RB40, RB41, RB42, RB43, RB44, RB45, RB46, RB47, RB48,
        RB49, RB50, RB51, RB52, RB53, RB54, RB55, RB56, RB57, RB58, RB59, RB60,
        RB61, RB62, RB63, RB64, RB65, RB66, RB67, RB68, RB69, RB70, RB71, RB72,
        RB73, RB74, RB75, RB76, RB77, RB78, RB79, RB80, RB81, RB82, RB83, RB84,
        RB85, RB86, RB87, RB88, RB89, RB90, RB91, RB92, RB93, RB94, RB95, RB96
    };

    /*PCR BARCODE*
     * All BC = All RB
     */
    inline std::unordered_map<std::string, SequenceInfo> get_trim_info()
    {
        std::unordered_map<std::string, SequenceInfo> trim_info;
        trim_info.try_emplace("SQK-LSK114",
                              "SQK-LSK114", LA_ADAPTER_5, LSK_TOP5END, LA_ADAPTER_3, LSK_TOP3END
        );
        trim_info.try_emplace("SQK-RAD114", "SQK-RAD114", RA_ADAPTER, RAD_TOP5END);
        trim_info.try_emplace("SQK-ULK114", "SQK-ULK114", RA_ADAPTER, RAD_TOP5END);
        trim_info.try_emplace("SQK-PCS114",
                              "SQK-PCS114",
                              SSPII,
                              PCS_TOP5END,
                              CRTA,
                              PCS_TOP3END,
                              CRTA_REV_COM,
                              PCS_BOT5END,
                              SSPII_REV_COM,
                              PCS_BOT3END
        );
        for (int i{0}; i < 96; i++) {
            if (i < 24) {
                std::string nbd24_name{fmt::format("SQK-NBD114.24-{}", i + 1)};
                std::string rbk24_name{fmt::format("SQK-RBK114.24-{}", i + 1)};
                std::string pcb24_name{fmt::format("SQK-PCB114.24-{}", i + 1)};
                trim_info.try_emplace(nbd24_name,
                                      nbd24_name,
                                      NA_ADAPTER_5 + NB_LEFT_FLANKING_5 + myutility::rev_com(NB_VEC[i]) +
                                      NB_RIGHT_FLANKING_5,
                                      NBD_TOP5END,
                                      NB_LEFT_FLANKING_3 + NB_VEC[i] + NB_RIGHT_FLANKING_3 + NA_ADAPTER_3,
                                      NBD_TOP3END
                );
                trim_info.try_emplace(rbk24_name,
                                      rbk24_name,
                                      RB_LEFT_FLANKING + RB_VEC[i] + RB_RIGTH_FLANKING,
                                      RBK_TOP5END
                );
                trim_info.try_emplace(pcb24_name,
                                      pcb24_name,
                                      RB_VEC[i] + SSPII, PCB_TOP5END, CRTA + myutility::rev_com(RB_VEC[i]),
                                      PCB_TOP3END,
                                      RB_VEC[i] + CRTA_REV_COM, PCB_BOT5END,
                                      SSPII_REV_COM + myutility::rev_com(RB_VEC[i]),
                                      PCB_BOT3END
                );
            }
            std::string nbd96_name{fmt::format("SQK-NBD114.96-{}", i + 1)};
            std::string rbk96_name{fmt::format("SQK-RBK114.96-{}", i + 1)};
            trim_info.try_emplace(nbd96_name,
                                  nbd96_name,
                                  NA_ADAPTER_5 + NB_LEFT_FLANKING_5 + myutility::rev_com(NB_VEC[i]) +
                                  NB_RIGHT_FLANKING_5,
                                  NBD_TOP5END,
                                  NB_LEFT_FLANKING_3 + NB_VEC[i] + NB_RIGHT_FLANKING_3 +
                                  NA_ADAPTER_3,
                                  NBD_TOP3END
            );
            trim_info.try_emplace(rbk96_name,
                                  rbk96_name,
                                  RB_LEFT_FLANKING + RB_VEC[i] + RB_RIGTH_FLANKING, RBK_TOP5END
            );
        }
        return trim_info;
    }
}

#endif //NANOFQ_ADAPTER_H
