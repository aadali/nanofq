#ifndef NANOFQ_ADAPTER_H
#define NANOFQ_ADAPTER_H

#include <string>
// https://nanoporetech.com/document/chemistry-technical-document#adapter-sequences
namespace barcode_info {
    /*
    * SQK-LSK114
    * LSK114 library reads structure
    *           |--->  LA_ADAPTER_5   <----| | insert Seq | |--->   LA_ADAPTER_3   <---|
    * 5-TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT-..............-AGCAATACGTAACTGAACGAAGTACAGG-3
    *
    * 3' end always is truncated
    */
    const std::string LA_ADAPTER_5 = "CCTGTACTTCGTTCAGTTACGTATTGCT";
    const std::string LA_ADAPTER_3 = "AGCAATACGTAACTGAACGAAGTACAGG"; // the first 15 or even less bases is enough, because 3'end always is truncated


    /*
     * SQK-NBD114-24; SQK-NBD114-96
     * NBD114-24/96 library reads structure
     * Example for Native Barcode01
     *           |NA_ADAPTER_5                |L_F_5   |Native Barcode01        |R_F_5   |insert Seq         |L_F_3   |Barcode01_rev_com       |R_F_3         |NA_ADAPTER_3
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
     * 3-...CCCAAABBBBAABBBBAABBBBAABBBBAAAGCAATATCAGCACCAACAGAAA .........TTTTTTTTTTTTTTT CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG...-5
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

    /*PCR BARCODE*
     * All BC = All RB
     */
}

#endif //NANOFQ_ADAPTER_H
