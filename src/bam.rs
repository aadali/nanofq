use crate::fastq::{DORADO_TRIM_LEADING_BASE_NUMBER, EachStats};
use crate::utils::{get_q2p_table, quit_with_error};
use rayon::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::index;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{self, FetchDefinition, HeaderView, IndexedReader, Read};
use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::ops::AddAssign;
use std::path::Path;
use std::sync::OnceLock;

static ENCODED_BASES_GC_COUNT: OnceLock<HashMap<u8, usize>> = OnceLock::new();
pub fn get_encoded_bases_gc_count_table() -> &'static HashMap<u8, usize> {
    ENCODED_BASES_GC_COUNT.get_or_init(|| {
        HashMap::from([
            (18, 1),  // AC: 18
            (20, 0),  // AG: 20
            (24, 0),  // AT: 24
            (17, 0),  // AA: 17
            (129, 0), // TA: 129
            (136, 0), // TT: 136
            (132, 1), // TG: 132
            (130, 1), // TC: 130
            (33, 1),  // CA: 33
            (40, 1),  // CT: 40
            (36, 2),  // CG: 36
            (34, 2),  // CC: 34
            (65, 1),  // GA: 65
            (66, 2),  // GC: 66
            (68, 2),  // GG: 68
            (72, 1),  // GT: 72
            (16, 0),  // A : 16
            (128, 0), // T : 128
            (32, 1),  // C : 32
            (64, 1),  // G : 64
        ])
    })
}
pub trait BamRecordStats {
    fn gc_count(&self) -> f32;
    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> Option<f64>;
    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> EachStats;
}
impl BamRecordStats for rust_htslib::bam::Record {
    fn gc_count(&self) -> f32 {
        let seq_len = self.qual().len();
        let gc_number: usize = self
            .seq()
            .encoded
            .iter()
            .map(|x| *get_encoded_bases_gc_count_table().get(x).unwrap_or(&0usize))
            .sum();
        gc_number as f32 / seq_len as f32
    }

    fn calculate_read_quality(&self, dont_use_dorado_quality: bool) -> Option<f64> {
        let quals = self.qual();
        if quals.len() == 0 {
            eprintln!(
                "Empty quality found for: {}",
                str::from_utf8(self.qname()).unwrap()
            );
            return None;
        }
        let real_seq_len = quals.len();
        if dont_use_dorado_quality {
            let avg_err_prob = quals
                .iter()
                .map(|x| get_q2p_table()[*x as usize + 33])
                .sum::<f64>()
                / real_seq_len as f64;
            let read_quality = avg_err_prob.log10() * -10.0;
            Some(read_quality)
        } else {
            let quality_tag_res = self.aux(b"qs");
            if quality_tag_res.is_ok() {
                let quality_tag = quality_tag_res.unwrap();
                let read_quality = match quality_tag {
                    Aux::Float(quality) => quality as f64,
                    Aux::Double(quality) => quality,
                    _ => {
                        quit_with_error(&format!(
                            "Parse qs tag for {:?} in tid: {}, at position: {}",
                            str::from_utf8(self.qname()).unwrap(),
                            self.tid(),
                            self.reference_start(),
                        ));
                    }
                };
                Some(read_quality)
            } else {
                let (seq_len, skip) = if real_seq_len > DORADO_TRIM_LEADING_BASE_NUMBER {
                    (
                        real_seq_len - DORADO_TRIM_LEADING_BASE_NUMBER,
                        DORADO_TRIM_LEADING_BASE_NUMBER,
                    )
                } else {
                    (real_seq_len, 0)
                };
                let avg_err_prob = quals
                    .iter()
                    .skip(skip)
                    .map(|x| get_q2p_table()[*x as usize + 33])
                    .sum::<f64>()
                    / seq_len as f64;
                Some(avg_err_prob.log10() * -10.0)
            }
        }
    }

    fn stats(&self, gc: bool, dont_use_dorado_quality: bool) -> EachStats {
        let len = self.qual().len();
        let read_quality = self.calculate_read_quality(dont_use_dorado_quality);
        let gc = if gc { Some(self.gc_count()) } else { None };
        if read_quality.is_none() {
            quit_with_error(&format!(
                "Empty quality found for: {}",
                str::from_utf8(self.qname()).unwrap()
            ))
        }
        (
            Box::new(str::from_utf8(self.qname()).unwrap().to_string()),
            len as u32,
            read_quality.unwrap() as f32,
            gc,
        )
    }
}

#[derive(Default, Debug)]
pub struct BasicBamStatistics {
    // only for Aligned Bam
    reads_mapped: usize,            // primary alignments reads number
    reads_unmapped: usize,          // unmapped reads number
    reads_mq0: usize,               // primary alignment reads with map quality == 0
    primary_alignment: usize,       // primary alignments
    supplementary_alignment: usize, // supplementary alignments
    secondary_alignment: usize,     // secondary alignments
    bases_mapped: usize,            // primary_alignment seq_len()
    bases_mapped_cigar: usize, //  Match + Ins of (primary_alignment and supplementary_alignment), more accurate
    mismatches: usize,         // NM tag
    error_rate: f64,           // mismatches / bases_mapped_cigar
    map_rate: f64,             // reads_mapped / (reads_mapped + reads_unmapped)
}
impl BasicBamStatistics {
    pub fn is_empty(&self) -> bool {
        self.reads_mapped == 0
            && self.reads_unmapped == 0
            && self.reads_mq0 == 0
            && self.reads_mq0 == 0
            && self.primary_alignment == 0
            && self.supplementary_alignment == 0
            && self.secondary_alignment == 0
            && self.bases_mapped == 0
            && self.bases_mapped_cigar == 0
            && self.bases_mapped_cigar == 0
            && self.mismatches == 0
    }
}
impl Display for BasicBamStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "#BamStatsSummary\n\
            ReadsMapped\t{}\n\
            ReadsUnmapped\t{}\n\
            ReadsMapQ0\t{}\t# primary alignments with map quality==0\n\
            PrimaryAlignment\t{}\n\
            SupplementaryAlignment\t{}\n\
            SecondaryAlignment\t{}\n\
            BasesMapped\t{}\t# primary alignments seq_len, ignore clipping\n\
            BasesMappedCigar\t{}\t# Match + Ins of (primary alignments and supplementary alignments), more accurate\n\
            Mismatches\t{}\t# NM tag\n\
            ErrorRate\t{:.6}\t# Mismatches / BasesMappedCigar\n\
            MapRate\t{:.6}\t# ReadsMapped / (ReadsMapped + ReadsUnmapped)\n",
            self.reads_mapped,
            self.reads_unmapped,
            self.reads_mq0,
            self.primary_alignment,
            self.supplementary_alignment,
            self.secondary_alignment,
            self.bases_mapped,
            self.bases_mapped_cigar,
            self.mismatches,
            self.error_rate,
            self.map_rate
        )
    }
}

impl AddAssign for BasicBamStatistics {
    fn add_assign(&mut self, rhs: Self) {
        self.reads_mapped += rhs.reads_mapped;
        self.reads_unmapped += rhs.reads_unmapped;
        self.reads_mq0 += rhs.reads_mq0;
        self.primary_alignment += rhs.primary_alignment;
        self.supplementary_alignment += rhs.supplementary_alignment;
        self.secondary_alignment += rhs.secondary_alignment;
        self.bases_mapped += rhs.bases_mapped;
        self.bases_mapped_cigar += rhs.bases_mapped_cigar;
        self.mismatches += rhs.mismatches;
        self.error_rate = self.mismatches as f64 / self.bases_mapped_cigar as f64;
        self.map_rate = self.reads_mapped as f64 / (self.reads_mapped + self.reads_unmapped) as f64;
    }
}
fn get_nm_aux(record: &bam::Record) -> usize {
    let nm_tag = record.aux(b"NM").unwrap_or(Aux::I32(0));
    match nm_tag {
        // NM tag in bam file is I32, but NM tag in sam file is U8, U16,
        // it maybe a bug in rust-htslib
        Aux::I8(nm) => nm as usize,
        Aux::U8(nm) => nm as usize,
        Aux::I16(nm) => nm as usize,
        Aux::U16(nm) => nm as usize,
        Aux::I32(nm) => nm as usize,
        Aux::U32(nm) => nm as usize,
        _ => {
            quit_with_error(&format!(
                "Parse NM tag failed for {}",
                str::from_utf8(record.qname()).unwrap()
            ));
        }
    }
}

fn get_cigar_base_length(cigar: &Cigar) -> usize {
    match cigar {
        Cigar::Match(m) | Cigar::Equal(m) | Cigar::Diff(m)=> *m as usize,
        Cigar::Ins(i) => *i as usize,
        _ => 0usize,
    }
}
fn bins_contig(contig_len: u64, bins_number: usize) -> Vec<(u64, u64)> {
    let mut edges = vec![];
    let step = contig_len / bins_number as u64;
    let mut left = 0u64;
    let mut right = 0u64;
    for i in 0..bins_number {
        if i == bins_number - 1 {
            right = contig_len;
            edges.push((left, right));
        } else {
            right = left + step;
            edges.push((left, right));
            left = right
        }
    }
    edges
}
fn bins(header: &'_ HeaderView, thread: usize) -> Vec<FetchDefinition<'_>> {
    let mut contigs_edges = vec![];
    let contig_numbers = header.target_count();
    for tig in 0..contig_numbers {
        let contig_len_opt = header.target_len(tig);
        if contig_len_opt.is_none() {
            quit_with_error(&format!("Get contig length for {} failed", tig))
        }
        let contig_len = contig_len_opt.unwrap();
        let edges = bins_contig(contig_len, thread);
        for edge in edges {
            contigs_edges.push(FetchDefinition::Region(
                tig as i32,
                edge.0 as i64,
                edge.1 as i64,
            ))
        }
    }
    contigs_edges.push(FetchDefinition::Unmapped);
    contigs_edges
}

#[derive(Debug, Eq, PartialEq)]
pub enum BamType {
    SAM,
    UnalignedBam,
    UnsortedBam,
    SortedUnindexedBam,
    IndexedBam,
}

pub fn check_bam_type(bam_file: &str) -> BamType {
    if !(bam_file.ends_with("bam") || bam_file.ends_with(".sam")) {
        quit_with_error("input file should be bam or sam file ends with \".bam\" or \".sam\"")
    }
    if bam_file.ends_with(".sam") {
        BamType::SAM
    } else {
        let bam_reader = rust_htslib::bam::Reader::from_path(bam_file)
            .expect("Make bam::Reader failed for check_bam_type");
        let header_view = bam_reader.header();
        if 0 == header_view.target_count() {
            BamType::UnalignedBam
        } else {
            let header_hashmap = bam::Header::from_template(header_view).to_hashmap();
            let hd_opt = header_hashmap.get(&"HD".to_string());
            if hd_opt.is_none() {
                BamType::UnsortedBam
            } else {
                let hd_linear_map = hd_opt.unwrap().first().unwrap();
                let so_opt = hd_linear_map.get(&"SO".to_string());
                match so_opt {
                    None => BamType::UnsortedBam,
                    Some(so) => {
                        if so == "coordinate" {
                            let bam_file_index_path =
                                std::path::PathBuf::from(&format!("{}.bai", bam_file));
                            if bam_file_index_path.exists() {
                                let bam_file_meta_res = Path::new(bam_file).metadata();
                                let bam_file_index_meta_res = bam_file_index_path.metadata();
                                if bam_file_index_meta_res.is_ok() && bam_file_meta_res.is_ok() {
                                    if bam_file_meta_res.unwrap().created().unwrap()
                                        < bam_file_index_meta_res.unwrap().created().unwrap()
                                    {
                                        BamType::IndexedBam
                                    } else {
                                        BamType::SortedUnindexedBam
                                    }
                                } else {
                                    BamType::SortedUnindexedBam
                                }
                            } else {
                                BamType::SortedUnindexedBam
                            }
                        } else {
                            BamType::UnsortedBam
                        }
                    }
                }
            }
        }
    }
}

fn stats_from_bam_reader<R>(
    bam_reader: &mut R,
    gc: bool,
    dont_use_dorado_quality: bool,
    region_start: i64,
) -> (BasicBamStatistics, Vec<EachStats>)
where
    R: bam::Read,
{
    let mut basic_bam_stats = BasicBamStatistics::default();
    let mut all_stats = vec![];
    let mut record = bam::Record::new();
    record.set_qname(b"InitBamRecord");
    while let Some(x) = bam_reader.read(&mut record) {
        if x.is_err() {
            quit_with_error(&format!(
                "Parse record failed: {}",
                str::from_utf8(record.qname()).unwrap()
            ));
        }

        if record.is_unmapped() {
            basic_bam_stats.reads_unmapped += 1;
            all_stats.push(record.stats(gc, dont_use_dorado_quality));
        } else {
            if record.pos() < region_start {
                continue;
            }
            if record.flags() & 0x900 == 0 {
                all_stats.push(record.stats(gc, dont_use_dorado_quality));
                basic_bam_stats.bases_mapped += record.seq_len();
                basic_bam_stats.primary_alignment += 1;
                basic_bam_stats.reads_mapped += 1;
                if record.mapq() == 0 {
                    basic_bam_stats.reads_mq0 += 1;
                }
                let nm = get_nm_aux(&record);
                basic_bam_stats.mismatches += nm;
                for cigar in record.cigar().iter() {
                    basic_bam_stats.bases_mapped_cigar += get_cigar_base_length(cigar);
                }
            } else {
                if record.is_supplementary() {
                    basic_bam_stats.supplementary_alignment += 1;
                    let nm = get_nm_aux(&record);
                    basic_bam_stats.mismatches += nm;
                    for cigar in record.cigar().iter() {
                        basic_bam_stats.bases_mapped_cigar += get_cigar_base_length(cigar);
                    }
                }
                if record.is_secondary() {
                    basic_bam_stats.secondary_alignment += 1;
                }
            }
        }
    }
    basic_bam_stats.error_rate =
        basic_bam_stats.mismatches as f64 / basic_bam_stats.bases_mapped_cigar as f64;
    basic_bam_stats.map_rate = basic_bam_stats.reads_mapped as f64
        / (basic_bam_stats.reads_mapped + basic_bam_stats.reads_unmapped) as f64;
    (basic_bam_stats, all_stats)
}

fn stats_indexed_bam_fetch(
    bam_reader: &mut IndexedReader,
    region: FetchDefinition,
    gc: bool,
    dont_use_dorado_quality: bool,
) -> (BasicBamStatistics, Vec<EachStats>) {
    let region_start = match &region {
        FetchDefinition::Region(_, region_start, _) => *region_start,
        FetchDefinition::Unmapped => i64::MIN,
        _ => quit_with_error("Bad FetchDefinition"),
    };
    let may_be_err_msg = format!("Fetch region: {:?} from IndexedReader failed", &region);
    let fetch_result = bam_reader.fetch(region);
    if fetch_result.is_err() {
        quit_with_error(&may_be_err_msg)
    }
    stats_from_bam_reader(bam_reader, gc, dont_use_dorado_quality, region_start)
}

pub fn index_bam(bam_file: &str, thread: usize) -> Result<(), anyhow::Error> {
    index::build(bam_file, None, index::Type::Bai, thread as u32)?;
    Ok(())
}

/// ```
/// for UnalignedBam/UnsortedBam/SAM file or all bam/sam from stdin
/// ```
pub fn stats_xam(
    bam_reader: &mut bam::Reader,
    thread: usize,
    gc: bool,
    dont_use_dorado_quality: bool,
) -> (BasicBamStatistics, Vec<EachStats>) {
    debug_assert!(thread > 0);
    bam_reader.set_threads(thread).unwrap();
    stats_from_bam_reader(bam_reader, gc, dont_use_dorado_quality, i64::MIN)
}

thread_local! {
    static INDEXED_BAM_READER: RefCell<IndexedReader> = panic!("!");
}

/// ```
/// For IndexedBam file. If it's a SortedUnindexBam, index it, turn it into IndexedBam
/// ```
pub fn stats_indexed_bam(
    bam_file: &str,
    thread: usize,
    gc: bool,
    dont_use_dorado_quality: bool,
) -> (BasicBamStatistics, Vec<EachStats>) {
    debug_assert_eq!(check_bam_type(bam_file), BamType::IndexedBam);
    let indexed_bam_reader =
        IndexedReader::from_path(bam_file).expect("Read indexed bam failed for stats_indexed_bam");
    let fetch_regions = bins(indexed_bam_reader.header(), thread);
    let files = (0..thread)
        .into_iter()
        .map(|_| bam_file.to_string())
        .collect::<Vec<String>>();

    rayon::ThreadPoolBuilder::new()
        .num_threads(thread)
        .start_handler(move |index| {
            let bam_file_string = files[index].clone();
            let mut reader = IndexedReader::from_path(&bam_file_string).expect("Open failed");
            reader.set_threads(thread).unwrap();
            INDEXED_BAM_READER.set(reader);
        })
        .thread_name(|x| format!("Thread: {x}"))
        .build_global()
        .expect("thread pool builder failed for stats_indexed_bam function");

    let result: Vec<_> = fetch_regions
        .into_par_iter()
        .map(|region| {
            INDEXED_BAM_READER.with_borrow_mut(|indexed_reader| {
                stats_indexed_bam_fetch(indexed_reader, region, gc, dont_use_dorado_quality)
            })
        })
        .collect();

    let mut basic_bam_stats = BasicBamStatistics::default();
    let mut all_stats = vec![];
    for x in result {
        basic_bam_stats += x.0;
        all_stats.extend(x.1);
    }
    (basic_bam_stats, all_stats)
}
