use crate::fastq::{FastqReader, NanoRead};
use crate::utils::{quit_with_error, rev_com};
use bio::bio_types::strand::ReqStrand;
use regex::Regex;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record as BamRecord};
use seq_io::fastq::Record;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

const FETCH_POS_PADDING: i64 = 50;
pub enum ReadNames<'a> {
    FromCli(&'a str),
    FromFile(&'a str),
}

pub enum ReadsInBam<'a> {
    Regions(Vec<(String, u32, u32)>),
    Names(ReadNames<'a>),
}

impl<'a> ReadsInBam<'a> {
    pub fn from_names_string(names: &'a str) -> Self {
        ReadsInBam::Names(ReadNames::FromCli(names))
    }

    pub fn from_names_file(names_file: &'a str) -> Self {
        ReadsInBam::Names(ReadNames::FromFile(names_file))
    }

    pub fn from_region_string(regions: &'a str) -> Self {
        let mut fetch_regions = vec![];
        let pat = Regex::new("[\\w-]+:[1-9]+[0-9]*-[1-9]+[0-9]*").unwrap();
        for region in regions.split(',') {
            if pat.is_match(region) {
                let x = region.split(':').collect::<Vec<&str>>();
                let contig = x[0];
                let y = x[1].split('-').collect::<Vec<&str>>();
                let start = y[0].parse::<u32>().unwrap_or_else(|x| {
                    quit_with_error(&format!(
                        "{x}, Bad region: {region}, pos should be less than 4294967296"
                    ))
                });
                let end = y[1].parse::<u32>().unwrap_or_else(|x| {
                    quit_with_error(&format!(
                        "{x}, Bad region: {region}, pos should be less than 4294967296"
                    ))
                });
                fetch_regions.push((contig.to_string(), start, end));
            } else {
                quit_with_error(&format!("Bad region format: {region}"))
            }
        }
        ReadsInBam::Regions(fetch_regions)
    }

    pub fn from_bed(bed: &str) -> Self {
        let mut fetch_region = vec![];
        let x = std::fs::read_to_string(bed).expect("Failed to read {bed}");
        for region in x.lines() {
            let line = region.split('\t').take(3).collect::<Vec<&str>>();
            match &line[..] {
                [contig, start, end] => {
                    let start = start.parse::<u32>().unwrap_or_else(|x| {
                        quit_with_error(&format!(
                            "{x}, Bad region: {region}, pos should be less than 4294967296"
                        ))
                    });
                    let end = end.parse::<u32>().unwrap_or_else(|x| {
                        quit_with_error(&format!(
                            "{x}, Bad region: {region}, pos should be less than 4294967296"
                        ))
                    });
                    fetch_region.push((contig.to_string(), start, end));
                }
                _ => quit_with_error("Bad bed file, 3 columns(Contig, Start, End) needed"),
            }
        }
        ReadsInBam::Regions(fetch_region)
    }
}

impl<'a> ReadNames<'a> {
    fn get_read_names(self) -> HashSet<String> {
        let mut read_names_set = HashSet::new();
        match &self {
            ReadNames::FromCli(names) => {
                for read_name in names.split(',') {
                    if read_name.len() == 0 {
                        continue;
                    }
                    for b in read_name.as_bytes() {
                        if !b.is_ascii() {
                            quit_with_error(&format!("Illegal byte found in name {}", read_name))
                        }
                        if b.is_ascii_whitespace() {
                            quit_with_error(&format!("White space found in name {}", read_name))
                        }
                    }
                    read_names_set.insert(read_name.to_string());
                }
            }
            ReadNames::FromFile(names_file) => {
                let f = File::open(&names_file).expect(&format!("Failed when open {names_file}"));
                let names_reader = BufReader::new(f);
                for (line_number, read_name_line_res) in names_reader.lines().enumerate() {
                    let read_name_line = read_name_line_res.expect(&format!(
                        "Get read name from line {} failed",
                        line_number + 1
                    ));
                    if read_name_line.len() == 0 {
                        continue;
                    }
                    for b in read_name_line.as_bytes() {
                        if !b.is_ascii() {
                            quit_with_error(&format!(
                                "Illegal byte found in line {}",
                                line_number + 1
                            ))
                        }
                        if b.is_ascii_whitespace() {
                            quit_with_error(&format!(
                                "White space found in line {}",
                                line_number + 1
                            ))
                        }
                    }
                    read_names_set.insert(read_name_line);
                }
            }
        }
        read_names_set
    }
}

///
/// Determines if a given BAM record does not contain hard clipping.
///
/// Hard clipping is indicated by the 'H' operation in the CIGAR string of a BAM record.
/// This function checks the first and last elements of the CIGAR string to determine
/// if there is any hard clipping present. If either the first or last CIGAR operation
/// is a hard clip, the function returns `false`, indicating that the record contains
/// hard clipping. Otherwise, it returns `true`.
///
fn record_has_no_hard(record: &BamRecord) -> bool {
    let record_cigar = record.cigar();
    let first_cigar = record_cigar.first();
    let last_cigar = record_cigar.last();
    match [first_cigar, last_cigar] {
        [Some(&Cigar::HardClip(_)), _] => false,
        [_, Some(&Cigar::HardClip(_))] => false,
        _ => true,
    }
}

type SaTag<'a> = (&'a str, i64, ReqStrand);

///
/// Return Supplementary Alignment(SA) info from SA tag, if a bam record has sa aux,
/// otherwise, return None
///
fn get_sa_tag_from_record(record: &'_ BamRecord) -> Option<Vec<SaTag<'_>>> {
    let sa_aux_res = record.aux(b"SA");
    match sa_aux_res {
        Ok(Aux::String(sa)) => {
            let mut sa_infos = vec![];
            for each_sa in sa[..sa.len() - 1].split(';') {
                let x = each_sa.split(',').take(3).collect::<Vec<&str>>();
                match &x[..] {
                    &[a, b, c] => sa_infos.push((
                        a,
                        b.parse::<i64>().unwrap(),
                        if c == "+" {
                            ReqStrand::Forward
                        } else {
                            ReqStrand::Reverse
                        },
                    )),
                    _ => quit_with_error(&format!("Bad SA tag format: {{{sa}}}")),
                }
            }
            Some(sa_infos)
        }
        _ => None,
    }
}

///
/// Retrieves a FASTQ record from BAM file using Supplementary Alignment (SA) tags.
///
/// # Arguments
/// * `bam_reader` - A mutable reference to an `IndexedReader` for reading BAM records.
/// * `sa_tags_opt` - An optional vector of `SaTag` structs. Each `SaTag` contains:
///     - The contig name as a string slice.
///     - The start position (1-based).
///     - The strand information.
/// * `read_name` - The name of the read to match, as a string slice.
///
/// # Returns
/// * `Option<String>` - The FASTQ record as a `String` if found, or `None` if not found.
///
/// # Behavior
/// - Iterates over each `SaTag` in `sa_tags_opt` to find matching reads.
/// - Uses the `fetch` method on `bam_reader` to limit the search within a specified range around the SA tag's start position.
/// - Checks each BAM record against the provided `read_name`, start position, and strand.
/// - Converts the first matching BAM record without hard clipping into a FASTQ format and returns it.
/// - If no matching record is found, or if `sa_tags_opt` is `None`, returns `None`.
///
/// # Notes
/// - The function assumes that positions in SA tags are 1-based, while BAM record positions are 0-based.
/// - Only primary alignments are considered; secondary alignments are skipped.
/// - The function checks for hard clipping in the BAM record before converting to FASTQ.
///
fn get_fastq_record_from_sa_tags(
    bam_reader: &mut IndexedReader,
    sa_tags_opt: Option<Vec<SaTag>>,
    read_name: &str,
) -> Option<String> {
    match sa_tags_opt {
        None => None,
        Some(sa_tags) => {
            for sa_tag in &sa_tags {
                let contig = &sa_tag.0;
                let start = sa_tag.1;
                let strand = sa_tag.2;
                let _ = bam_reader.fetch((
                    contig,
                    start - FETCH_POS_PADDING,
                    start + FETCH_POS_PADDING,
                ));
                let mut record = BamRecord::new();
                while let Some(_) = bam_reader.read(&mut record) {
                    if record.is_secondary() {
                        continue;
                    }
                    let name = std::str::from_utf8(record.qname()).unwrap();
                    let pos = record.pos();
                    // position returned by record.pos() is 0-based
                    // but position in SA tag is 1-based
                    if read_name == name && pos + 1 == start && strand == record.strand() {
                        if record_has_no_hard(&record) {
                            let fastq_record = from_bam_record_to_fastq_record(&record, name);
                            return Some(fastq_record);
                        }
                    }
                }
            }
            None
        }
    }
}

///
/// Converts a BAM record into a FASTQ format string.
/// This function should be used that bam record has full original seq length
///
/// # Arguments
///
/// * `record` - A reference to the BAM record to be converted.
/// * `read_name` - The name of the read to be used in the FASTQ header.
///
/// # Returns
///
/// A `String` that represents the FASTQ formatted data.
///
fn from_bam_record_to_fastq_record(record: &BamRecord, read_name: &str) -> String {
    let mut read_fastq = String::new();
    let seq = record.seq().as_bytes();
    let seq = std::str::from_utf8(&seq).unwrap();
    let qual = record.qual();
    read_fastq.push_str(&format!("@{read_name}\n"));
    if record.strand() == ReqStrand::Forward {
        let qual_str = qual.iter().map(|x| (*x + 33) as char).collect::<String>();
        read_fastq.push_str(&format!("{seq}\n+\n{}", qual_str));
    } else {
        let qual_str = qual
            .iter()
            .rev()
            .map(|x| (*x + 33) as char)
            .collect::<String>();
        read_fastq.push_str(&format!("{}\n+\n{}", rev_com(seq), qual_str,));
    }
    read_fastq
}

///
/// Attempts to find and return a FASTQ record for a given BAM record in indexed bam reader.
///
/// # Arguments
/// * `record` - A reference to the BAM record from which to extract the FASTQ information.
/// * `bam_reader` - A mutable reference to an `IndexedReader` used to read additional BAM records if necessary.
/// * `name` - The name of the sequence, used in constructing the FASTQ record.
///
/// # Returns
/// * `Option<String>` - The FASTQ record as a `String` if successfully found, or `None` if not applicable or not found.
///
/// # Details
/// - If the `record` is unmapped, primary, or has no hard clipping in its CIGAR string, the function constructs the FASTQ record directly from the BAM record.
/// - For supplementary alignments (typically indicated by the presence of hard clipping in the CIGAR), the function attempts to retrieve the full sequence and quality scores using SA tags from the BAM record, requiring access to the `bam_reader` to fetch related records.
/// - The function checks specific flags and conditions on the BAM record to determine the appropriate method for extracting or reconstructing the FASTQ record.
/// - The returned FASTQ record, when present, is formatted according to the standard FASTQ format, including the sequence name, sequence, and quality scores.
fn find_fastq_record(
    record: &BamRecord,
    bam_reader: &mut IndexedReader,
    name: &str,
) -> Option<String> {
    if record.is_unmapped() || record.flags() & 0x900 == 0 || record_has_no_hard(&record) {
        // when record is unmapped or primary or alignment that is not secondary and no hard clipping in cigar
        // full seq and quality will be found
        let fastq_record = from_bam_record_to_fastq_record(&record, name);
        Some(fastq_record)
    } else {
        // for supplementary alignment,
        // typically, one hard clipping will be found at first or last cigar
        let sa_tags = get_sa_tag_from_record(&record);
        let fastq_record_opt = get_fastq_record_from_sa_tags(bam_reader, sa_tags, name);
        fastq_record_opt
    }
}
///
/// Reads and processes records from a BAM file, writing the resulting FASTQ records to a provided writer.
///
/// # Arguments
///
/// * `reads_in_bam` - A struct that specifies either a genomic region (contig, start, end) or a set of read names to fetch from the BAM file.
/// * `bam_reader` - A mutable reference to an `IndexedReader` for reading from the primary BAM file.
/// * `bam_reader2` - A mutable reference to a second `IndexedReader`, used for cross-referencing or fetching additional data.
/// * `w` - A mutable reference to a type implementing the `Write` trait, where the output FASTQ records will be written.
///
/// # Behavior
///
/// - If `reads_in_bam` is of variant `Region`, it fetches all records within the specified genomic region, skipping secondary alignments and duplicates. For each unique read name, it attempts to find a corresponding FASTQ record using `find_fastq_record` with `bam_reader2`. Found FASTQ records are written to `w`.
/// - If `reads_in_bam` is of variant `Names`, it iterates over all records in the BAM file, stopping once all specified read names have been processed. It skips secondary alignments and, for each matching read name, tries to find and write the FASTQ record, removing the name from the list upon processing.
/// - Errors during BAM parsing lead to an immediate termination with an error message.
/// - The function uses a `HashSet` to track found read names, ensuring no duplicates are processed.
///
pub fn reads_from_bam<W: Write>(
    reads_in_bam: ReadsInBam,
    bam_reader: &mut IndexedReader,
    bam_reader2: &mut IndexedReader,
    w: &mut W,
) {
    match reads_in_bam {
        ReadsInBam::Regions(regions) => {
            let header = bam_reader.header().to_owned();
            for region in regions {
                let (contig, start, end) = region;
                if header.tid(contig.as_bytes()).is_none() {
                    quit_with_error(&format!("Could found contig: {}", &contig))
                }
                let _ = bam_reader
                    .fetch((contig.as_bytes(), start, end))
                    .expect(&format!("Fetched region ({contig}, {start}, {end}) failed"));
                let mut record = BamRecord::new();
                let mut found_read_names = HashSet::new();
                while let Some(record_res) = bam_reader.read(&mut record) {
                    match record_res {
                        Ok(_) => {
                            if record.is_secondary() {
                                continue;
                            }
                            let name = std::str::from_utf8(record.qname()).unwrap();
                            if found_read_names.contains(name) {
                                continue;
                            }
                            let fastq_record_opt = find_fastq_record(&record, bam_reader2, name);
                            match fastq_record_opt {
                                Some(fastq_record) => {
                                    writeln!(w, "{}", fastq_record).unwrap();
                                    found_read_names.insert(name.to_string());
                                }
                                None => {
                                    eprintln!(
                                        "Fastq Record with name: {name}, full length seq couldn't be found"
                                    )
                                }
                            }
                        }
                        Err(_) => quit_with_error("Failure parsing bam file"),
                    }
                }
            }
        }

        ReadsInBam::Names(names) => {
            let mut read_names = names.get_read_names();
            bam_reader.fetch(FetchDefinition::All).unwrap();
            let mut record = BamRecord::new();
            while let Some(record_res) = bam_reader.read(&mut record) {
                if read_names.is_empty() {
                    return;
                }
                match record_res {
                    Ok(_) => {
                        if record.is_secondary() {
                            continue;
                        }
                        let name = std::str::from_utf8(record.qname()).unwrap();
                        if read_names.contains(name) {
                            // everytime code run to here, one read name will be absolutely removed from read_names: HashSet<String>,
                            // even the fastq record not found in this bam
                            let fastq_record_opt = find_fastq_record(&record, bam_reader2, name);
                            match fastq_record_opt {
                                Some(fastq_record) => writeln!(w, "{fastq_record}").unwrap(),
                                None => {
                                    eprintln!(
                                        "Fastq Record with name: {name}, full length seq couldn't be found"
                                    )
                                }
                            }
                            read_names.remove(name);
                        }
                    }
                    Err(_) => quit_with_error("Failure parsing bam file"),
                }
            }
            if !read_names.is_empty() {
                for name in &read_names {
                    eprintln!("Fastq Record with name: {name}, full length seq couldn't be found")
                }
            }
        }
    }
}

///
/// Reads FASTQ records from a reader and writes them to a writer if their names are in the provided list.
///
/// # Arguments
///
/// * `read_names` - A `ReadNames` struct containing the list of read names to be processed.
/// * `reader` - A mutable reference to a `FastqReader<R>` where `R: Read`, used for reading FASTQ records.
/// * `writer` - A mutable reference to a `W: Write` type, used for writing selected FASTQ records.
///
/// # Behavior
///
/// - Iterates over FASTQ records from the given reader.
/// - Checks if the record's name is in the `read_names` list.
/// - If a match is found, the record is written to the writer and removed from the `read_names` list.
/// - The function stops processing once all specified read names have been found or the end of the reader is reached.
/// - After processing, if there are any read names left in the `read_names` list that were not found, it prints an error message for each missing name.
///
pub fn reads_from_fastq<R: io::Read, W: Write>(
    read_names: ReadNames,
    reader: &mut FastqReader<R>,
    writer: &mut W,
) {
    let mut read_names = read_names.get_read_names();
    while let Some(fastq_record_res) = reader.next() {
        if read_names.is_empty() {
            return;
        }
        match fastq_record_res {
            Ok(fastq_record) => {
                let id = fastq_record
                    .head()
                    .split(|x: &u8| x.is_ascii_whitespace())
                    .next()
                    .expect("Get read id failed");
                let name = std::str::from_utf8(id)
                    .unwrap_or_else(|err| quit_with_error(&format!("{err}")));
                if read_names.remove(name) {
                    NanoRead::write(&fastq_record, writer).unwrap()
                }
            }
            Err(err) => quit_with_error(&format!("{err}")),
        }
    }
    if !read_names.is_empty() {
        for not_found_name in read_names {
            eprintln!("couldn't found {not_found_name}")
        }
    }
}

// pub fn test_from_bam() {
//     // let bam_file = "/home/a/test_data/hac__no_trim__emit_moves.bam";
//     let bam_file = "/home/a/test_data/hac__5mC__no_trim__emit_moves.call.sorted.bam";
//     let read_names_str = "f4d4f0ea-8cd7-4563-b2e5-c703a1ce7636,xyasae,bc46205a-066e-4fe9-b9a0-d604d83e7fe9,f85bd896-cb7d-4f9e-9e4f-15f5b1f09184,1550195c-25a3-45e9-9b0c-52c0cc83e3b5,,bfe06368-8fb9-4ef6-a86f-d4b985346b11,cd1c81a2-d835-413e-9671-1f90053bfa20,6f5f8f84-428c-489f-87b7-b9855855bccc,579fa65b-7812-4ad4-b0e8-2320c6b3220a,2f0bca36-f577-4f21-a37f-899d84373f06,daba5f75-c2ef-4cd7-a23d-8e5bebde5642,af0b6021-b939-4a3a-a5dd-cf7d01c2d3da,2c94aca0-d843-4392-9f64-567fc746f0f6";
//     let read_names_file = "/home/a/test_data/read_names.txt";
//     let reads_in_bam_name = ReadsInBam::Names(ReadNames::FromCli(read_names_str));
//     let reads_in_bam_region = ReadsInBam::Regions(vec![("chr21".to_string(), 28336100, 28336900)]);
//     let mut bam_reader1 = IndexedReader::from_path(bam_file).unwrap();
//     let _ = bam_reader1.fetch(("chr99, 9023, 120324"));
//     println!("hello world, start");
//     for x in bam_reader1.records() {
//         let x = x.unwrap();
//         println!("name: {} => {}", x.pos(), x.tid())
//
//     }
//     return;
//     // let _ = bam_reader1.fetch(("chr21", 23423123,28901129)).unwrap();
//     // let _ = bam_reader1.fetch(("chr21", 24322123,25001129)).unwrap();
//     // let _ = bam_reader1.fetch(("chr21", 24324119,25001129)).unwrap();
//     // 60d29c51-85f6-433b-8b4f-5db0a0bb8b9f    20      24278753        24324119
//     // b7651a60-260a-4a1d-967e-6111d85d81d5    20      24282908        24329856
//     let mut x = vec![];
//     for record in bam_reader1.records() {
//         let record = record.unwrap();
//         let name = std::str::from_utf8(record.qname()).unwrap();
//         x.push(format!(
//             "{}\t{}\t{}\t{}",
//             std::str::from_utf8(record.qname()).unwrap(),
//             record.tid(),
//             record.pos(),
//             record.reference_end()
//         ));
//     }
//     let mut bam_reader2 = IndexedReader::from_path(bam_file).unwrap();
//     let outfile1 = "/home/a/test_data/bam1.fastq";
//     let mut writer1 = std::fs::File::create(outfile1).unwrap();
//     reads_from_bam(
//         // reads_in_bam_name,
//         reads_in_bam_region,
//         &mut bam_reader1,
//         &mut bam_reader2,
//         &mut writer1,
//     );
// }

#[cfg(test)]
mod test {
    use regex::Regex;
    #[test]
    fn t() {
        let pat = Regex::new("[\\w-]+:[1-9]+[0-9]*-[1-9]+[0-9]*").unwrap();
        let r = pat.is_match("chr21-a:990-12340");
        assert!(r)
    }

    #[test]
    fn t2() {
        let x = "99999999999";
        let y = x.parse::<u32>().unwrap();
    }
}
