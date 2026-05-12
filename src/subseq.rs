use crate::input_type::{InputType, check_input_type};
use crate::utils::{quit_with_error, rev_com};
use bio::bio_types::strand::ReqStrand;
use clap::ArgMatches;
use needletail::parser::{LineEnding, write_fastq};
use needletail::{FastxReader, Sequence, parse_fastx_file};
use regex::Regex;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{FetchDefinition, IndexedReader, Read, Record as BamRecord};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

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
                            let fastq_record = brecord2frecord(&record, name);
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
fn brecord2frecord(record: &BamRecord, read_name: &str) -> String {
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
        let fastq_record = brecord2frecord(&record, name);
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
    writer: &mut W,
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
                                    writeln!(writer, "{}", fastq_record).unwrap();
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
                                Some(fastq_record) => writeln!(writer, "{fastq_record}").unwrap(),
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
pub fn reads_from_fastq<W: Write>(
    read_names: ReadNames,
    reader: &mut Box<dyn FastxReader>,
    writer: &mut W,
) {
    let mut read_names = read_names.get_read_names();
    let mut read_idx = 1;
    while let Some(Ok(record)) = reader.next() {
        if read_names.is_empty() {
            return;
        }
        let mut headers = record.id().splitn(2, |x| x.is_ascii_whitespace());
        let name = headers
            .next()
            .expect(&format!("Parse read name failed at {read_idx}th record"));
        let name =
            str::from_utf8(name).expect(&format!("Parse read name failed at {read_idx}th record"));
        if read_names.remove(name) {
            write_fastq(
                record.id(),
                record.sequence(),
                record.qual(),
                writer,
                LineEnding::Unix,
            )
            .expect(&format!(
                "Failed write {read_idx}th fastq record into output"
            ))
        }
        read_idx += 1;
    }
    if !read_names.is_empty() {
        for not_found_name in read_names {
            eprintln!("couldn't found {not_found_name}")
        }
    }
}

pub fn run_subseq(subseq_cmd: &ArgMatches) {
    let input = subseq_cmd.get_one::<String>("input").unwrap();
    let output = subseq_cmd.get_one::<String>("output");
    let names = subseq_cmd.get_one::<String>("names");
    let names_file = subseq_cmd.get_one::<String>("names_file");
    let region = subseq_cmd.get_one::<String>("region");
    let bed = subseq_cmd.get_one::<String>("bed");
    let input_t = check_input_type(input);
    let names_sum = names.is_some() as u8 + names_file.is_some() as u8;
    let names_region_sum = names_sum + region.is_some() as u8 + bed.is_some() as u8;
    // let mut writer = BufWriter::new(File::create(output.unwrap()).expect(
    //     &format!("Failed to create {}", output.unwrap())
    // ));
    match input_t {
        InputType::OneFastqGzippedFile | InputType::OneFastqFile => {
            if names_sum != 1 {
                quit_with_error(
                    "For fastq[.gz] input, JUST one of [\"names\", \"names_file\"] must be specified",
                )
            }

            let read_names = if names.is_some() {
                ReadNames::FromCli(names.unwrap())
            } else {
                ReadNames::FromFile(names_file.unwrap())
            };
            let mut fastq_reader =
                parse_fastx_file(input).expect(&format!("Failed to read {input}"));
            match output {
                None => reads_from_fastq(read_names, &mut fastq_reader, &mut io::stdout()),
                Some(output) => {
                    let mut writer = BufWriter::new(
                        File::create(output).expect(&format!("Failed to create {output}")),
                    );
                    reads_from_fastq(read_names, &mut fastq_reader, &mut writer)
                }
            }
        }
        InputType::IndexedBam => {
            if names_region_sum != 1 {
                quit_with_error(
                    "For indexed bam, JUST one of [\"names\", \"names_file\", \"region\", \"bed\"] must be specified",
                )
            }
            let reads_in_bam = if names.is_some() {
                ReadsInBam::from_names_string(&names.unwrap())
            } else if names_file.is_some() {
                ReadsInBam::from_names_file(&names_file.unwrap())
            } else if region.is_some() {
                ReadsInBam::from_region_string(&region.unwrap())
            } else {
                ReadsInBam::from_bed(&bed.unwrap())
            };
            let mut bam_reader1 = IndexedReader::from_path(input).unwrap();
            bam_reader1.set_threads(4).unwrap();
            let mut bam_reader2 = IndexedReader::from_path(input).unwrap();
            bam_reader2.set_threads(2).unwrap();
            match output {
                None => reads_from_bam(
                    reads_in_bam,
                    &mut bam_reader1,
                    &mut bam_reader2,
                    &mut io::stdout(),
                ),
                Some(output) => {
                    let mut writer =
                        File::create(output).expect(&format!("Failed to create {output}"));
                    reads_from_bam(
                        reads_in_bam,
                        &mut bam_reader1,
                        &mut bam_reader2,
                        &mut writer,
                    );
                }
            }
        }
        _ => quit_with_error("Only fastq[.gz] or indexed bam file supported, check your input"),
    }
}
