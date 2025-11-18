use crate::alignment::{LocalAligner, Scores};
use crate::fastq::{FastqReader, NanoRead};
use crate::trim::adapter::TrimConfig;
use crate::trim::trim_seq;
use crate::utils::rev_com;
use ansi_term::Color;
use bio::alphabets::dna;
use flate2::read::MultiGzDecoder;
use seq_io::fastq::{OwnedRecord, Record};
use std::fmt::Display;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

pub struct SubOwnedRecord {
    owned_record: OwnedRecord,
    start: usize,
    end: usize,
    rev_com: bool,
}

impl Record for SubOwnedRecord {
    fn head(&self) -> &[u8] {
        self.owned_record.head()
    }

    fn seq(&self) -> &[u8] {
        &self.owned_record.seq()[self.start..=self.end]
    }

    fn qual(&self) -> &[u8] {
        &self.owned_record.qual()[self.start..=self.end]
    }
}

impl SubOwnedRecord {
    pub fn write_rev_com(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
        write!(
            writer,
            "@{}\n{}\n+\n{}\n",
            format!(
                "{}_rev_com {}",
                self.owned_record.id()?,
                self.desc().unwrap_or(Ok("")).unwrap()
            ),
            rev_com(unsafe { std::str::from_utf8_unchecked(self.seq()) }),
            unsafe {
                String::from_utf8_unchecked(
                    self.qual().iter().rev().map(|x| *x).collect::<Vec<u8>>(),
                )
            }
        )?;
        Ok(())
    }

    pub fn write_rev_com_fa(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
        write!(
            writer,
            ">{}_rev_com\n{}\n",
            self.owned_record.id()?,
            rev_com(unsafe { std::str::from_utf8_unchecked(self.seq()) }),
        )?;
        Ok(())
    }

    pub fn write_fa(&self, writer: &mut dyn Write) -> Result<(), anyhow::Error> {
        write!(writer, ">{}\n{}\n", self.owned_record.id()?, unsafe {
            std::str::from_utf8_unchecked(self.seq())
        })?;
        Ok(())
    }
}

fn get_suitable_amplicon_by_align<R: Read>(
    fastq_reader: &mut FastqReader<R>,
    end5: &str,
    end3: &str,
    rev_com_end5: &str,
    rev_com_end3: &str,
    est_len: usize,
    range: u8,
) -> Vec<SubOwnedRecord> {
    let mut candidate_amplicon = vec![];
    let end_align_para = (180usize, 0.9, 0.9);
    let primer_cfg = TrimConfig {
        kit_name: "primer",
        end5: Some((end5, end_align_para)),
        end3: Some((end3, end_align_para)),
        rev_com_end5: Some((rev_com_end5, end_align_para)),
        rev_com_end3: Some((&rev_com_end3, end_align_para)),
    };
    let scores = Scores {
        match_: 3,
        mismatch: -3,
        gap_open: -5,
        gap_extend: -1,
    };
    let mut local_aligner = LocalAligner::new((200, 200), scores);
    let low = (est_len as f64 * ((100 - range) as f64 / 100.0)) as usize;
    let up = (est_len as f64 * ((100 + range) as f64 / 100.0)) as usize;
    loop {
        if let Some(ref_record_res) = fastq_reader.next() {
            match ref_record_res {
                Ok(ref_record) => {
                    if ref_record.seq().len() < low || ref_record.seq().len() > up {
                        continue;
                    }
                    let (trim_from, trim_to, _, is_rev_com) = trim_seq(
                        &primer_cfg,
                        ref_record.seq(),
                        &format!(
                            "{}: {}",
                            ref_record.id().expect("parse into read id error"),
                            ref_record.seq().len()
                        ),
                        &mut local_aligner,
                        false,
                        est_len / 2,
                        false,
                    );
                    if trim_from == 0 || trim_to == 0 || trim_to == ref_record.seq().len() {
                        continue;
                    }
                    candidate_amplicon.push(SubOwnedRecord {
                        owned_record: ref_record.to_owned_record(),
                        start: trim_from,
                        end: trim_to - 1,
                        rev_com: is_rev_com,
                    })
                }
                Err(err) => {
                    eprintln!("{:?}", err);
                    std::process::exit(1);
                }
            }
        } else {
            break;
        }
    }
    candidate_amplicon
}
pub fn get_candidate_amplicon<P: AsRef<Path>>(
    input_fastq: Option<P>,
    fwd_primer: &str,
    rev_primer: &str,
    est_len: usize,
    range: u8,
) -> Result<Vec<SubOwnedRecord>, anyhow::Error> {
    let end5 = fwd_primer;
    let revcom_rev_primer = dna::revcomp(rev_primer.as_bytes());
    let end3 = std::str::from_utf8(&revcom_rev_primer)
        .expect("Get reverse complementary sequence of rev primer failed");
    let rev_com_end5 = rev_primer;
    let revcom_fwd_primer = dna::revcomp(fwd_primer.as_bytes());
    let rev_com_end3 = std::str::from_utf8(&revcom_fwd_primer)
        .expect("Get reverse complementary sequence of fwd primer failed");
    let sub_owned_record_vec = if input_fastq.is_none() {
        let mut fastq_reader = FastqReader::new(std::io::stdin());
        get_suitable_amplicon_by_align(
            &mut fastq_reader,
            end5,
            end3,
            rev_com_end5,
            rev_com_end3,
            est_len,
            range,
        )
    } else {
        let input_fastq_path = PathBuf::from(input_fastq.unwrap().as_ref());
        debug_assert!(input_fastq_path.is_file());
        if input_fastq_path.to_str().unwrap().ends_with(".gz") {
            let mut fastq_reader =
                FastqReader::new(MultiGzDecoder::new(std::fs::File::open(&input_fastq_path)?));
            get_suitable_amplicon_by_align(
                &mut fastq_reader,
                end5,
                end3,
                rev_com_end5,
                rev_com_end3,
                est_len,
                range,
            )
        } else {
            let mut fastq_reader = FastqReader::new(std::fs::File::open(&input_fastq_path)?);
            get_suitable_amplicon_by_align(
                &mut fastq_reader,
                end5,
                end3,
                rev_com_end5,
                rev_com_end3,
                est_len,
                range,
            )
        }
    };
    if sub_owned_record_vec.len() < 3 {
        eprintln!("{}", Color::Red.paint("few fastq records got"));
        std::process::exit(1)
    }
    Ok(sub_owned_record_vec)
}

pub fn filter_candidate_amplicon(
    candidate_amplicon: Vec<SubOwnedRecord>,
    number: usize,
) -> Vec<SubOwnedRecord>
// Retain the data within half-standard deviation
{
    let mut final_amplicon = vec![];
    let total_len = candidate_amplicon
        .iter()
        .fold(0usize, |sum, each_sub_owned_record| {
            each_sub_owned_record.seq().len() + sum
        });
    let mean_len = total_len as f64 / candidate_amplicon.len() as f64;
    let std_len = (candidate_amplicon
        .iter()
        .fold(0.0f64, |sum, each_sub_owned_record| {
            sum + (each_sub_owned_record.seq().len() as f64 - mean_len).powf(2.0)
        })
        / (candidate_amplicon.len() as f64))
        .sqrt();
    let (len_lower, len_upper) = (mean_len - 0.5 * std_len, mean_len + 0.5 * std_len);
    for each_candidate_amplicon in candidate_amplicon {
        if (each_candidate_amplicon.qual().len() as f64) > len_lower
            && (each_candidate_amplicon.seq().len() as f64) < len_upper
        {
            final_amplicon.push(each_candidate_amplicon);
        }
    }
    final_amplicon.sort_by(|first, second| {
        first
            .calculate_read_quality(false)
            .partial_cmp(&second.calculate_read_quality(false))
            .unwrap()
    });
    final_amplicon.into_iter().take(number).collect()
}

pub fn write_final_amplicon(
    final_amplicon: Vec<SubOwnedRecord>,
    fq_writer: &mut dyn Write,
    fa_writer: &mut dyn Write,
) -> Result<(), anyhow::Error> {
    for each_amplicon in &final_amplicon {
        if each_amplicon.rev_com {
            each_amplicon.write_rev_com(fq_writer)?;
            each_amplicon.write_rev_com_fa(fa_writer)?;
        } else {
            NanoRead::write(each_amplicon, fq_writer)?;
            each_amplicon.write_fa(fa_writer)?;
        }
    }
    Ok(())
}

pub fn mafft_msa<P: AsRef<Path> + Display>(
    msa_input_path: &P,
    msa_output_path: &P,
    mafft_path: &str,
) -> Result<(), anyhow::Error> {
    let cmd = format!(
        "{}  --auto --thread 4 {} > {}",
        mafft_path, msa_input_path, msa_output_path
    );

    let msa_res = std::process::Command::new("/bin/bash")
        .arg("-c")
        .arg(&cmd)
        .output();
    match msa_res {
        Ok(msa) => {
            if !msa.status.success() {
                eprintln!("{}", ansi_term::Color::Red.paint(cmd));
                eprintln!(
                    "{}",
                    ansi_term::Color::Red.paint(std::str::from_utf8(&msa.stderr)?)
                );
                eprintln!(
                    "{}",
                    ansi_term::Color::Red.paint("Multiple Sequence Alignment by mafft failed")
                );
                std::process::exit(1);
            }
        }
        Err(err) => {
            eprintln!("{}", ansi_term::Color::Red.paint(format!("{:?}", err)));
            std::process::exit(1);
        }
    }
    Ok(())
}

fn get_index_of_base(base: &u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b't' | b'T' => 1,
        b'g' | b'G' => 2,
        b'c' | b'C' => 3,
        b'-' => 4,
        other => {
            eprintln!(
                "Expect A/T/C/G/a/t/c/g/-, but found {} in msa result file",
                *other as char
            );
            std::process::exit(1);
        }
    }
}

pub fn get_consensus_from_msa<P: AsRef<Path>>(
    file_path: P,
    read_id: &str,
) -> Result<String, anyhow::Error> {
    let msa_reader = bio::io::fasta::Reader::new(std::fs::File::open(file_path)?);
    let mut records = msa_reader.records();
    let first_record = records
        .next()
        .expect("Empty fasta file found")
        .expect("Bad fasta format found");
    let align_len = first_record.seq().len();
    let mut v = vec![[0usize, 0, 0, 0, 0]; align_len];
    for (idx, base) in first_record.seq().iter().enumerate() {
        v[idx][get_index_of_base(base)] += 1;
    }
    loop {
        if let Some(each_record_res) = records.next() {
            let each_record = each_record_res?;
            for (idx, base) in each_record.seq().iter().enumerate() {
                v[idx][get_index_of_base(base)] += 1;
            }
        } else {
            break;
        }
    }
    let mut consensus_seq = String::new();
    for each_pos in &v {
        let max_count = each_pos.iter().max().unwrap();
        if max_count == &each_pos[0] {
            consensus_seq.push('A')
        } else if max_count == &each_pos[1] {
            consensus_seq.push('T')
        } else if max_count == &each_pos[2] {
            consensus_seq.push('G')
        } else if max_count == &each_pos[3] {
            consensus_seq.push('C')
        } else {
            continue;
        }
    }
    Ok(format!(">{}\n{}\n", read_id, consensus_seq))
}
