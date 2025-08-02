use crate::alignment::LocalAligner;
use std::cell::RefCell;
use std::path::{Path, PathBuf};

thread_local! {
    static LOCAL_ALIGNER: RefCell<LocalAligner>= RefCell::new(LocalAligner::default());
}

fn collect_fastq_dir(path: &Path) -> Result<Vec<PathBuf>, anyhow::Error> {
    assert!(path.is_dir());
    let all_fqs = path
        .read_dir()
        .expect(&format!(
            "read directory: {} failed",
            path.to_str().unwrap()
        ))
        .filter_map(|x| {
            if let Ok(fs) = x {
                let fs_path = fs.path();
                let fs_path_str = fs_path.to_str()?;
                return if fs_path_str.ends_with(".fastq")
                    || fs_path_str.ends_with(".fq")
                    || fs_path_str.ends_with(".fastq.gz")
                    || fs_path_str.ends_with(".fq.gz")
                {
                    Some(fs_path)
                } else {
                    None
                };
            }
            None
        })
        .collect::<Vec<PathBuf>>();
    Ok(all_fqs)
}

mod sub_run {
    pub mod stats {
        use crate::fastq::{EachStats, FastqReader, NanoRead};
        use crate::run::collect_fastq_dir;
        use flate2::bufread::MultiGzDecoder;
        use rayon::prelude::*;
        use seq_io::fastq::{RecordSet, RefRecord};
        use std::any::Any;
        use std::fs::File;
        use std::io::{BufReader, Read};
        use std::path::Path;
        use std::sync::mpsc;
        use std::sync::mpsc::Receiver;
        use std::thread;

        fn stats_receiver(receiver: Receiver<RecordSet>, gc: bool) -> Vec<EachStats> {
            let mut all_stats: Vec<EachStats> = vec![];
            for record_set in receiver {
                let record_vec = record_set.into_iter().collect::<Vec<RefRecord>>();
                all_stats.extend(
                    record_vec
                        .into_par_iter()
                        .map(|x| x.stats(gc))
                        .collect::<Vec<EachStats>>(),
                );
            }
            all_stats
        }

        pub fn stats<R>(reader: R, thread: usize, gc: bool) -> Vec<EachStats>
        where
            R: Read + Send + Any,
        {
            if thread == 1 {
                FastqReader::new(reader).stats(gc)
            } else {
                let (sender, receiver) = mpsc::sync_channel(1000);
                let _ = thread::spawn(move || {
                    let mut reader = FastqReader::<R>::new(reader);
                    loop {
                        let mut record_set = RecordSet::default();
                        if reader.read_record_set(&mut record_set).is_none() {
                            break;
                        }
                        if sender.send(record_set).is_err() {
                            break;
                        }
                    }
                    Result::<(), anyhow::Error>::Ok(())
                });
                stats_receiver(receiver, gc)
            }
        }

        pub fn stats_fastq_dir(path: &Path, thread: usize, gc: bool) -> Vec<EachStats> {
            let fastqs = collect_fastq_dir(path).unwrap();
            fastqs
                .into_par_iter()
                .map(|fq| {
                    if fq.to_str().unwrap().ends_with(".gz") {
                        stats(
                            MultiGzDecoder::new(BufReader::new(File::open(fq).unwrap())),
                            thread,
                            gc,
                        )
                    } else {
                        stats(BufReader::new(File::open(fq).unwrap()), thread, gc)
                    }
                })
                .flatten()
                .collect::<Vec<EachStats>>()
        }
    }
    pub mod filter {
        use crate::fastq::{FastqReader, FilterOption, NanoRead};
        use crate::run::collect_fastq_dir;
        use flate2::bufread::MultiGzDecoder;
        use rayon::prelude::*;
        use seq_io::fastq::{RecordSet, RefRecord};
        use std::any::Any;
        use std::fs::File;
        use std::io::{BufReader, BufWriter, Read, Write};
        use std::path::Path;
        use std::sync::mpsc;
        use std::sync::mpsc::Receiver;
        use std::thread;

        fn filter_receiver(
            receiver: Receiver<RecordSet>,
            fo: &FilterOption,
            writer: &mut dyn Write,
            retain_failed: bool,
            failed_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error> {
            for record_set in receiver {
                let record_vec = record_set.into_iter().collect::<Vec<RefRecord>>();
                let vec2 = record_vec
                    .par_iter()
                    .map(|x| x.is_passed(fo))
                    .collect::<Vec<bool>>();
                if retain_failed {
                    for (ref_record, is_passed) in record_vec.iter().zip(&vec2) {
                        if *is_passed {
                            NanoRead::write(ref_record, writer)?;
                        } else {
                            NanoRead::write(ref_record, failed_writer)?;
                        }
                    }
                } else {
                    for (ref_record, is_passed) in record_vec.iter().zip(&vec2) {
                        if *is_passed {
                            NanoRead::write(ref_record, writer)?;
                        }
                    }
                }
            }
            Ok(())
        }

        pub fn filter<R>(
            reader: R,
            thread: usize,
            writer: &mut dyn Write,
            fo: &FilterOption,
            retain_failed: bool,
            failed_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error>
        where
            R: Read + Send + Any,
        {
            if thread == 1 {
                FastqReader::new(reader).filter(writer, fo, retain_failed, failed_writer)?
            } else {
                let (sender, receiver) = mpsc::sync_channel(1000);
                let _ = thread::spawn(move || {
                    let mut reader = FastqReader::<R>::new(reader);
                    loop {
                        let mut record_set = RecordSet::default();
                        if reader.read_record_set(&mut record_set).is_none() {
                            break;
                        }
                        if sender.send(record_set).is_err() {
                            break;
                        }
                    }
                });
                filter_receiver(receiver, fo, writer, retain_failed, failed_writer)?
            }
            Ok(())
        }

        pub fn filter_fastq_dir(
            path: &Path,
            thread: usize,
            writer: &mut dyn Write,
            fo: &FilterOption,
            retain_failed: bool,
            failed_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error> {
            let fastqs = collect_fastq_dir(path).unwrap();
            for fq in fastqs {
                if fq.to_str().unwrap().ends_with(".gz") {
                    filter(
                        MultiGzDecoder::new(BufReader::new(File::open(fq)?)),
                        thread,
                        writer,
                        fo,
                        retain_failed,
                        failed_writer,
                    )?
                } else {
                    filter(
                        BufReader::new(File::open(fq)?),
                        thread,
                        writer,
                        fo,
                        retain_failed,
                        failed_writer,
                    )?;
                }
            }
            Ok(())
        }
    }

    pub mod trim {
        use crate::fastq::{FastqReader, NanoRead};
        use crate::run::{LOCAL_ALIGNER, collect_fastq_dir};
        use crate::trim::adapter::TrimConfig;
        use flate2::bufread::MultiGzDecoder;
        use rayon::prelude::*;
        use seq_io::fastq::{Record, RecordSet, RefRecord};
        use std::any::Any;
        use std::fs::File;
        use std::io::{BufReader, BufWriter, Read, Write};
        use std::path::Path;
        use std::sync::mpsc;
        use std::sync::mpsc::Receiver;
        use std::thread;

        fn trim_receiver(
            receiver: Receiver<RecordSet>,
            trim_cfg: &TrimConfig,
            writer: &mut dyn Write,
            min_len: usize,
            pretty_log: bool,
            log_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error> {
            for record_set in receiver {
                let record_vec = record_set.into_iter().collect::<Vec<RefRecord>>();
                let trimmed_result: Vec<_> = record_vec
                    .par_iter()
                    .map(|each_ref_record| {
                        LOCAL_ALIGNER.with_borrow_mut(|local_aligner| {
                            each_ref_record.trim(trim_cfg, local_aligner, min_len, pretty_log, true)
                        })
                    })
                    .collect();
                if pretty_log {
                    trimmed_result.iter().zip(record_vec.iter()).for_each(
                        |((trimmed_info_opt, log_string), ref_record)| {
                            if let Some((sub_seq, sub_qual)) = trimmed_info_opt {
                                write!(
                                    writer,
                                    "@{}\n{}\n+\n{}\n",
                                    unsafe { std::str::from_utf8_unchecked(ref_record.head()) },
                                    unsafe { std::str::from_utf8_unchecked(sub_seq) },
                                    unsafe { std::str::from_utf8_unchecked(sub_qual) }
                                )
                                .expect(&format!(
                                    "write trimmed fastq into output file failed for {}",
                                    std::str::from_utf8(ref_record.head()).unwrap()
                                ));
                            }
                            write!(
                                log_writer,
                                "{}",
                                log_string.as_deref().expect(&format!(
                                    "write trimmed fastq into output file failed for {}",
                                    std::str::from_utf8(ref_record.head()).unwrap()
                                ))
                            )
                            .unwrap();
                        },
                    )
                } else {
                    trimmed_result.iter().zip(record_vec.iter()).for_each(
                        |((trimmed_info_opt, _), ref_record)| {
                            if let Some((sub_seq, sub_qual)) = trimmed_info_opt {
                                write!(
                                    writer,
                                    "@{}\n{}\n+\n{}\n",
                                    unsafe { std::str::from_utf8_unchecked(ref_record.head()) },
                                    unsafe { std::str::from_utf8_unchecked(sub_seq) },
                                    unsafe { std::str::from_utf8_unchecked(sub_qual) }
                                )
                                .expect(&format!(
                                    "write trimmed fastq into output file failed for {}",
                                    std::str::from_utf8(ref_record.head()).unwrap()
                                ));
                            }
                        },
                    )
                }
            }
            Ok(())
        }

        pub fn trim<R>(
            reader: R,
            writer: &mut dyn Write,
            trim_cfg: &TrimConfig,
            min_len: usize,
            pretty_log: bool,
            log_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error>
        where
            R: Read + Send + Any,
        {
            let (sender, receiver) = mpsc::sync_channel(1000);
            let _ = thread::spawn(move || {
                let mut reader = FastqReader::<R>::new(reader);
                loop {
                    let mut record_set = RecordSet::default();
                    if reader.read_record_set(&mut record_set).is_none() {
                        break;
                    }
                    if sender.send(record_set).is_err() {
                        break;
                    }
                }
            });
            trim_receiver(receiver, trim_cfg, writer, min_len, pretty_log, log_writer)?;
            Ok(())
        }

        pub fn trim_fastq_dir(
            path: &Path,
            writer: &mut dyn Write,
            trim_cfg: &TrimConfig,
            min_len: usize,
            pretty_log: bool,
            log_writer: &mut BufWriter<File>,
        ) -> Result<(), anyhow::Error> {
            let fastqs = collect_fastq_dir(path)?;
            for fq in fastqs {
                if fq.to_str().unwrap().ends_with(".gz") {
                    trim(
                        MultiGzDecoder::new(BufReader::new(File::open(fq)?)),
                        writer,
                        trim_cfg,
                        min_len,
                        pretty_log,
                        log_writer,
                    )?
                } else {
                    trim(
                        BufReader::new(File::open(fq)?),
                        writer,
                        trim_cfg,
                        min_len,
                        pretty_log,
                        log_writer,
                    )?
                }
            }
            Ok(())
        }
    }
}

pub mod run_entry {
    use crate::alignment::{LocalAligner, Scores};
    use crate::amplicon::{
        QueryAmpMode, filter_candidate_amplicon, get_candidate_amplicon, get_consensus_from_msa,
        mafft_msa, write_final_amplicon,
    };
    use crate::fastq::{FilterOption};
    use crate::run::LOCAL_ALIGNER;
    use crate::run::sub_run::filter::{filter, filter_fastq_dir};
    use crate::run::sub_run::stats::{stats, stats_fastq_dir};
    use crate::run::sub_run::trim::{trim, trim_fastq_dir};
    use crate::summary::{make_plot, write_stats, write_summary};
    use crate::trim::adapter::{TrimConfig, get_trim_cfg};
    use crate::utils::rev_com;
    use clap::ArgMatches;
    use flate2::bufread::MultiGzDecoder;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write };
    use std::path::{Path};

    pub fn run_stats(stats_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
        let input = stats_cmd.get_one::<String>("input");
        let output = stats_cmd.get_one::<String>("output");
        let summary = stats_cmd.get_one::<String>("summary").unwrap();
        let topn = stats_cmd.get_one::<u16>("topn").unwrap();
        let quality = stats_cmd.get_one::<Vec<f64>>("quality").unwrap();
        let lengths = stats_cmd.get_one::<Vec<usize>>("length");
        let gc = stats_cmd.get_flag("gc");
        let thread = stats_cmd.get_one::<u16>("thread").unwrap();
        let plot = stats_cmd.get_one::<String>("plot");
        let python = stats_cmd.get_one::<String>("python").unwrap();
        let quantile = stats_cmd.get_one::<f64>("quantile").unwrap();
        let format = stats_cmd
            .get_many::<String>("format")
            .unwrap()
            .collect::<Vec<&String>>();

        rayon::ThreadPoolBuilder::new()
            .num_threads(*thread as usize)
            .build_global()?;

        let mut stats_result = match input {
            // None => stats_result = stats_stdin(*thread as usize, gc),
            None => stats(std::io::stdin(), *thread as usize, gc),
            Some(input) => {
                let input_path = Path::new(input);
                if input_path.is_file() {
                    if input_path.to_str().unwrap().ends_with(".gz") {
                        stats(
                            MultiGzDecoder::new(BufReader::new(File::open(input)?)),
                            *thread as usize,
                            gc,
                        )
                    } else {
                        stats(BufReader::new(File::open(input)?), *thread as usize, gc)
                    }
                } else {
                    stats_fastq_dir(input_path, *thread as usize, gc)
                }
            }
        };

        let tmp_stats_outfile = format!("/tmp/NanofqStatsTmpResult_{}.tsv", uuid::Uuid::new_v4());
        match output {
            None => {
                write_stats(&stats_result, &mut std::io::stdout(), gc)?;
                if plot.is_some() {
                    let mut writer = std::fs::File::create(&tmp_stats_outfile)?;
                    write_stats(&stats_result, &mut writer, gc)?;
                }
            }
            Some(output_file) => write_stats(
                &stats_result,
                &mut std::io::BufWriter::new(File::create(output_file).unwrap()),
                gc,
            )?,
        }
        let basic_stats =
            write_summary(&mut stats_result, lengths, quality, *topn as usize, summary);
        let formats = format
            .iter()
            .map(|x| (**x).clone())
            .collect::<Vec<String>>();
        if plot.is_some() {
            if output.is_none() {
                make_plot(
                    &basic_stats,
                    *quantile,
                    plot.unwrap(),
                    &formats,
                    python,
                    &tmp_stats_outfile,
                )?;
            } else {
                make_plot(
                    &basic_stats,
                    *quantile,
                    plot.unwrap(),
                    &formats,
                    python,
                    output.unwrap(),
                )?;
            }
        }
        Ok(())
    }

    pub fn run_filter(filter_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
        let input = filter_cmd.get_one::<String>("input");
        let output = filter_cmd.get_one::<String>("output");
        let min_len = filter_cmd.get_one::<usize>("min_len").unwrap();
        let max_len = filter_cmd.get_one::<usize>("max_len").unwrap();
        let min_qual = filter_cmd.get_one::<f64>("min_qual").unwrap();
        let max_qual = filter_cmd.get_one::<f64>("max_qual").unwrap();
        let gc = filter_cmd.get_flag("gc");
        let min_gc = filter_cmd.get_one::<f64>("min_gc").unwrap();
        let max_gc = filter_cmd.get_one::<f64>("max_gc").unwrap();
        let thread = filter_cmd.get_one::<u16>("thread").unwrap();
        let failed_fq_path = filter_cmd.get_one::<String>("retain_failed");
        let filter_option = FilterOption {
            min_len: *min_len,
            max_len: *max_len,
            min_qual: *min_qual,
            max_qual: *max_qual,
            gc,
            min_gc: *min_gc,
            max_gc: *max_gc,
            retain_failed: failed_fq_path,
        };
        let failed_retain = if failed_fq_path.is_none() {
            false
        } else {
            true
        };
        let mut failed_writer = filter_option.set_failed_fastq_file()?.unwrap();
        rayon::ThreadPoolBuilder::new()
            .num_threads(*thread as usize)
            .build_global()?;

        let mut writer: Box<dyn Write> = if output.is_none() {
            Box::new(BufWriter::new(std::io::stdout()))
        } else {
            Box::new(BufWriter::new(File::create(output.unwrap())?))
        };

        match input {
            None => {
                filter(
                    std::io::stdin(),
                    *thread as usize,
                    &mut writer,
                    &filter_option,
                    failed_retain,
                    &mut failed_writer,
                )?;
            }
            Some(input_path) => {
                let ends_with_gz = input_path.ends_with(".gz");
                let input_path = Path::new(input_path);
                if input_path.is_file() {
                    if ends_with_gz {
                        let reader = MultiGzDecoder::new(BufReader::new(File::open(input_path)?));
                        filter(
                            reader,
                            *thread as usize,
                            &mut writer,
                            &filter_option,
                            failed_retain,
                            &mut failed_writer,
                        )?
                    } else {
                        let reader = BufReader::new(File::open(input_path)?);
                        filter(
                            reader,
                            *thread as usize,
                            &mut writer,
                            &filter_option,
                            failed_retain,
                            &mut failed_writer,
                        )?
                    }
                } else {
                    filter_fastq_dir(
                        input_path,
                        *thread as usize,
                        &mut writer,
                        &filter_option,
                        failed_retain,
                        &mut failed_writer,
                    )?
                }
            }
        }
        Ok(())
    }
    pub fn run_trim(trim_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
        let input = trim_cmd.get_one::<String>("input");
        let output = trim_cmd.get_one::<String>("output");
        let primers = trim_cmd.get_one::<String>("primers");
        let kit = trim_cmd.get_one::<String>("kit");
        let min_len = *trim_cmd.get_one::<u32>("min_len").unwrap() as usize;
        let thread = *trim_cmd.get_one::<u16>("thread").unwrap();
        let match_ = *trim_cmd.get_one::<i32>("match").unwrap();
        let mismatch = *trim_cmd.get_one::<i32>("mismatch").unwrap();
        let gap_open = *trim_cmd.get_one::<i32>("gap_opened").unwrap();
        let gap_extend = *trim_cmd.get_one::<i32>("gap_extend").unwrap();
        let end5_len = (
            trim_cmd.value_source("end5_len").unwrap(),
            *trim_cmd.get_one::<u32>("end5_len").unwrap() as usize,
        );
        let end5_pct = (
            trim_cmd.value_source("end5_align_pct").unwrap(),
            *trim_cmd.get_one::<f64>("end5_align_pct").unwrap(),
        );
        let end5_align_ident = (
            trim_cmd.value_source("end5_align_ident").unwrap(),
            *trim_cmd.get_one::<f64>("end5_align_ident").unwrap(),
        );

        let end3_len = (
            trim_cmd.value_source("end3_len").unwrap(),
            *trim_cmd.get_one::<u32>("end3_len").unwrap() as usize,
        );
        let end3_pct = (
            trim_cmd.value_source("end3_align_pct").unwrap(),
            *trim_cmd.get_one::<f64>("end3_align_pct").unwrap(),
        );
        let end3_align_ident = (
            trim_cmd.value_source("end3_align_ident").unwrap(),
            *trim_cmd.get_one::<f64>("end3_align_ident").unwrap(),
        );

        let end5_len_rc = (
            trim_cmd.value_source("end5_len_rc").unwrap(),
            *trim_cmd.get_one::<u32>("end5_len_rc").unwrap() as usize,
        );
        let end5_pct_rc = (
            trim_cmd.value_source("end5_align_pct_rc").unwrap(),
            *trim_cmd.get_one::<f64>("end5_align_pct_rc").unwrap(),
        );
        let end5_align_ident_rc = (
            trim_cmd.value_source("end5_align_ident_rc").unwrap(),
            *trim_cmd.get_one::<f64>("end5_align_ident_rc").unwrap(),
        );

        let end3_len_rc = (
            trim_cmd.value_source("end3_len_rc").unwrap(),
            *trim_cmd.get_one::<u32>("end3_len_rc").unwrap() as usize,
        );
        let end3_pct_rc = (
            trim_cmd.value_source("end3_align_pct_rc").unwrap(),
            *trim_cmd.get_one::<f64>("end3_align_pct_rc").unwrap(),
        );
        let end3_align_ident_rc = (
            trim_cmd.value_source("end3_align_ident_rc").unwrap(),
            *trim_cmd.get_one::<f64>("end3_align_ident_rc").unwrap(),
        );
        let scores = Scores {
            match_,
            mismatch,
            gap_open,
            gap_extend,
        };
        let rev_com_used = !trim_cmd.get_flag("rev_com_not_used");
        let mut real_primers = String::new();
        let mut fwd_primer_rc = String::new();
        let mut rev_primer_rc = String::new();
        let mut all_kit_trim_cfg = get_trim_cfg();
        let mut trim_cfg = TrimConfig::default();
        let ref_trim_cfg = if let Some(primers) = primers {
            real_primers = primers.clone();
            let fields: Vec<_> = real_primers
                .splitn(2, ",")
                .into_iter()
                .filter(|x| x.len() > 1)
                .collect();
            let fwd_primer = fields[0];
            let rev_primer = fields[1];
            fwd_primer_rc = rev_com(fwd_primer);
            rev_primer_rc = rev_com(rev_primer);
            trim_cfg = TrimConfig {
                kit_name: "customer",
                end5: Some((fwd_primer, (end5_len.1, end5_pct.1, end5_align_ident.1))),
                end3: Some((&rev_primer_rc, (end3_len.1, end3_pct.1, end3_align_ident.1))),
                rev_com_end5: if rev_com_used {
                    Some((
                        &rev_primer,
                        (end5_len_rc.1, end5_pct_rc.1, end5_align_ident_rc.1),
                    ))
                } else {
                    None
                },
                rev_com_end3: if rev_com_used {
                    Some((
                        &fwd_primer_rc,
                        (end3_len_rc.1, end3_pct_rc.1, end3_align_ident_rc.1),
                    ))
                } else {
                    None
                },
            };
            &trim_cfg
        } else {
            let kit_str = kit.unwrap().as_str();
            let kit_trim_cfg = all_kit_trim_cfg.get_mut(kit_str).unwrap();
            kit_trim_cfg.update(
                end5_len,
                end5_pct,
                end5_align_ident,
                end3_len,
                end3_pct,
                end3_align_ident,
                end5_len_rc,
                end5_pct_rc,
                end5_align_ident_rc,
                end3_len_rc,
                end3_pct_rc,
                end3_align_ident_rc,
            );
            kit_trim_cfg
        };
        let (row, col) = ref_trim_cfg.get_dim();
        let log = trim_cmd.get_one::<String>("log");
        let pretty_log = if log.is_some() { true } else { false };
        let mut log_writer = match log {
            None => BufWriter::new(File::create("/tmp/NanoFqTrimmed.log")?),
            Some(log_file) => BufWriter::new(File::create(log_file)?),
        };
        if pretty_log {
            write!(log_writer, "{}", ref_trim_cfg.get_info())?
        }
        if thread == 1 {
            LOCAL_ALIGNER.with_borrow_mut(|local_aligner: &mut LocalAligner| {
                local_aligner.update((row + 1, col + 1), scores)
            });
        }
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread as usize)
            .start_handler(move |_| {
                LOCAL_ALIGNER.with_borrow_mut(|local_aligner: &mut LocalAligner| {
                    local_aligner.update((row + 1, col + 1), scores)
                })
            })
            .build_global()?;
        let mut writer: Box<dyn Write> = match output {
            None => Box::new(BufWriter::new(std::io::stdout())),
            Some(output_file_path) => Box::new(BufWriter::new(File::create(output_file_path)?)),
        };

        match input {
            None => trim(
                std::io::stdin(),
                &mut writer,
                ref_trim_cfg,
                min_len,
                pretty_log,
                &mut log_writer,
            )?,
            Some(input_path) => {
                let ends_with_gz = input_path.ends_with(".gz");
                let input_path = Path::new(input_path);
                if input_path.is_file() {
                    if ends_with_gz {
                        let reader = MultiGzDecoder::new(BufReader::new(File::open(input_path)?));
                        trim(
                            reader,
                            &mut writer,
                            ref_trim_cfg,
                            min_len,
                            pretty_log,
                            &mut log_writer,
                        )?
                    } else {
                        let reader = BufReader::new(File::open(input_path)?);
                        trim(
                            reader,
                            &mut writer,
                            ref_trim_cfg,
                            min_len,
                            pretty_log,
                            &mut log_writer,
                        )?
                    }
                } else {
                    trim_fastq_dir(
                        input_path,
                        &mut writer,
                        ref_trim_cfg,
                        min_len,
                        pretty_log,
                        &mut log_writer,
                    )?
                }
            }
        }
        Ok(())
    }

    pub fn run_amplicon(amplicon_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
        let input = amplicon_cmd.get_one::<String>("input");
        let output = amplicon_cmd.get_one::<String>("output").unwrap();
        let fwd = amplicon_cmd.get_one::<String>("fwd").unwrap();
        let rev = amplicon_cmd.get_one::<String>("rev").unwrap();
        let est_len = *amplicon_cmd.get_one::<u32>("len").unwrap() as usize;
        let mafft = amplicon_cmd.get_one::<String>("mafft").unwrap();
        let number = *amplicon_cmd.get_one::<u32>("number").unwrap() as usize;
        let name = amplicon_cmd.get_one::<String>("name").unwrap();
        let x = amplicon_cmd.get_one::<String>("mode").unwrap();
        let mode = if x == "find" {
            QueryAmpMode::Find
        } else {
            QueryAmpMode::Align
        };
        let output_dir = Path::new(output);
        if !output_dir.exists() {
            std::fs::create_dir_all(&output_dir)?;
        }
        let candidate_amplicon = get_candidate_amplicon(input, fwd, rev, est_len, &mode)?;
        let final_amplicon = filter_candidate_amplicon(candidate_amplicon, number);
        let fq_file_path = output_dir.join(format!("{}_amplicon.fastq", name));
        let fa_file_path = output_dir.join(format!("{}_amplicon.fasta", name));
        let mafft_output_path = output_dir.join(format!("{}_amplicon.fasta", name));
        {
            let mut fq_writer = std::io::BufWriter::new(std::fs::File::create(&fq_file_path)?);
            let mut fa_writer = std::io::BufWriter::new(std::fs::File::create(&fa_file_path)?);
            write_final_amplicon(final_amplicon, &mut fq_writer, &mut fa_writer)?;
        }
        mafft_msa(&fa_file_path, &mafft_output_path, mafft)?;
        get_consensus_from_msa(&mafft_output_path, name)?;
        Ok(())
    }
}
