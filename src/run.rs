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
        use crate::utils::quit_with_error;
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

        fn stats_receiver(
            receiver: Receiver<RecordSet>,
            gc: bool,
            dont_use_dorado_quality: bool,
        ) -> Vec<EachStats> {
            let mut all_stats: Vec<EachStats> = vec![];
            for record_set in receiver {
                let record_vec = record_set.into_iter().collect::<Vec<RefRecord>>();
                all_stats.extend(
                    record_vec
                        .into_par_iter()
                        .filter_map(|x| x.stats(gc, dont_use_dorado_quality))
                        .collect::<Vec<EachStats>>(),
                );
            }
            all_stats
        }

        pub fn stats<R>(
            reader: R,
            thread: usize,
            gc: bool,
            dont_use_dorado_quality: bool,
        ) -> Vec<EachStats>
        where
            R: Read + Send + Any,
        {
            if thread == 1 {
                FastqReader::new(reader).stats(gc, dont_use_dorado_quality)
            } else {
                let (sender, receiver) = mpsc::sync_channel(1000);
                let _ = thread::spawn(move || {
                    let mut reader = FastqReader::<R>::new(reader);
                    loop {
                        let mut record_set = RecordSet::default();
                        let read_set_result = reader.read_record_set(&mut record_set);
                        if read_set_result.is_none() {
                            break;
                        } else {
                            match read_set_result.unwrap() {
                                Ok(_) => {}
                                Err(err) => quit_with_error(&err.to_string()),
                            }
                        }
                        let x = sender.send(record_set);
                        if x.is_err() {
                            break;
                        }
                    }
                    Result::<(), anyhow::Error>::Ok(())
                });
                stats_receiver(receiver, gc, dont_use_dorado_quality)
            }
        }

        pub fn stats_fastq_dir(
            path: &Path,
            thread: usize,
            gc: bool,
            dont_use_dorado_quality: bool,
        ) -> Vec<EachStats> {
            let fastqs = collect_fastq_dir(path).unwrap();
            fastqs
                .into_par_iter()
                .map(|fq| {
                    if fq.to_str().unwrap().ends_with(".gz") {
                        stats(
                            MultiGzDecoder::new(BufReader::new(File::open(fq).unwrap())),
                            thread,
                            gc,
                            dont_use_dorado_quality,
                        )
                    } else {
                        stats(
                            BufReader::new(File::open(fq).unwrap()),
                            thread,
                            gc,
                            dont_use_dorado_quality,
                        )
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
        use std::io::{BufReader, Read, Write};
        use std::path::Path;
        use std::sync::mpsc;
        use std::sync::mpsc::Receiver;
        use std::thread;

        fn filter_receiver(
            receiver: Receiver<RecordSet>,
            fo: &FilterOption,
            writer: &mut dyn Write,
            retain_failed: bool,
            failed_writer: &mut dyn Write,
        ) -> Result<Vec<(String, u32, f32)>, anyhow::Error> {
            let mut this_receiver_stats = Vec::<(String, u32, f32)>::new();
            for record_set in receiver {
                let record_vec = record_set.into_iter().collect::<Vec<RefRecord>>();
                let vec2 = record_vec
                    .par_iter()
                    .map(|x| x.is_passed(fo))
                    .collect::<Vec<(bool, String, u32, f32)>>();
                if retain_failed {
                    for (ref_record, (is_passed, read_name, read_len, read_qual)) in
                        record_vec.iter().zip(vec2)
                    {
                        if is_passed {
                            this_receiver_stats.push((read_name, read_len, read_qual));
                            NanoRead::write(ref_record, writer)?;
                        } else {
                            NanoRead::write(ref_record, failed_writer)?;
                        }
                    }
                } else {
                    for (ref_record, (is_passed, read_name, read_len, read_qual)) in
                        record_vec.iter().zip(vec2)
                    {
                        if is_passed {
                            this_receiver_stats.push((read_name, read_len, read_qual));
                            NanoRead::write(ref_record, writer)?;
                        }
                    }
                }
            }
            Ok(this_receiver_stats)
        }

        pub fn filter<R>(
            reader: R,
            thread: usize,
            writer: &mut dyn Write,
            fo: &FilterOption,
            retain_failed: bool,
            failed_writer: &mut dyn Write,
        ) -> Result<Vec<(String, u32, f32)>, anyhow::Error>
        where
            R: Read + Send + Any,
        {
            let mut filter_stats = Vec::<(String, u32, f32)>::new();
            if thread == 1 {
                let x =
                    FastqReader::new(reader).filter(writer, fo, retain_failed, failed_writer)?;
                return Ok(x);
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
                let receiver_res =
                    filter_receiver(receiver, fo, writer, retain_failed, failed_writer)?;
                filter_stats.extend(receiver_res);
            }
            Ok(filter_stats)
        }

        pub fn filter_fastq_dir(
            path: &Path,
            thread: usize,
            writer: &mut dyn Write,
            fo: &FilterOption,
            retain_failed: bool,
            failed_writer: &mut dyn Write,
        ) -> Result<Vec<(String, u32, f32)>, anyhow::Error> {
            let mut filter_stats = Vec::<(String, u32, f32)>::new();
            let fastqs = collect_fastq_dir(path).unwrap();
            for fq in fastqs {
                if fq.to_str().unwrap().ends_with(".gz") {
                    let this_file_filter_res = filter(
                        MultiGzDecoder::new(BufReader::new(File::open(fq)?)),
                        thread,
                        writer,
                        fo,
                        retain_failed,
                        failed_writer,
                    )?;
                    filter_stats.extend(this_file_filter_res);
                } else {
                    let this_file_filter_res = filter(
                        BufReader::new(File::open(fq)?),
                        thread,
                        writer,
                        fo,
                        retain_failed,
                        failed_writer,
                    )?;
                    filter_stats.extend(this_file_filter_res);
                }
            }
            Ok(filter_stats)
        }
    }


}

pub mod run_entry {
    use crate::amplicon::{
        filter_candidate_amplicon, get_candidate_amplicon, get_consensus_from_msa, mafft_msa,
        write_final_amplicon,
    };
    use crate::fastq::{FastqReader, FilterOption, NanoRead};
    use crate::input_type::{InputType, check_input_type};
    use crate::run::sub_run::filter::{filter, filter_fastq_dir};
    use crate::run::{LOCAL_ALIGNER, collect_fastq_dir};
    use crate::utils::{quit_with_error, rev_com};
    use clap::ArgMatches;
    use flate2::bufread::MultiGzDecoder;
    use seq_io::fastq::Record;
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Read, Write};
    use std::path::Path;


    pub fn run_filter(filter_cmd: &ArgMatches) -> Result<(), anyhow::Error> {
        let input = filter_cmd.get_one::<String>("input");
        let output = filter_cmd.get_one::<String>("output");
        let min_len = filter_cmd.get_one::<usize>("min_len").unwrap();
        let max_len = filter_cmd.get_one::<usize>("max_len").unwrap();
        let dont_use_dorado_quality = filter_cmd.get_flag("dont_use_dorado_quality");
        let min_qual = filter_cmd.get_one::<f64>("min_qual").unwrap();
        let max_qual = filter_cmd.get_one::<f64>("max_qual").unwrap();
        let gc = filter_cmd.get_flag("gc");
        let min_gc = filter_cmd.get_one::<f64>("min_gc").unwrap();
        let max_gc = filter_cmd.get_one::<f64>("max_gc").unwrap();
        let thread = filter_cmd.get_one::<u16>("thread").unwrap();
        let max_bases = filter_cmd.get_one::<u64>("max_bases");
        let failed_fq_path = filter_cmd.get_one::<String>("retain_failed");
        let filter_option = FilterOption {
            min_len: *min_len as u32,
            max_len: *max_len as u32,
            dont_use_dorado_quality,
            min_qual: *min_qual as f32,
            max_qual: *max_qual as f32,
            gc,
            min_gc: *min_gc as f32,
            max_gc: *max_gc as f32,
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

        let tmp_filter_fastq = format!("./.{}.tmp.filter.fastq", uuid::Uuid::new_v4());
        let mut filter_stats = {
            let mut tmp_writer = BufWriter::new(File::create(&tmp_filter_fastq)?);
            let input_t = check_input_type(input.unwrap());
            match input_t {
                // InputType::FastqFromStdin => filter(
                //     std::io::stdin(),
                //     *thread as usize,
                //     &mut tmp_writer,
                //     &filter_option,
                //     failed_retain,
                //     &mut failed_writer,
                // )?,
                InputType::OneFastqFile => {
                    let reader = BufReader::new(File::open(input.unwrap())?);
                    filter(
                        reader,
                        *thread as usize,
                        &mut tmp_writer,
                        &filter_option,
                        failed_retain,
                        &mut failed_writer,
                    )?
                }
                InputType::OneFastqGzippedFile => {
                    let reader = MultiGzDecoder::new(BufReader::new(File::open(input.unwrap())?));
                    filter(
                        reader,
                        *thread as usize,
                        &mut tmp_writer,
                        &filter_option,
                        failed_retain,
                        &mut failed_writer,
                    )?
                }
                InputType::DirectoryContainFastqsOrFastqsGzipped => filter_fastq_dir(
                    std::path::Path::new(input.unwrap()),
                    *thread as usize,
                    &mut tmp_writer,
                    &filter_option,
                    failed_retain,
                    &mut failed_writer,
                )?,
                _ => quit_with_error(
                    "Bad input type, only 1). fastq from stdin; 2). fastq[.gz] file; 3). directory containing some fastq[.gz] files supported. Check --input",
                ),
            }
        };

        fn mv_tmp_to_final(old_path: &str, final_path: &str) -> Result<(), anyhow::Error> {
            let cmd_process = std::process::Command::new("mv")
                .arg(old_path)
                .arg(final_path)
                .output();
            match cmd_process {
                Ok(output) => {
                    if !output.status.success() {
                        println!("{}", std::str::from_utf8(&output.stderr)?);
                        std::process::exit(1);
                    }
                }
                Err(error) => {
                    println!("Get final filtered fastq failed");
                    println!("{:?}", error);
                }
            }
            Ok(())
        }

        match max_bases {
            None => {
                if output.is_none() {
                    let mut tmp_reader =
                        std::io::BufReader::new(std::fs::File::open(&tmp_filter_fastq)?);
                    let mut buf = [0; 1024 * 8];
                    loop {
                        let bytes_size = tmp_reader.read(&mut buf)?;
                        if bytes_size == 0 {
                            break;
                        }
                        std::io::stdout().write_all(&buf[..bytes_size])?;
                    }
                    std::fs::remove_file(&tmp_filter_fastq)?;
                } else {
                    mv_tmp_to_final(&tmp_filter_fastq, output.unwrap())?;
                }
            }
            Some(target_bases_count) => {
                let target_bases_count = *target_bases_count as usize;
                let mut total_filter_bases = 0usize;
                for each in &filter_stats {
                    total_filter_bases += each.1 as usize
                }
                if total_filter_bases < target_bases_count {
                    if output.is_none() {
                        let mut tmp_reader =
                            std::io::BufReader::new(std::fs::File::open(&tmp_filter_fastq)?);
                        let mut buf = [0; 1024 * 8];
                        loop {
                            let bytes_size = tmp_reader.read(&mut buf)?;
                            if bytes_size == 0 {
                                break;
                            }
                            std::io::stdout().write_all(&buf[..bytes_size])?
                        }
                        std::fs::remove_file(tmp_filter_fastq)?;
                    } else {
                        mv_tmp_to_final(&tmp_filter_fastq, output.unwrap())?;
                    }
                } else {
                    let mut read_names = HashSet::<String>::new();
                    filter_stats.sort_by(|x, y| y.2.partial_cmp(&x.2).unwrap());
                    let mut total_retain_bases = 0usize;
                    for each in filter_stats {
                        total_retain_bases += each.1 as usize;
                        read_names.insert(each.0);
                        if total_retain_bases >= target_bases_count {
                            break;
                        }
                    }
                    let tmp_reader =
                        std::io::BufReader::new(std::fs::File::open(&tmp_filter_fastq)?);
                    let mut tmp_fastq_reader = FastqReader::new(tmp_reader);
                    let mut writer: Box<dyn Write> = if output.is_none() {
                        Box::new(BufWriter::new(std::io::stdout()))
                    } else {
                        Box::new(BufWriter::new(File::create(output.unwrap())?))
                    };
                    loop {
                        if let Some(each_record_res) = tmp_fastq_reader.next() {
                            let each_record = each_record_res?;
                            if read_names.contains(each_record.id()?) {
                                NanoRead::write(&each_record, &mut writer)?
                            } else {
                                NanoRead::write(&each_record, &mut failed_writer)?
                            }
                        } else {
                            break;
                        }
                    }
                    std::fs::remove_file(&tmp_filter_fastq)?;
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
        let range = amplicon_cmd.get_one::<u8>("range").unwrap();
        let mafft = amplicon_cmd.get_one::<String>("mafft").unwrap();
        let number = *amplicon_cmd.get_one::<u32>("number").unwrap() as usize;
        let name = amplicon_cmd.get_one::<String>("name").unwrap();
        let output_dir = Path::new(output);
        if !output_dir.exists() {
            std::fs::create_dir_all(&output_dir)?;
        }
        let input_t = check_input_type(input.unwrap());
        let candidate_amplicon = match input_t {
            // InputType::FastqFromStdin => {
            //     get_candidate_amplicon(Option::<&String>::None, fwd, rev, est_len, *range)?
            // }
            InputType::OneFastqGzippedFile | InputType::OneFastqFile => {
                get_candidate_amplicon(input, fwd, rev, est_len, *range)?
            }
            InputType::DirectoryContainFastqsOrFastqsGzipped => {
                let fqs = collect_fastq_dir(std::path::Path::new(input.unwrap()))?;
                debug_assert!(fqs.len() > 0);
                fqs.iter()
                    .map(|fq| {
                        get_candidate_amplicon(Some(fq), fwd, rev, est_len, *range).expect(
                            &format!(
                                "Get candidate amplicon from {} failed",
                                fq.to_str().unwrap()
                            ),
                        )
                    })
                    .flatten()
                    .collect::<Vec<_>>()
            }
            _ => quit_with_error(
                "Bad input type, only 1). fastq from stdin; 2). fastq[.gz] file; 3). directory containing some fastq[.gz] files supported. Check --input",
            ),
        };
        let final_amplicon = filter_candidate_amplicon(candidate_amplicon, number);
        let fq_file_path = output_dir.join("candidate_amplicon.fastq");
        let fa_file_path = output_dir.join("candidate_amplicon.fasta");
        let mafft_output_path = output_dir.join("msa.fasta");
        {
            let mut fq_writer = std::io::BufWriter::new(std::fs::File::create(&fq_file_path)?);
            let mut fa_writer = std::io::BufWriter::new(std::fs::File::create(&fa_file_path)?);
            write_final_amplicon(final_amplicon, &mut fq_writer, &mut fa_writer)?;
        }
        mafft_msa(
            &fa_file_path.to_str().unwrap(),
            &mafft_output_path.to_str().unwrap(),
            mafft,
        )?;
        let consensus_fasta_content = get_consensus_from_msa(&mafft_output_path, name)?;
        std::fs::write(
            output_dir.join(format!("{}.fasta", name)),
            consensus_fasta_content,
        )?;
        Ok(())
    }

}
