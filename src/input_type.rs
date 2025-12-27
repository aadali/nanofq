use crate::utils::quit_with_error;
use rust_htslib::bam::{self, Read};
use std::path::Path;

pub enum InputType {
    FastqFromStdin,
    DirectoryContainFastqsOrFastqsGzipped,
    OneFastqFile,
    OneFastqGzippedFile,
    OneBamOrSamFromStdin,
    OneSamFile,
    UnsortedBam,
    UnalignedBam,
    SortedUnindexedBam,
    IndexedBam,
    WrongType(String),
}

fn check_bam_type(bam_file: &str) -> InputType {
    let bam_reader_res = bam::Reader::from_path(bam_file);
    if bam_reader_res.is_err() {
        quit_with_error(&format!(
            "Check bam type failed for {} when making bam::Reader",
            bam_file
        ));
    }
    let bam_reader = bam_reader_res.unwrap();
    let header_view = bam_reader.header();
    if header_view.target_count() == 0 {
        return InputType::UnalignedBam;
    } else {
        let header_hashmap = bam::Header::from_template(header_view).to_hashmap();
        let hd_opt = header_hashmap.get(&"HD".to_string());
        if hd_opt.is_none() {
            return InputType::UnsortedBam;
        } else {
            let hd_linear_map = hd_opt.unwrap().first().unwrap();
            let so_opt = hd_linear_map.get(&"SO".to_string());
            match so_opt {
                None => InputType::UnsortedBam,
                Some(so) => {
                    if so == "coordinate" {
                        let bai_index_path = std::path::PathBuf::from(&format!("{}.bai", bam_file));
                        let csi_index_path = std::path::PathBuf::from(&format!("{}.csi", bam_file));
                        if bai_index_path.exists() || csi_index_path.exists() {
                            let bam_file_meta_res = Path::new(bam_file).metadata();
                            let bai_index_meta_res = bai_index_path.metadata();
                            let csi_index_meta_res = csi_index_path.metadata();
                            if bam_file_meta_res.is_ok() && bai_index_meta_res.is_ok() {
                                if bam_file_meta_res.unwrap().created().unwrap()
                                    < bai_index_meta_res.unwrap().created().unwrap()
                                {
                                    return InputType::IndexedBam;
                                } else {
                                    return InputType::SortedUnindexedBam;
                                }
                            } else if bam_file_meta_res.is_ok() && csi_index_meta_res.is_ok() {
                                if bam_file_meta_res.unwrap().created().unwrap()
                                    < csi_index_meta_res.unwrap().created().unwrap()
                                {
                                    return InputType::IndexedBam;
                                } else {
                                    return InputType::SortedUnindexedBam;
                                }
                            } else {
                                return InputType::SortedUnindexedBam;
                            }
                        } else {
                            return InputType::SortedUnindexedBam;
                        }
                    } else {
                        return InputType::UnsortedBam;
                    }
                }
            }
        }
    }
}

pub fn check_input_type<P: AsRef<Path> + ToString>(p: Option<P>, bam: bool) -> InputType {
    if p.is_none() {
        if bam {
            InputType::OneBamOrSamFromStdin
        } else {
            InputType::FastqFromStdin
        }
    } else {
        let input_fn = p.unwrap().to_string();
        let input_path = std::path::PathBuf::from(&input_fn);
        if !input_path.exists() {
            quit_with_error(&format!("No such file or directory: {}", &input_fn))
        }
        if input_path.is_file() {
            if input_fn.ends_with(".fastq") || input_fn.ends_with(".fq") {
                InputType::OneFastqFile
            } else if input_fn.ends_with(".fastq.gz") || input_fn.ends_with(".fq.gz") {
                InputType::OneFastqGzippedFile
            } else if input_fn.ends_with(".sam") {
                InputType::OneSamFile
            } else if input_fn.ends_with(".bam") {
                check_bam_type(&input_fn)
            } else {
                InputType::WrongType("Bad file extension name".to_string())
            }
        } else if input_path.is_dir() {
            let mut count = 0;
            let read_dir_res = input_path.read_dir();
            if read_dir_res.is_err() {
                quit_with_error(&format!("Open directory: {} failed", input_fn))
            }
            for entry in read_dir_res.unwrap() {
                if let Ok(entry) = entry {
                    let p = entry.path();
                    let p = p.to_str().unwrap();
                    if p.ends_with(".fastq")
                        || p.ends_with(".fq")
                        || p.ends_with(".fastq.gz")
                        || p.ends_with(".fa.gz")
                    {
                        count += 1;
                    }
                }
            }
            if count == 0 {
                quit_with_error(&format!(
                    "No fastq or fastq.gz file found in directory: {}",
                    input_fn
                ))
            } else {
                InputType::DirectoryContainFastqsOrFastqsGzipped
            }
        } else {
            InputType::WrongType("Unknow input type".to_string())
        }
    }
}
