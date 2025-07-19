use super::fastq::EachStats;
use rayon::prelude::*;
use std::io::Write;
use std::iter::Sum;

fn get_nx(stats_vec: &Vec<EachStats>, total_length: usize, x: f64) -> Option<usize> {
    let total_length = total_length as f64;
    let mut current_total_length = 0usize;
    for each_stats in stats_vec {
        current_total_length += each_stats.1;
        if current_total_length as f64 / total_length > x {
            return Some(each_stats.1);
        }
    }
    None
}

fn get_length_than_n(stats_vec: &Vec<EachStats>, n: usize) -> (usize, usize) {
    let mut current_total_length = 0usize;
    let mut current_reads_number = 0usize;
    for each_stats in stats_vec {
        let this_length = each_stats.1;
        if this_length > n {
            current_total_length += this_length;
            current_reads_number += 1;
        }
    }
    (current_total_length, current_reads_number)
}

fn get_quality_than_n(stats_vec: &Vec<EachStats>, n: f64) -> (usize, usize) {
    let mut current_total_length = 0usize;
    let mut current_reads_number = 0usize;
    for each_stats in stats_vec {
        let this_quality = each_stats.2;
        if this_quality.1 > n {
            current_total_length += each_stats.1;
            current_reads_number += 1;
        }
    }
    (current_total_length, current_reads_number)
}

struct LengthQuality((usize, f64));

impl Default for LengthQuality {
    fn default() -> Self {
        LengthQuality((0usize, 0.0f64))
    }
}

impl Sum for LengthQuality {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut length_quality = LengthQuality::default();
        for ele in iter {
            length_quality.0.0 += ele.0.0;
            length_quality.0.1 += ele.0.1;
        }
        length_quality
    }
}

pub fn get_summary(
    stats_vec: &mut Vec<EachStats>,
    read_lengths: Option<&Vec<usize>>,
    read_qualities: &Vec<f64>,
    n: usize,
) -> String {
    let mut contents = String::new();
    let total_reads = stats_vec.len();
    let sum_length_err_prob = stats_vec
        .into_par_iter()
        .fold(
            || LengthQuality::default(),
            |sum, element| LengthQuality((sum.0.0 + element.1, sum.0.1 + element.2.0)),
        )
        .sum::<LengthQuality>();

    let (total_bases, all_reads_avg_prob) = sum_length_err_prob.0;
    let mean_read_qual = (all_reads_avg_prob / total_reads as f64).log10() * -10.0f64;
    if total_bases / 1_000_000_000 > 1 {
        contents.push_str(&format!(
            "BaseNumber\t{:.9}Gb\n",
            total_bases as f64 / 1_000_000_000.0
        ))
    } else {
        contents.push_str(&format!(
            "BaseNumber\t{:.6}Mb\n",
            total_bases as f64 / 1_000_000.0
        ));
    }
    stats_vec.par_sort_by_key(|x| -(x.1 as isize));

    let mut top_n_lengths_reads =
        format!("#Top {n} longest reads\nnth\tReadName\tReadLen\tReadQuality\n");
    stats_vec
        .iter()
        .take(n)
        .enumerate()
        .for_each(|idx_each_stats| {
            top_n_lengths_reads.push_str(&format!(
                "{}\t{}\t{}\t{:.2}\n",
                idx_each_stats.0 + 1,
                idx_each_stats.1.0,
                idx_each_stats.1.1,
                idx_each_stats.1.2.1
            ))
        });

    contents.push_str(&format!("ReadsNumber\t{}\n", total_reads));
    contents.push_str(&format!(
        "N10\t{}\n",
        get_nx(stats_vec, total_bases, 0.1).expect("Calculate N10 Failed")
    ));
    contents.push_str(&format!(
        "N50\t{}\n",
        get_nx(stats_vec, total_bases, 0.5).expect("Calculate N50 Failed")
    ));
    contents.push_str(&format!(
        "N90\t{}\n",
        get_nx(stats_vec, total_bases, 0.9).expect("Calculate N90 Failed")
    ));
    let read_len_quantile25 = *&stats_vec[(total_reads as f64 * 0.75) as usize].1;
    let read_len_quantile50 = if total_reads % 2 == 0 {
        (&stats_vec[total_reads / 2 - 1].1 + &stats_vec[total_reads / 2 + 1].1) / 2
    } else {
        *&stats_vec[(total_reads + 1) / 2].1
    };
    let read_len_quantile75 = *&stats_vec[(total_reads as f64 * 0.25) as usize].1;
    contents.push_str(&format!("ReadLengthQuantile25\t{read_len_quantile25}\n"));
    contents.push_str(&format!("ReadLengthQuantile50\t{read_len_quantile50}\n"));
    contents.push_str(&format!("ReadLengthQuantile75\t{read_len_quantile75}\n"));
    let mean_read_length = total_bases as f64 / total_reads as f64;
    contents.push_str(&format!("ReadMenaLen\t{:.2}\n", mean_read_length));
    let len_std = (stats_vec.iter().fold(0.0, |sum, item| {
        (item.1 as f64 - mean_read_length).powf(2.0) + sum
    }) / stats_vec.len() as f64)
        .sqrt();
    contents.push_str(&format!("ReadLenStd\t{:.2}\n", len_std));

    // stats specified read length
    if let Some(read_lengths) = read_lengths {
        contents.push_str(
            "#ReadLength > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
        );
        read_lengths.iter().for_each(|each_length| {
            let (bases_number, reads_number) = get_length_than_n(stats_vec, *each_length);
            let reads_info = format!(
                "{}({:.2}%); {:.6}Mb({:.2}%)",
                reads_number,
                reads_number as f64 / total_reads as f64 * 100.0,
                bases_number as f64 / 1_000_000.0,
                bases_number as f64 / total_bases as f64 * 100.0
            );
            contents.push_str(&format!("ReadLength > {each_length}\t{reads_info}\n"))
        });
    }

    // stats_vec decreased by read quality
    stats_vec.par_sort_by(|first, second| second.2.1.partial_cmp(&first.2.1).unwrap());
    let read_qual_quantile25 = *&stats_vec[(total_reads as f64 * 0.75) as usize].2.1;
    let read_qual_quantile50 = if total_reads % 2 == 0 {
        (&stats_vec[total_reads / 2 - 1].2.1 + &stats_vec[total_reads / 2 + 1].2.1) / 2.0
    } else {
        *&stats_vec[(total_reads + 1) / 2].2.1
    };
    let read_qual_quantile75 = *&stats_vec[(total_reads as f64 * 0.25) as usize].2.1;
    contents.push_str(&format!(
        "ReadQualityQuantile25\t{:.2}\n",
        read_qual_quantile25
    ));
    contents.push_str(&format!(
        "ReadQualityQuantile50\t{:.2}\n",
        read_qual_quantile50
    ));
    contents.push_str(&format!(
        "ReadQualityQuantile75\t{:.2}\n",
        read_qual_quantile75
    ));
    contents.push_str(&format!("ReadMeanQuality\t{:.2}\n", mean_read_qual));
    contents.push_str(
        "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
    );
    read_qualities.iter().for_each(|each_qual| {
        let (bases_number, reads_number) = get_quality_than_n(stats_vec, *each_qual);
        let reads_info = format!(
            "{}({:.2}%); {:.6}Mb({:.2})%",
            reads_number,
            reads_number as f64 / total_reads as f64 * 100.0,
            bases_number as f64 / 1_000_000.0,
            bases_number as f64 / total_bases as f64 * 100.0
        );
        contents.push_str(&format!("ReadQuality > {each_qual}\t{reads_info}\n"))
    });

    contents.push_str(&top_n_lengths_reads);
    contents.push_str(&format!(
        "#Top {n} highest quality reads\nnth\tReadName\tReadLen\tReadQuality\n"
    ));
    stats_vec
        .iter()
        // .rev()
        .take(n)
        .enumerate()
        .for_each(|idx_each_stats| {
            contents.push_str(&format!(
                "{}\t{}\t{}\t{:.2}\n",
                idx_each_stats.0 + 1,
                idx_each_stats.1.0,
                idx_each_stats.1.1,
                idx_each_stats.1.2.1
            ))
        });
    contents
}

pub fn write_summary(
    stats_vec: &mut Vec<EachStats>,
    read_lengths: Option<&Vec<usize>>,
    read_qvalues: &Vec<f64>,
    n: usize,
    output: &String,
) {
    let summary_info = get_summary(stats_vec, read_lengths, read_qvalues, n);
    std::fs::write(output, &summary_info).expect(&format!(
        "write summary info into {output}. The info is:\n{summary_info}"
    ));
}
pub fn write_stats<W: Write>(
    stats_vec: &Vec<EachStats>,
    output: &mut W,
    gc: bool,
) -> Result<(), anyhow::Error> {
    let mut content = String::new();
    if gc {
        for each_stats in stats_vec {
            content.push_str(&format!(
                "{}\t{}\t{:.4}\t{:.4}\n",
                each_stats.0,
                each_stats.1,
                each_stats.2.1,
                each_stats.3.unwrap()
            ));
        }
    } else {
        stats_vec.iter().for_each(|each_stats| 
            content.push_str(
                &format!(
                    "{}\t{}\t{:.4}\n",
                    each_stats.0, each_stats.1, each_stats.2.1
                    
                )
            ));
    }
    write!(output, "{}", content)?;
    output.flush()?;
    Ok(())
}

pub fn plot() {}
