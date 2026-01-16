use super::fastq::EachStats;
use crate::bam::BasicBamStatistics;
use ansi_term;
use polars::prelude::*;
use rayon::prelude::*;
use std::cmp::max_by_key;
use std::collections::HashMap;
use std::io::Write;
use uuid;

fn stats_vec_to_dataframe(stats_vec: Vec<EachStats>) -> PolarsResult<DataFrame> {
    let lengths = stats_vec.iter().map(|x| x.1 as u32).collect::<Series>();
    let qualities = stats_vec.iter().map(|x| x.2).collect::<Series>();
    let gc = stats_vec.iter().map(|x| x.3).collect::<Series>();
    let names = stats_vec.into_iter().map(|x| *x.0).collect::<Series>();
    df!(
        "names" => names,
        "lengths" => lengths,
        "qualities" => qualities,
        "gc" => gc,
    )
}

#[derive(Default, Debug)]
pub struct BasicStatistics {
    reads_number: usize,
    bases_number: usize,
    median_qual: f32,
    mode_qual: f32,
    max_qual: f32,
    min_qual: f32,
    mean_qual: f32,
    n50: u32,
    min_len: u32,
    max_len: u32,
    mean_len: f32,
    std_len: f32,
}

fn get_nx(stats_vec: &Vec<EachStats>, total_length: usize, x: f64) -> Option<u32> {
    let total_length = total_length as f64;
    let mut current_total_length = 0usize;
    for each_stats in stats_vec {
        current_total_length += each_stats.1 as usize;
        if current_total_length as f64 / total_length > x {
            return Some(each_stats.1);
        }
    }
    None
}

fn get_nx2<I>(lengths: I, total_length: usize, x: f64) -> Option<usize>
where
    I: IntoIterator<Item = Option<u32>>,
{
    for x in lengths {
        println!("{x:?}")
    }
    Some(24)
}

fn get_length_than_n(stats_vec: &Vec<EachStats>, n: u32) -> (usize, usize) {
    let mut current_total_length = 0usize;
    let mut current_reads_number = 0usize;
    for each_stats in stats_vec {
        let this_length = each_stats.1;
        if this_length > n {
            current_total_length += this_length as usize;
            current_reads_number += 1;
        }
    }
    (current_total_length, current_reads_number)
}

fn get_quality_than_n(stats_vec: &Vec<EachStats>, n: f32) -> (usize, usize) {
    let mut current_total_length = 0usize;
    let mut current_reads_number = 0usize;
    for each_stats in stats_vec {
        let read_quality = each_stats.2;
        if read_quality > n {
            current_total_length += each_stats.1 as usize;
            current_reads_number += 1;
        }
    }
    (current_total_length, current_reads_number)
}

fn get_read_qual_mode(stats_vec: &Vec<EachStats>) -> f32 {
    let mut counter = HashMap::new();
    let epsilon = 0.01;
    for each_stats in stats_vec {
        let key = (each_stats.2 / epsilon).round() as usize;
        counter
            .entry(key)
            .and_modify(|count| *count += 1)
            .or_insert(1usize);
    }
    let res = counter.iter().fold((0usize, 0usize), |x, y| {
        max_by_key(x, (*y.0, *y.1), |a| a.1)
    });
    res.0 as f32 / 100.0
}
pub fn get_summary2(
    stats_vec: Vec<EachStats>,
    read_lengths: Option<&Vec<u32>>,
    read_qualities: &[f64],
    n: usize,
    basic_bam_statistics: &BasicBamStatistics,
) -> (String, BasicStatistics) {
    let mut contents: String = String::default();
    let stats_df =
        stats_vec_to_dataframe(stats_vec).expect("Covert stats vector ot DataFrame failed");
    let lengths_col = stats_df["lengths"].as_series().unwrap().u32().unwrap();
    get_nx2(lengths_col, 2932, 9.32);
    let stats_lazy = stats_df.lazy();
    // names, lengths, qualities, gc
    stats_lazy
        .clone()
        .select([
            col("lengths").count().alias("ReadsNumber"),
            col("lengths").sum().alias("BaseNumber"),
            col("lengths")
                .quantile(lit(0.25), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadLengthQuantile25"),
            col("lengths")
                .quantile(lit(0.5), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadLengthQuantile50"),
            col("lengths")
                .quantile(lit(0.75), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadLengthQuantile75"),
            col("lengths")
                .mean()
                .round(2, RoundMode::default())
                .alias("ReadMeanLen"),
            col("lengths")
                .std(0)
                .round(2, RoundMode::default())
                .alias("ReadLenStd"),
            col("qualities")
                .quantile(lit(0.25), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadQualityQuantile25"),
            col("qualities")
                .quantile(lit(0.5), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadQualityMedia"),
            col("qualities")
                .quantile(lit(0.75), QuantileMethod::Linear)
                .round(2, RoundMode::default())
                .alias("ReadQualityQuantile75"),
            col("qualities")
                .mean()
                .round(2, RoundMode::default())
                .alias("ReadMeanQuality"),
        ])
        .collect()
        .expect("Collect summary dataframe failed");

    let topn_qualities = stats_lazy
        .clone()
        .sort(
            ["qualities"],
            SortMultipleOptions::default().with_order_descending(true),
        )
        .collect()
        .unwrap()
        .head(Some(n))
        .insert_column(
            0,
            Series::new("nth".into(), (1..=n as u32).collect::<Vec<u32>>()),
        )
        .unwrap();

    let topn_lengths = stats_lazy
        .clone()
        .sort(
            ["lengths"],
            SortMultipleOptions::default().with_order_descending(true),
        )
        .collect()
        .unwrap()
        .head(Some(n))
        .insert_column(
            0,
            Series::new("nth".into(), (1..=n as u32).collect::<Vec<u32>>()),
        )
        .unwrap();
    (contents, BasicStatistics::default())
}
pub fn get_summary(
    stats_vec: &mut Vec<EachStats>,
    read_lengths: Option<&Vec<usize>>,
    read_qualities: &[f64],
    n: usize,
    basic_bam_stats: &BasicBamStatistics,
) -> (String, BasicStatistics) {
    let mut contents = String::new();
    let total_reads = stats_vec.len();
    let total_bases = stats_vec
        .iter()
        .map(|each_stats| each_stats.1 as usize)
        .sum::<usize>();
    let mean_read_qual =
        stats_vec.iter().map(|each_stats| each_stats.2).sum::<f32>() / stats_vec.len() as f32;
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
    let max_read_len = stats_vec[0].1;
    let min_read_len = stats_vec.iter().last().unwrap().1;

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
                idx_each_stats.1.2
            ))
        });

    contents.push_str(&format!("ReadsNumber\t{}\n", total_reads));
    contents.push_str(&format!(
        "N10\t{}\n",
        get_nx(stats_vec, total_bases, 0.1).expect("Calculate N10 Failed")
    ));
    let n50 = get_nx(stats_vec, total_bases, 0.5).expect("Calculate N50 Failed");
    contents.push_str(&format!("N50\t{}\n", n50));
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
    let mean_read_length = total_bases as f32 / total_reads as f32;
    contents.push_str(&format!("ReadMeanLen\t{:.2}\n", mean_read_length));
    let len_std = (stats_vec.iter().fold(0.0, |sum, item| {
        (item.1 as f32 - mean_read_length).powf(2.0) + sum
    }) / stats_vec.len() as f32)
        .sqrt();
    contents.push_str(&format!("ReadLenStd\t{:.2}\n", len_std));

    // stats specified read length
    if let Some(read_lengths) = read_lengths {
        contents.push_str(
            "#ReadLength > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
        );
        read_lengths.iter().for_each(|each_length| {
            let (bases_number, reads_number) = get_length_than_n(stats_vec, *each_length as u32);
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
    stats_vec.par_sort_by(|first, second| second.2.partial_cmp(&first.2).expect("NAN was found"));
    let max_read_qual = stats_vec[0].2;
    let min_read_qual = stats_vec.iter().last().unwrap().2;
    let read_qual_quantile25 = *&stats_vec[(total_reads as f64 * 0.75) as usize].2;
    let read_qual_quantile50 = if total_reads % 2 == 0 {
        (&stats_vec[total_reads / 2 - 1].2 + &stats_vec[total_reads / 2 + 1].2) / 2.0
    } else {
        *&stats_vec[(total_reads + 1) / 2].2
    };
    let read_qual_quantile75 = *&stats_vec[(total_reads as f64 * 0.25) as usize].2;
    contents.push_str(&format!(
        "ReadQualityQuantile25\t{:.2}\n",
        read_qual_quantile25
    ));
    contents.push_str(&format!("ReadQualityMedian\t{:.2}\n", read_qual_quantile50));
    contents.push_str(&format!(
        "ReadQualityQuantile75\t{:.2}\n",
        read_qual_quantile75
    ));
    contents.push_str(&format!("ReadMeanQuality\t{:.2}\n", mean_read_qual));
    contents.push_str(
        "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
    );
    read_qualities.iter().for_each(|each_qual| {
        let (bases_number, reads_number) = get_quality_than_n(stats_vec, *each_qual as f32);
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
                idx_each_stats.1.2
            ))
        });
    if !basic_bam_stats.is_empty() {
        contents.push_str(&basic_bam_stats.to_string())
    }
    (
        contents,
        BasicStatistics {
            reads_number: total_reads,
            bases_number: total_bases,
            median_qual: read_qual_quantile50,
            mode_qual: get_read_qual_mode(stats_vec),
            max_qual: max_read_qual,
            min_qual: min_read_qual,
            mean_qual: mean_read_qual,
            n50: n50,
            min_len: min_read_len,
            max_len: max_read_len,
            mean_len: mean_read_length,
            std_len: len_std,
        },
    )
}

pub fn write_summary(
    stats_vec: &mut Vec<EachStats>,
    read_lengths: Option<&Vec<usize>>,
    read_qvalues: &[f64],
    n: usize,
    basic_bam_stats: &BasicBamStatistics,
    output: &str,
) -> BasicStatistics {
    let (summary_info, basic_stats) =
        get_summary(stats_vec, read_lengths, read_qvalues, n, basic_bam_stats);
    std::fs::write(output, &summary_info).expect(&format!(
        "write summary info into {output}. The info is:\n{summary_info}"
    ));
    basic_stats
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
                each_stats.2,
                each_stats.3.unwrap()
            ));
        }
    } else {
        stats_vec.iter().for_each(|each_stats| {
            content.push_str(&format!(
                "{}\t{}\t{:.4}\n",
                each_stats.0, each_stats.1, each_stats.2
            ))
        });
    }
    write!(output, "{}", content)?;
    output.flush()?;
    Ok(())
}

pub fn make_plot(
    basic_statistics: &BasicStatistics,
    quan: f64,
    prefix: &str,
    format: &Vec<String>,
    python: &str,
    stats_file: &str,
) -> Result<(), anyhow::Error> {
    let formats = format.join(",");
    // use uuid to uniq each script of process if this subcommand called by different process at same time
    let script = &format!("/tmp/NanofqStatsPlot_{}.py", uuid::Uuid::new_v4());
    let cmd = format!(
        "PYTHON3 {} --input {} --quan {:.2} --n50 {} --len_bins 100 --qual_bins 100 --mode_qual {:.2} --prefix {} --format {}",
        script, stats_file, quan, basic_statistics.n50, basic_statistics.mode_qual, prefix, formats
    );
    std::fs::write(script, PYTHON_PLOT_SCRIPT)?;
    let cmd_result = std::process::Command::new(python)
        .arg(script)
        .arg("--input")
        .arg(stats_file)
        .arg("--quan")
        .arg(format!("{:2}", quan))
        .arg("--n50")
        .arg(format!("{}", basic_statistics.n50))
        .arg("--len_bins")
        .arg("100")
        .arg("--qual_bins")
        .arg("100")
        .arg("--mode_qual")
        .arg(format!("{:.2}", basic_statistics.mode_qual))
        .arg("--prefix")
        .arg(prefix)
        .arg("--format")
        .arg(formats)
        .output();
    match cmd_result {
        Ok(output) => {
            if !output.status.success() {
                println!("status: {}", output.status);
                eprintln!("std_err: {}", std::str::from_utf8(&output.stderr)?);
                eprintln!(
                    "{}\n{}\n{}",
                    ansi_term::Color::Yellow.paint("Stats finished but make plot failed. You can use this command to make plot:"),
                    ansi_term::Color::Green.paint(cmd),
                    ansi_term::Color::Yellow.paint("Replace PYTHON3 with your own python3 path; Matplotlib is needed")
                );
                std::process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("{}", ansi_term::Color::Red.paint(format!("{:?}", e)));
            std::process::exit(1);
        }
    }
    Ok(())
}

const PYTHON_PLOT_SCRIPT: &str = r#"
from matplotlib import pyplot as plt
import numpy as np
import argparse

ALLOWED_FORMATS = ["jpg", "pdf", "svg", "png"]


def get_arguments():
    parser = argparse.ArgumentParser("plot long reads length and quality")
    parser.add_argument("--input", help="the tsv file, col2 is length of read and col3 is quality of read",
                        required=True)
    parser.add_argument("--quan", type=float,
                        help="the shortest ratio and longest ratio of reads will not be rendered on figure, should be in range(0.0, 1.0)",
                        default=0.01)
    parser.add_argument("--n50", type=int, help="n50 of all reads", required=True)
    parser.add_argument("--mode_qual", type=float, help="the mode quality of reads", required=True)
    parser.add_argument("--len_bins", type=int, help="how many bins used for read lengths distribution", default=100)
    parser.add_argument("--qual_bins", type=int, help="how many bins used for read quality distribution", default=100)
    parser.add_argument("--prefix", type=str,
                        help="the prefix of output figures path, the last word will be used as figure name",
                        required=True)
    parser.add_argument("--format", type=str,
                        help="what format you need, can be one of jpg, pdf, svg, png or tow or more of them separated by coma",
                        default="pdf")
    return parser.parse_args()


def read_stats_tsv(file_path: str):
    lens_quals = []
    with open(file_path, 'r') as infile:
        for line in infile:
            fields = line.strip().split("\t")
            read_len = int(fields[1])
            read_qual = float(fields[2])
            lens_quals.append((read_len, read_qual))
    return lens_quals


def set_labels(ax: plt.axes, axis: str, unit: str):
    assert axis in ['x', 'y']
    assert unit in ['Mb', 'k']
    if axis == 'x':
        ticks = ax.get_xticks()
        if unit == "Mb":
            labels = [f"{i / 1000000}{unit}" for i in ticks]
        else:
            labels = [f"{i / 1000}{unit}" for i in ticks]
    else:
        ticks = ax.get_yticks()
        if unit == "Mb":
            labels = [f"{i / 1000000}{unit}" for i in ticks]
        else:
            labels = [f"{i / 1000}{unit}" for i in ticks]
    if ticks[0] < 0:
        ticks = ticks[1:-1]
        labels = labels[1:-1]

    if axis == 'x':
        ax.xaxis.set_ticks(ticks)
        ax.set_xticklabels(labels)
    else:
        ax.yaxis.set_ticks(ticks)
        ax.set_yticklabels(labels)
    return ticks, labels


def plot(lens_quals: list[tuple],
         len_quan: float,
         n50: int,
         mode_qual: float,
         len_bins: int,
         qual_bins: int,
         prefix: str,
         formats: list[str]):
    fig, axes = plt.subplots(2, 2, figsize=(20, 12),
                             # layout="constrained",
                             sharex="col")
    assert 0 <= len_quan <= 1
    ax00, ax01, ax10, ax11 = axes.flatten()
    lengths = list(map(lambda x: x[0], lens_quals))
    quals = list(map(lambda x: x[1], lens_quals))
    lower_lmt_len, upper_lmt_len = np.quantile(lengths, [len_quan, 1 - len_quan])
    lengths2 = list(filter(lambda x: lower_lmt_len < x < upper_lmt_len, lengths))
    read_count, read_len_bins = np.histogram(lengths2, len_bins)
    bases_count, _ = np.histogram(lengths2, len_bins, weights=lengths2)
    ax00.bar(read_len_bins[:-1], bases_count, width=np.diff(read_len_bins))
    ax10.bar(read_len_bins[:-1], read_count, width=np.diff(read_len_bins))

    read_count_q, read_qual_bins = np.histogram(quals, qual_bins)
    base_count_q, _ = np.histogram(quals, qual_bins, weights=lengths)
    ax01.bar(read_qual_bins[1:], base_count_q, width=np.diff(read_qual_bins))
    ax11.bar(read_qual_bins[1:], read_count_q, width=np.diff(read_qual_bins))
    for each_ax in [ax00, ax10, ax01, ax11]:
        ymin, ymax = each_ax.get_ylim()
        each_ax.set_ylim((ymin - ymax) * 0.02, ymax)

    for each_ax in [ax00, ax01]:
        set_labels(each_ax, 'y', 'Mb')

    for each_ax in [ax10, ax11]:
        set_labels(each_ax, 'y', 'k')

    for each_ax in [ax00, ax10]:
        xmin, xmax = each_ax.get_xlim()
        if xmin < n50 < xmax:
            each_ax.axvline(n50, 0, 1, color="black", linewidth=0.8, linestyle="dashed")
            if each_ax == ax00:
                ax00.annotate(f"N50={n50}", xy=(n50, 1), xycoords=("data", "axes fraction"),
                              ha='left', va="bottom", rotation=30, color="red", fontsize=10)

    for each_ax in [ax01, ax11]:
        each_ax.axvline(mode_qual, 0, 1, color="black", linewidth=1, linestyle="dashed")
        if each_ax == ax01:
            each_ax.annotate(f"ModeQual={mode_qual:.2f}", xy=(mode_qual, 1), xycoords=('data', 'axes fraction'),
                             ha='left', va='bottom', rotation=30, color='red', fontsize=10)
    set_labels(ax10, 'x', 'k')
    ax00.set_ylabel("NumberOfBases", fontsize=15)
    ax01.set_ylabel("NumberOfBases", fontsize=15)
    ax10.set_xlabel("ReadLength", fontsize=15)
    ax11.set_xlabel("ReadQuality", fontsize=15)
    ax10.set_ylabel("ReadCount", fontsize=15)
    ax11.set_ylabel("ReadCount", fontsize=15)
    name = prefix.strip().split("/")[-1]
    if len(name) == 0:
        fig_name = "ReadLengthAndQualityDistribution"
    else:
        fig_name = f"{name}:ReadLengthAndQualityDistribution"
    name = prefix.split("/")[-1]
    fig.suptitle(fig_name, fontsize=30, fontweight="bold")
    for f in formats:
        assert f in ALLOWED_FORMATS, f"{f} is not supported"
        fig.savefig(f"{prefix}__ReadLengthAndQualityDistribution.{f}")


if __name__ == '__main__':
    arguments = get_arguments()
    input_file = arguments.input
    quan = arguments.quan
    n50 = arguments.n50
    mode_qual = arguments.mode_qual
    len_bins = arguments.len_bins
    qual_bins = arguments.qual_bins
    prefix = arguments.prefix
    formats = arguments.format
    for f in formats.strip().split(","):
        if f not in ALLOWED_FORMATS:
            print(f"{f} is not supported format figure")
            exit(1)

    read_stats = read_stats_tsv(input_file)
    plot(read_stats,
         len_quan=quan,
         n50=n50,
         mode_qual=mode_qual,
         len_bins=len_bins,
         qual_bins=qual_bins,
         prefix=prefix,
         formats=formats.strip().split(","))
"#;
