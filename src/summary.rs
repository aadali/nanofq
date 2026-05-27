use crate::bam::BasicBamStatistics;
use crate::fastq2::RecordEachStats;
use ansi_term;
use rayon::prelude::*;
use statrs::statistics::{Data, Distribution, Max, Median, Min, OrderStatistics, Statistics};
use std::cmp::Reverse;
use uuid;


#[derive(Default, Debug)]
pub struct BasicStatistics {
    reads_number: usize,
    bases_number: usize,
    n10: u32,
    n50: u32,
    n90: u32,

    min_len: u32,
    max_len: u32,
    mean_len: f32,
    std_len: f32,
    quantile25_len: f32,
    median_len: f32,
    quantile75_len: f32,

    min_qual: f32,
    max_qual: f32,
    mean_qual: f32,
    std_qual: f32,
    quantile25_qual: f32,
    median_qual: f32,
    quantile75_qual: f32,
}

impl BasicStatistics {
    fn basic_info(&self) -> String {
        let mut contents = String::new();
        if self.bases_number / 1_000_000_000 > 1 {
            contents.push_str(&format!(
                "BasesNumber:\t{:.9}Gb\n",
                self.bases_number as f64 / 1_000_000_000.0
            ))
        } else {
            contents.push_str(&format!(
                "BasesNumber:\t{:.6}Mb\n",
                self.bases_number as f64 / 1_000_000.0
            ));
        }
        contents.push_str(&format!("ReadsNumber:\t{}\n", self.reads_number));

        contents.push_str(&format!("N10:\t{}\n", self.n10));
        contents.push_str(&format!("N50:\t{}\n", self.n50));
        contents.push_str(&format!("N90:\t{}\n", self.n90));

        contents.push_str(&format!("ReadMinLen:\t{}\n", self.min_len));
        contents.push_str(&format!("ReadMaxLen:\t{}\n", self.max_len));
        contents.push_str(&format!("ReadMeanLen:\t{:.2}\n", self.mean_len));
        contents.push_str(&format!("ReadStdLen:\t{:.2}\n", self.std_len));
        contents.push_str(&format!("ReadLenQuan25:\t{:.2}\n", self.quantile25_len));
        contents.push_str(&format!("ReadMedianLen:\t{:.2}\n", self.median_len));
        contents.push_str(&format!("ReadLenQuan75:\t{:.2}\n", self.quantile75_len));

        contents.push_str(&format!("ReadMinQual:\t{}\n", self.min_qual));
        contents.push_str(&format!("ReadMaxQual:\t{}\n", self.max_qual));
        contents.push_str(&format!("ReadMeanQual:\t{:.2}\n", self.mean_qual));
        contents.push_str(&format!("ReadStdQual:\t{:.2}\n", self.std_qual));
        contents.push_str(&format!("ReadQualQuan25:\t{:.2}\n", self.quantile25_qual));
        contents.push_str(&format!("ReadMedianQual:\t{:.2}\n", self.median_qual));
        contents.push_str(&format!("ReadQualQuan75:\t{:.2}\n", self.quantile75_qual));
        contents
    }
}

fn get_n10_n50_n90(all_stats: &Vec<RecordEachStats>, total_length: f64) -> (u32, u32, u32) {
    let mut current_total_length = 0f64;
    let n10 = total_length * 0.10;
    let n50 = total_length * 0.50;
    let n90 = total_length * 0.90;
    let mut find_n10 = false;
    let mut find_n50 = false;
    let mut find_n90 = false;
    let mut n10_length: u32 = 0;
    let mut n50_length: u32 = 0;
    let mut n90_length: u32 = 0;
    let mut all_stats_iter = all_stats.iter();
    while let Some(each_stats) = all_stats_iter.next() {
        if find_n10 && find_n50 && find_n90 {
            break;
        }
        current_total_length += each_stats.length as f64;
        if !find_n10 && current_total_length > n10 {
            find_n10 = true;
            n10_length = each_stats.length;
        }
        if !find_n50 && current_total_length > n50 {
            find_n50 = true;
            n50_length = each_stats.length;
        }
        if !find_n90 && current_total_length > n90 {
            find_n90 = true;
            n90_length = each_stats.length;
        }
    }
    (n10_length, n50_length, n90_length)
}

pub fn get_summary(
    all_stats: Vec<RecordEachStats>,
    read_lengths: Option<&[u32]>,
    read_qualities: &[f64],
    n: usize,
    basic_bam_statistics: &BasicBamStatistics,
) -> (String, BasicStatistics) {
    let basic: BasicStatistics;
    let mut topn_length_contents = String::from(&format!(
        "#Top {n} longest reads\nnth\tReadName\tReadLen\tReadQuality\n"
    ));
    let mut all_stats = all_stats;
    all_stats.sort_by_key(|x| Reverse(x.length));
    for i in 0usize..(*[n as usize, all_stats.len()].iter().min().unwrap()) {
        let this_stats = &all_stats[i];
        topn_length_contents.push_str(&format!(
            "{}\t{}\t{}\t{:.2}\n",
            i, this_stats.name, this_stats.length, this_stats.qual
        ))
    }
    let mut lengths = Data::new(
        all_stats
            .iter()
            .map(|x| x.length as f64)
            .collect::<Vec<_>>(),
    );
    let reads_number = all_stats.len();
    let total_length = lengths.iter().sum::<f64>();
    let (n10, n50, n90) = get_n10_n50_n90(&all_stats, total_length);
    let min_len = lengths.min() as u32;
    let max_len = lengths.max() as u32;
    let mean_len = lengths.mean().unwrap();
    let std_len = lengths.iter().population_std_dev();
    let len_quantile_25 = lengths.quantile(0.25);
    let len_quantile_median = lengths.median();
    let len_quantile_75 = lengths.quantile(0.75);

    let mut sub_reads_info = String::new();
    if read_lengths.is_some() {
        sub_reads_info.push_str(
            "#ReadLength > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
        );
        let reads_infos = read_lengths
            .unwrap()
            .into_par_iter()
            .map(|each_length| {
                let longer_reads = all_stats
                    .iter()
                    .take_while(|each_stats| each_stats.length > *each_length)
                    .map(|each_stats| each_stats.length as usize)
                    .collect::<Vec<_>>();
                let longer_reads_number = longer_reads.len();
                let longer_bases_number = longer_reads.iter().sum::<usize>();
                let mut longer_reads_info = String::from(&format!("ReadLength > {each_length}\t"));
                // let longer_reads_info =
                longer_reads_info.push_str(&format!(
                    "{}({:.2}%); {:.6}Mb({:.2}%)\n",
                    longer_reads_number,
                    longer_bases_number as f64 / reads_number as f64 * 100.0,
                    longer_bases_number as f64 / 1_000_000.0,
                    longer_bases_number as f64 / total_length as f64 * 100.0
                ));
                longer_reads_info
            })
            .collect::<Vec<_>>();
        for each_sub_reads_info in reads_infos {
            sub_reads_info.push_str(&each_sub_reads_info)
        }
    }

    let mut topn_quality_contents = String::from(&format!(
        "#Top {n} highest quality reads\nnth\tReadName\tReadLen\tReadQuality\n"
    ));
    all_stats.sort_unstable_by(|x, y| y.qual.partial_cmp(&x.qual).unwrap());
    for i in 0usize..(*[n as usize, all_stats.len()].iter().min().unwrap()) {
        let this_stats = &all_stats[i];
        topn_quality_contents.push_str(&format!(
            "{}\t{}\t{}\t{:.2}\n",
            i, this_stats.name, this_stats.length, this_stats.qual
        ))
    }
    sub_reads_info.push_str(
        "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
    );
    let reads_infos = read_qualities
        .into_par_iter()
        .map(|each_qual| {
            let better_reads = all_stats
                .iter()
                .take_while(|each_stats| each_stats.qual as f64 > *each_qual)
                .map(|each_stats| each_stats.length as usize)
                .collect::<Vec<_>>();
            let better_reads_number = better_reads.len();
            let better_bases_number = better_reads.iter().sum::<usize>();
            let mut longer_reads_info = String::from(&format!("ReadQuality > {each_qual}\t"));
            // let longer_reads_info =
            longer_reads_info.push_str(&format!(
                "{}({:.2}%); {:.6}Mb({:.2}%)\n",
                better_reads_number,
                better_reads_number as f64 / reads_number as f64 * 100.0,
                better_bases_number as f64 / 1_000_000.0,
                better_bases_number as f64 / total_length as f64 * 100.0
            ));
            longer_reads_info
        })
        .collect::<Vec<_>>();
    for each_sub_read_info in reads_infos {
        sub_reads_info.push_str(&each_sub_read_info);
    }
    let mut qualities = Data::new(
        all_stats
            .iter()
            .map(|x| x.qual as f64)
            .collect::<Vec<f64>>(),
    );
    let min_qual = qualities.min();
    let max_qual = qualities.max();
    let mean_qual = qualities.mean().unwrap();
    let std_qual = qualities.iter().population_std_dev();
    let qual_quantile_25 = qualities.quantile(0.25);
    let qual_quantile_median = qualities.median();
    let qual_quantile_75 = qualities.quantile(0.75);
    basic = BasicStatistics {
        reads_number,
        bases_number: total_length as usize,
        n10,
        n50,
        n90,
        min_len,
        max_len,
        mean_len: mean_len as f32,
        std_len: std_len as f32,
        quantile25_len: len_quantile_25 as f32,
        median_len: len_quantile_median as f32,
        quantile75_len: len_quantile_75 as f32,
        min_qual: min_qual as f32,
        max_qual: max_qual as f32,
        mean_qual: mean_qual as f32,
        std_qual: std_qual as f32,
        quantile25_qual: qual_quantile_25 as f32,
        median_qual: qual_quantile_median as f32,
        quantile75_qual: qual_quantile_75 as f32,
    };
    let mut contents: String = String::default();
    contents.push_str(&basic.basic_info());
    contents.push_str(&sub_reads_info);
    contents.push_str(&topn_length_contents);
    contents.push_str(&topn_quality_contents);
    if !basic_bam_statistics.is_empty() {
        contents.push_str(&basic_bam_statistics.to_string())
    }
    (contents, basic)
}

pub fn write_summary(
    all_stats: Vec<RecordEachStats>,
    read_lengths: Option<&[u32]>,
    read_qvalues: &[f64],
    n: usize,
    basic_bam_stats: &BasicBamStatistics,
    output: &str,
) -> BasicStatistics {
    let (summary_info, basic_stats) =
        get_summary(all_stats, read_lengths, read_qvalues, n, basic_bam_stats);
    std::fs::write(output, &summary_info).expect(&format!(
        "write summary info into {output}. The info is:\n{summary_info}"
    ));
    basic_stats
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
        "PYTHON3 {} --input {} --quan {:.2} --n50 {} --len_bins 100 --qual_bins 100 --median_qual {:.2} --prefix {} --format {}",
        script,
        stats_file,
        quan,
        basic_statistics.n50,
        basic_statistics.median_qual,
        prefix,
        formats
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
        .arg("--median_qual")
        .arg(format!("{:.2}", basic_statistics.median_qual))
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
    parser.add_argument("--median_qual", type=float, help="the median quality of reads", required=True)
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
         median_qual: float,
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
        each_ax.axvline(median_qual, 0, 1, color="black", linewidth=1, linestyle="dashed")
        if each_ax == ax01:
            each_ax.annotate(f"MedianQual={median_qual:.2f}", xy=(median_qual, 1), xycoords=('data', 'axes fraction'),
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
    median_qual = arguments.median_qual
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
         median_qual=median_qual,
         len_bins=len_bins,
         qual_bins=qual_bins,
         prefix=prefix,
         formats=formats.strip().split(","))
"#;
