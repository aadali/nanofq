use crate::bam::BasicBamStatistics;
use crate::fastq::RecordEachStats;
use crate::utils::format_counts;
use ndhistogram::axis::Uniform;
use ndhistogram::{Histogram, ndhistogram};
use plotly::color::NamedColor;
use plotly::common::Marker;
use plotly::layout::themes::BuiltinTheme;
use plotly::layout::{Axis, GridPattern, HoverMode, LayoutGrid};
use plotly::{Bar,  Layout, Plot};
use rayon::prelude::*;
use statrs::statistics::{Data, Distribution, Max, Median, Min, OrderStatistics, Statistics};
use std::cmp::Reverse;
use std::io::Write;

struct SubReadsInfo {
    class: String,
    reads_count: usize,
    reads_percent: f64,
    bases_count: usize,
    bases_percent: f64,
}

impl SubReadsInfo {
    fn new(
        class: String,
        reads_count: usize,
        reads_percent: f64,
        bases_count: usize,
        bases_percent: f64,
    ) -> Self {
        SubReadsInfo {
            class,
            reads_count,
            reads_percent,
            bases_count,
            bases_percent,
        }
    }
}

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

        contents.push_str(&format!("ReadMinQual:\t{:.2}\n", self.min_qual));
        contents.push_str(&format!("ReadMaxQual:\t{:.2}\n", self.max_qual));
        contents.push_str(&format!("ReadMeanQual:\t{:.2}\n", self.mean_qual));
        contents.push_str(&format!("ReadStdQual:\t{:.2}\n", self.std_qual));
        contents.push_str(&format!("ReadQualQuan25:\t{:.2}\n", self.quantile25_qual));
        contents.push_str(&format!("ReadMedianQual:\t{:.2}\n", self.median_qual));
        contents.push_str(&format!("ReadQualQuan75:\t{:.2}\n", self.quantile75_qual));
        contents
    }
    fn dict_basic_info(&self) -> Vec<(&'static str, String)> {
        let mut items = vec![];
        items.push(("ReadsNumber", self.reads_number.to_string()));
        if self.bases_number / 1_000_000_000 > 1 {
            items.push((
                "BasesNumber",
                format!("{:.9}Gb", self.bases_number as f64 / 1_000_000_000.0),
            ))
        } else {
            items.push((
                "BasesNumber",
                format!("{:.6}Mb", self.bases_number as f64 / 1_000_000.0),
            ));
        }
        items.push(("N10", self.n10.to_string()));
        items.push(("N50", self.n50.to_string()));
        items.push(("N90", self.n90.to_string()));
        items.push(("MinReadLen", self.min_len.to_string()));
        items.push(("MaxReadLen", self.max_len.to_string()));
        items.push(("MeanReadLen", format!("{:.2}", self.mean_len)));
        items.push(("StdReadLen", format!("{:.2}", self.std_len)));
        items.push(("LowQuarReadLen", format!("{:.2}", self.quantile25_len)));
        items.push(("MediaReadLen", format!("{:.2}", self.median_len)));
        items.push(("UpQuarReadLen", format!("{:.2}", self.quantile75_len)));
        items.push(("MinReadQual", format!("{:.2}", self.min_qual)));
        items.push(("MaxReadQual", format!("{:.2}", self.max_qual)));
        items.push(("MeanReadQual", format!("{:.2}", self.mean_qual)));
        items.push(("StdReadQual", format!("{:.2}", self.std_qual)));
        items.push(("LowQuarReadQual", format!("{:.2}", self.quantile25_qual)));
        items.push(("MedianQuarReadQual", format!("{:.2}", self.median_qual)));
        items.push(("UpQuarReadQual", format!("{:.2}", self.quantile75_qual)));
        items
    }
}

struct PlotInfo {
    positions: Vec<f64>,
    height: Vec<usize>,
    hover_template_array: Vec<String>,
}

pub struct SummaryStats<'a> {
    all_stats: Vec<RecordEachStats>,
    read_lengths: Option<&'a [u32]>,
    read_qualities: &'a [f64],
    use_gc: bool,
    n: usize,
}
impl<'a> SummaryStats<'a> {
    pub fn new(
        all_stats: Vec<RecordEachStats>,
        read_lengths: Option<&'a [u32]>,
        read_qualities: &'a [f64],
        use_gc: bool,
        n: usize,
    ) -> Self {
        SummaryStats {
            all_stats,
            read_lengths,
            read_qualities,
            use_gc,
            n,
        }
    }

    fn get_reads_and_bases(&self) -> (usize, usize) {
        let total_bases = self
            .all_stats
            .iter()
            .map(|x| x.length as usize)
            .sum::<usize>();
        let total_reads = self.all_stats.len();
        (total_reads, total_bases)
    }

    fn get_n10_n50_n90(read_lengths: &Vec<f64>, total_length: f64) -> (u32, u32, u32) {
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
        let mut read_lengths_iter = read_lengths.iter();
        while let Some(each_len) = read_lengths_iter.next() {
            if find_n10 && find_n50 && find_n90 {
                break;
            }
            current_total_length += *each_len;
            if !find_n10 && current_total_length > n10 {
                find_n10 = true;
                n10_length = *each_len as u32;
            }
            if !find_n50 && current_total_length > n50 {
                find_n50 = true;
                n50_length = *each_len as u32;
            }
            if !find_n90 && current_total_length > n90 {
                find_n90 = true;
                n90_length = *each_len as u32;
            }
        }
        (n10_length, n50_length, n90_length)
    }
    fn get_basic_stats(&self, total_reads: usize, total_length: usize) -> BasicStatistics {
        let mut all_lengths = self
            .all_stats
            .iter()
            .map(|x| x.length as f64)
            .collect::<Vec<_>>();
        // let total_length = all_lengths.iter().sum::<f64>();
        all_lengths.par_sort_by_key(|x| Reverse((*x * 1000.0) as usize));
        let (n10, n50, n90) = Self::get_n10_n50_n90(&all_lengths, total_length as f64);
        let mut lengths = Data::new(all_lengths);
        // let total_reads = lengths.len();
        let min_len = lengths.min() as u32;
        let max_len = lengths.max() as u32;
        let mean_len = lengths.mean().unwrap();
        let std_len = lengths.iter().population_std_dev();
        let len_quantile_25 = lengths.quantile(0.25);
        let len_quantile_median = lengths.median();
        let len_quantile_75 = lengths.quantile(0.75);

        let mut qualities = Data::new(
            self.all_stats
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
        BasicStatistics {
            reads_number: total_reads,
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
        }
    }

    fn get_quality_sub_reads_info(
        &mut self,
        total_reads: usize,
        total_bases: usize,
    ) -> (Vec<SubReadsInfo>, Vec<RecordEachStats>) {
        self.all_stats
            .par_sort_unstable_by(|x, y| y.qual.partial_cmp(&x.qual).unwrap());
        let mut topn_records = vec![];
        for i in 0usize..(*[self.n, self.all_stats.len()].iter().min().unwrap()) {
            topn_records.push(self.all_stats.get(i).unwrap().clone())
        }
        let sub_reads_infos = self
            .read_qualities
            .into_par_iter()
            .map(|each_qual| {
                let better_reads = self
                    .all_stats
                    .iter()
                    .take_while(|each_stats| {
                        (each_stats.qual * 1000.0) as u32 > ((*each_qual * 1000.0) as u32)
                    })
                    .map(|each_stats| each_stats.length as usize)
                    .collect::<Vec<_>>();
                let better_reads_number = better_reads.len();
                let better_bases_number = better_reads.iter().sum::<usize>();
                SubReadsInfo::new(
                    format!("ReadQuality >= {each_qual}"),
                    better_reads_number,
                    better_reads_number as f64 / (total_reads as f64),
                    better_bases_number,
                    better_bases_number as f64 / (total_bases as f64),
                )
            })
            .collect::<Vec<_>>();
        (sub_reads_infos, topn_records)
    }

    fn get_length_sub_reads_info(
        &mut self,
        total_reads: usize,
        total_bases: usize,
    ) -> (Option<Vec<SubReadsInfo>>, Vec<RecordEachStats>) {
        self.all_stats
            .par_sort_by_key(|x| (Reverse(x.length), Reverse((x.qual * 1000.0) as u32)));
        let mut topn_records = vec![];
        for i in 0usize..(*[self.n, self.all_stats.len()].iter().min().unwrap()) {
            topn_records.push(self.all_stats.get(i).unwrap().clone())
        }

        if self.read_lengths.is_some() {
            let sub_reads_infos = self
                .read_lengths
                .unwrap()
                .into_par_iter()
                .map(|each_len| {
                    let longer_reads = self
                        .all_stats
                        .iter()
                        .take_while(|each_stats| each_stats.length >= *each_len)
                        .map(|each_stats| each_stats.length as usize)
                        .collect::<Vec<_>>();
                    let longer_reads_number = longer_reads.len();
                    let longer_base_number = longer_reads.iter().sum::<usize>();
                    SubReadsInfo::new(
                        format!("ReadLength >= {each_len}"),
                        longer_reads_number,
                        longer_reads_number as f64 / (total_reads as f64),
                        longer_base_number,
                        longer_base_number as f64 / (total_bases as f64),
                    )
                })
                .collect::<Vec<_>>();
            (Some(sub_reads_infos), topn_records)
        } else {
            (None, topn_records)
        }
    }

    fn histogram(&self, bins: usize, length_quantile: f64) -> Vec<PlotInfo> {
        let mut lengths = Data::new(
            self.all_stats
                .iter()
                .map(|x| x.length as f64)
                .collect::<Vec<_>>(),
        );
        let l = lengths.quantile(1.0 - length_quantile) as u32;
        let filtered_stats = self
            .all_stats
            .iter()
            .filter_map(|x| if x.length <= l { Some(x) } else { None })
            .collect::<Vec<_>>();

        let min_len = filtered_stats
            .iter()
            .min_by_key(|x| x.length)
            .unwrap()
            .length;

        let max_len = filtered_stats
            .iter()
            .max_by_key(|x| x.length)
            .unwrap()
            .length;

        let min_qual = self
            .all_stats
            .iter()
            .min_by_key(|x| (x.qual * 10000.0) as usize)
            .unwrap()
            .qual;
        let max_qual = self
            .all_stats
            .iter()
            .max_by_key(|x| (x.qual * 10000.0) as usize)
            .unwrap()
            .qual
            + 0.5;
        let mut len_hist =
            ndhistogram!(Uniform::new(bins, min_len as f64, max_len as f64).unwrap());
        let mut weighted_len_hist = len_hist.clone();
        let mut qual_hist =
            ndhistogram!(Uniform::new(bins, min_qual as f64, max_qual as f64).unwrap());
        let mut weighted_qual_hist = qual_hist.clone();
        for each_stats in filtered_stats {
            len_hist.fill(&(each_stats.length as f64));
            weighted_len_hist.fill_with(&(each_stats.length as f64), each_stats.length as f64);
        }
        for each_stats in self.all_stats.iter() {
            qual_hist.fill(&(each_stats.qual as f64));
            weighted_qual_hist.fill_with(&(each_stats.qual as f64), each_stats.length as f64);
        }
        let mut plot_infos = vec![];
        let mut lengths_plot_infos = vec![];
        for (idx, hist) in [len_hist, weighted_len_hist].iter().enumerate() {
            let range = "LengthRange";
            let (count, accum_count) = if idx == 0 {
                ("ReadsCount", "AccumReadsCount")
            } else {
                ("BasesCount", "AccumBasesCount")
            };
            let mut bar_positions = vec![];
            let mut bar_height = vec![];
            let mut current_accum = 0usize;
            let mut hover_template_array = vec![];
            for item in hist.iter() {
                if item.bin.start().is_none() || item.bin.end().is_none() {
                    continue;
                }
                let value = *item.value as usize;
                current_accum += value;
                let start = item.bin.start().unwrap();
                let end = item.bin.end().unwrap();
                bar_positions.push((start + end) / 2.0);
                bar_height.push(value);
                let (start_str, end_str) = (
                    (start.floor() as usize).to_string(),
                    (end.ceil() as usize).to_string(),
                );
                hover_template_array.push(format!(
                    "<b>{range}</b>: [{start_str}, {end_str})<br>\
                    <b>{count}</b>: {}<br>\
                    <b>{accum_count}</b>: {}\
                    <extra></extra>",
                    format_counts(value),
                    format_counts(current_accum)
                ));
            }
            let plot_info = PlotInfo {
                positions: bar_positions,
                height: bar_height,
                hover_template_array,
            };
            lengths_plot_infos.push(plot_info);
        }

        let mut quals_plot_infos = vec![];
        for (idx, hist) in [qual_hist, weighted_qual_hist].iter().enumerate() {
            let range = "QualRange";
            let (count, accum_count) = if idx == 0 {
                ("ReadsCount", "AccumReadsCount")
            } else {
                ("BasesCount", "AccumBasesCount")
            };
            let mut bar_positions = vec![];
            let mut bar_height = vec![];
            let mut current_accum = 0usize;
            let mut hover_template_array = vec![];
            for item in hist.iter() {
                if item.bin.start().is_none() || item.bin.end().is_none() {
                    continue;
                }
                let value = *item.value as usize;
                current_accum += value;
                let start = item.bin.start().unwrap();
                let end = item.bin.end().unwrap();
                bar_positions.push((start + end) / 2.0);
                bar_height.push(value);
                let (start_str, end_str) = (
                    format!("{}", (start * 100.0).trunc() / 100.0),
                    format!("{}", (end * 100.0).trunc() / 100.0),
                );
                hover_template_array.push(format!(
                    "<b>{range}</b>: [{start_str}, {end_str})<br>\
                    <b>{count}</b>: {}<br>\
                    <b>{accum_count}</b>: {}\
                    <extra></extra>",
                    format_counts(value),
                    format_counts(current_accum)
                ));
            }
            let plot_info = PlotInfo {
                positions: bar_positions,
                height: bar_height,
                hover_template_array,
            };
            quals_plot_infos.push(plot_info);
        }
        for (p1, p2) in lengths_plot_infos
            .into_iter()
            .zip(quals_plot_infos.into_iter())
        {
            plot_infos.push(p1);
            plot_infos.push(p2);
        }
        plot_infos
    }

    fn plot_html_div(&self, bins: usize, length_quantile: f64) -> String {
        let plot_infos = self.histogram(bins, length_quantile);
        let mut plot = Plot::new();
        for (idx, plot_info) in plot_infos.into_iter().enumerate() {
            let trace = Bar::new(plot_info.positions, plot_info.height)
                .hover_template_array(plot_info.hover_template_array)
                .x_axis(format!("x{}", idx + 1))
                .y_axis(format!("y{}", idx + 1))
                .show_legend(false)
                .marker(Marker::new().color(NamedColor::RoyalBlue));
            plot.add_trace(trace);
        }
        let layout = Layout::new()
            .hover_mode(HoverMode::XUnified)
            .grid(
                LayoutGrid::new()
                    .rows(2)
                    .columns(2)
                    .x_gap(0.1)
                    .y_gap(0.1)
                    .pattern(GridPattern::Independent),
            )
            .template(BuiltinTheme::PlotlyWhite.build())
            .x_axis(
                Axis::new()
                    // .line_color("black")
                    .show_grid(false), // .line_width(1),
            )
            .y_axis(
                Axis::new()
                    .title("ReadsCount")
                    .line_color("black")
                    .line_width(1),
            )
            .y_axis2(Axis::new().line_color("black").line_width(1))
            .x_axis3(
                Axis::new()
                    .title("ReadLength")
                    .line_color("black")
                    .show_grid(false)
                    .line_width(1),
            )
            .y_axis3(
                Axis::new()
                    .title("BasesCount")
                    .line_color("black")
                    .line_width(1),
            )
            .y_axis4(Axis::new().line_color("black").line_width(1))
            .x_axis4(
                Axis::new()
                    .title("ReadQuality")
                    .line_color("black")
                    .show_grid(false)
                    .line_width(1),
            )
            .height(800);
        plot.set_layout(layout);
        plot.to_html()
    }

    pub fn save_all_stats(&self, name: &str, out_file: &str) {
        let output_file = std::fs::File::create(out_file)
            .expect(&format!("Failed to create and open {out_file}"));
        let mut writer = std::io::BufWriter::new(output_file);
        let _ = writeln!(
            &mut writer,
            "#{name} stats result generated on {}",
            chrono::Local::now().format("%Y-%m-%d %H:%M").to_string()
        );
        for x in &self.all_stats {
            writeln!(
                &mut writer,
                "{}\t{}\t{:.6}{}",
                x.name,
                x.length,
                x.qual,
                if self.use_gc {
                    format!("\t{:.2}", x.gc.unwrap())
                } else {
                    "".to_string()
                }
            )
            .unwrap()
        }
    }

    pub fn write_to_html_file(
        &mut self,
        name: &str,
        bins: usize,
        length_quantile: f64,
        html_file: &str,
    ) {
        let home_page = "https://github.com/aadali/nanofq";
        let mut html = format!(
            r#"
        <!DOCTYPE html>
<html lang="zh">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{name} Stats Summary Report</title>
    <style>
        body {{
            font-family: "Segoe UI", "Microsoft YaHei", Arial, sans-serif;
            margin: 0;
            padding: 0;
            background-color: #f4f6f9;
        }}
        .header {{
            background: linear-gradient(135deg, #1a5276, #2E86AB);
            color: white;
            padding: 25px 40px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        .header-content {{
            max-width: 1400px;
            margin: 0 auto;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .header h1 {{ margin: 0; font-size: 28px; }}
        .header .meta {{ font-size: 14px; opacity: 0.9; }}

        .container {{
            max-width: 1400px;
            margin: 30px auto;
            padding: 0 20px;
        }}

        .section {{
            background: white;
            margin-bottom: 30px;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.08);
            overflow: hidden;
        }}
        .section-header {{
            background: #2E86AB;
            color: white;
            padding: 12px 25px;
            font-size: 18px;
            font-weight: 600;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        th, td {{
            padding: 12px 20px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 15px;
            padding: 25px;
        }}
        .stat-card {{
            background: #f8f9fa;
            padding: 15px 20px;
            border-radius: 8px;
            border-left: 5px solid #2E86AB;
        }}
        .stat-card .label {{ font-size: 13px; color: #666; }}
        .stat-card .value {{ font-size: 20px; font-weight: 600; color: #1a5276; margin-top: 5px; }}

        .chart-section {{
            padding: 25px;
        }}

        .plotly-graph-div {{
            width: 100% !important;
        }}
    </style>
</head>
<body>
    <!-- Header -->
    <div class="header">
        <div class="header-content">
            <div>
                <h1>🧬 {name} Stats Sequencing Report</h1>
                <div class="meta">Generated on {} • by <a href="{home_page}" target="_blank" style="color:#ffffff">nanofq</a></div>
            </div>
        </div>
    </div>

    <div class="container">

        <!-- Basic Statistics -->
        <div class="section">
            <div class="section-header">📊 Basic Statistics</div>
            <div class="stats-grid">
        "#,
            chrono::Local::now().format("%Y-%m-%d %H:%M").to_string()
        );

        let (total_reads, total_bases) = self.get_reads_and_bases();
        let basic_stats = self.get_basic_stats(total_reads, total_bases);
        for (key, value) in basic_stats.dict_basic_info() {
            html.push_str(&format!(
                r#"
                <div class="stat-card">
                    <div class="label">{key}</div>
                    <div class="value">{value}</div>
                </div>
            "#
            ))
        }
        html.push_str("</div>\n</div>\n");

        let (sub_reads_infos_opt, topn_length) =
            self.get_length_sub_reads_info(total_reads, total_bases);
        if self.read_lengths.is_some() {
            html.push_str(r#"
        <div class="section">
            <div class="section-header">📈 Read Length Distribution</div>
            <table>
                <thead>
                    <tr><th>Length Threshold</th><th>Reads (Percent)</th><th>Bases (Percent)</th></tr>
                </thead>
                <tbody>
            "#);
            for sub_reads_info in sub_reads_infos_opt.unwrap() {
                html.push_str(&format!(
                    r#"
                    <tr>
                        <td>{}</td>
                        <td>{}</td>
                        <td>{}</td>
                    </tr>

                "#,
                    sub_reads_info.class,
                    format!(
                        "{}({:.2}%)",
                        sub_reads_info.reads_count,
                        sub_reads_info.reads_percent * 100.0
                    ),
                    format!(
                        "{}({:.2}%)",
                        sub_reads_info.bases_count,
                        sub_reads_info.bases_percent * 100.0
                    ),
                ));
            }
            html.push_str(
                r#"
                </tbody>
                </table>
            </div>
            "#,
            );
        }
        let (sub_reads_infos, topn_quals) =
            self.get_quality_sub_reads_info(total_reads, total_bases);
        html.push_str(r#"
        <div class="section">
            <div class="section-header">📈 Read Quality Distribution</div>
            <table>
                <thead>
                    <tr><th>Quality Threshold</th><th>Reads (Percent)</th><th>Bases (Percent)</th></tr>
                </thead>
                <tbody>
        "#);
        for sub_reads_info in sub_reads_infos {
            html.push_str(&format!(
                r#"
                    <tr>
                        <td>{}</td>
                        <td>{}</td>
                        <td>{}</td>
                    </tr>

                "#,
                sub_reads_info.class,
                format!(
                    "{}({:.2}%)",
                    sub_reads_info.reads_count,
                    sub_reads_info.reads_percent * 100.0
                ),
                format!(
                    "{}({:.2}%)",
                    sub_reads_info.bases_count,
                    sub_reads_info.bases_percent * 100.0
                ),
            ));
        }
        html.push_str(
            r#"
                </tbody>
            </table>
        </div>

        "#,
        );
        html.push_str(&format!(
            r#"
        <div class="section">
            <div class="section-header">⭐ Top {} Longest Reads</div>
            <table>
                <thead>
                    <tr><th>#</th><th>Read Name</th><th>Length (bp)</th><th>Quality</th>{}</tr>
                </thead>
                <tbody>
        "#,
            self.n,
            if self.use_gc {
                "<th>GCContent</th>"
            } else {
                ""
            }
        ));
        for (idx, each_stats) in topn_length.iter().enumerate() {
            html.push_str(&format!(
                r#"
                                <tr>
                        <td>{}</td>
                        <td style="font-family: monospace; font-size: 12px;">{}</td>
                        <td>{}</td>
                        <td>{:.2}</td>
                        {}
                    </tr>
            "#,
                idx + 1,
                each_stats.name,
                each_stats.length,
                each_stats.qual,
                if self.use_gc {
                    format!("<td>{:.2}</td>", each_stats.gc.unwrap())
                } else {
                    "".to_string()
                }
            ))
        }

        html.push_str(&format!(
            r#"
                        </tbody>
            </table>
        </div>

        <div class="section">
            <div class="section-header">⭐ Top {} Highest Quality Reads</div>
            <table>
                <thead>
                    <tr><th>#</th><th>Read Name</th><th>Length (bp)</th><th>Quality</th>{}</tr>
                </thead>
                <tbody>
        "#,
            self.n,
            if self.use_gc {
                "<th>GCContent</th>"
            } else {
                ""
            }
        ));
        for (idx, each_stats) in topn_quals.iter().enumerate() {
            html.push_str(&format!(
                r#"
                                <tr>
                        <td>{}</td>
                        <td style="font-family: monospace; font-size: 12px;">{}</td>
                        <td>{}</td>
                        <td>{:.2}</td>
                        {}
                    </tr>
            "#,
                idx + 1,
                each_stats.name,
                each_stats.length,
                each_stats.qual,
                if self.use_gc {
                    format!("<td>{:.2}</td>", each_stats.gc.unwrap())
                } else {
                    "".to_string()
                }
            ))
        }
        html.push_str(
            r#"
                        </tbody>
            </table>
        </div>

        <div class="section">
            <div class="section-header">📊 Read Length And Quality Distribution</div>
            <div class="chart-section">
        "#,
        );
        let fig_html_string = self.plot_html_div(bins, length_quantile);
        html.push_str(&format!(
            r#"
                <div style="max-width: 1400px; margin: 0 auto;">
                {fig_html_string}
                </div>
            "#
        ));
        html.push_str(
            r#"
                    </div>
        </div>
    </div>
</body>
</html>

        "#,
        );
        std::fs::write(html_file, html).expect(&format!(
            "Failed to write contents into report file: {html_file}"
        ));
    }

    pub fn write_summary_to_text(
        &mut self,
        name: &str,
        basic_bam_stats: &BasicBamStatistics,
        summary_file: &str,
    ) {
        let (total_reads, total_bases) = self.get_reads_and_bases();
        let basic_stats = self.get_basic_stats(total_reads, total_bases);
        let mut contents = format!("AnalysisName:\t{name}\n");
        contents.push_str(&basic_stats.basic_info());
        let (sub_lengths_reads_infos_opt, topn_length) =
            self.get_length_sub_reads_info(total_reads, total_bases);
        let (sub_quals_reads_infos, topn_qual) =
            self.get_quality_sub_reads_info(total_reads, total_bases);
        if self.read_lengths.is_some() {
            contents.push_str(
                "#ReadLength > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
            );
            for sub_reads_info in sub_lengths_reads_infos_opt.unwrap() {
                contents.push_str(&format!(
                    "{}\t{}({:.2}%); {:.6}Mb({:.2}%)\n",
                    sub_reads_info.class,
                    sub_reads_info.reads_count,
                    sub_reads_info.reads_percent,
                    sub_reads_info.bases_count as f64 / 1_000_000.0,
                    sub_reads_info.bases_percent,
                ))
            }
        }
        contents.push_str(
            "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
        );
        for sub_reads_info in sub_quals_reads_infos {
            contents.push_str(&format!(
                "{}\t{}({:.2}%); {:.6}Mb({:.2}%)\n",
                sub_reads_info.class,
                sub_reads_info.reads_count,
                sub_reads_info.reads_percent,
                sub_reads_info.bases_count as f64 / 1_000_000.0,
                sub_reads_info.bases_percent,
            ))
        }
        contents.push_str(&format!(
            "#Top {} longest reads\nnth\tReadName\tReadLen\tReadQuality{}\n",
            self.n,
            if self.use_gc { "\tGCContent" } else { "" }
        ));

        for (idx, each_stats) in topn_length.iter().enumerate() {
            contents.push_str(&format!(
                "{}\t{}\t{}\t{:.2}{}\n",
                idx,
                each_stats.name,
                each_stats.length,
                each_stats.qual,
                if self.use_gc {
                    format!("\t{:.2}", each_stats.gc.unwrap())
                } else {
                    "".to_string()
                }
            ))
        }
        contents.push_str(&format!(
            "#Top {} highest quality reads\nnth\tReadName\tReadLen\tReadQuality{}\n",
            self.n,
            if self.use_gc { "\tGCContent" } else { "" }
        ));

        for (idx, each_stats) in topn_qual.iter().enumerate() {
            contents.push_str(&format!(
                "{}\t{}\t{}\t{:.2}{}\n",
                idx,
                each_stats.name,
                each_stats.length,
                each_stats.qual,
                if self.use_gc {
                    format!("\t{:.2}", each_stats.gc.unwrap())
                } else {
                    "".to_string()
                }
            ))
        }
        if !basic_bam_stats.is_empty() {
            contents.push_str(&basic_bam_stats.to_string());
        }
        std::fs::write(summary_file, &contents).expect(&format!(
            "Failed to write summary info into {summary_file}. The summary info is:\n{contents}"
        ))
    }
}

// pub fn get_summary(
//     all_stats: Vec<RecordEachStats>,
//     read_lengths: Option<&[u32]>,
//     read_qualities: &[f64],
//     use_gc: bool,
//     n: usize,
//     basic_bam_statistics: &BasicBamStatistics,
// ) -> (String, BasicStatistics) {
//     let basic: BasicStatistics;
//
//     // get topn length reads info
//     let mut topn_length_contents = String::from(&format!(
//         "#Top {n} longest reads\nnth\tReadName\tReadLen\tReadQuality{}\n",
//         if use_gc { "\tGCContent" } else { "" }
//     ));
//     let mut all_stats = all_stats;
//     all_stats.par_sort_by_key(|x| (Reverse(x.length), Reverse((x.qual * 1000.0) as u32)));
//     for i in 0usize..(*[n as usize, all_stats.len()].iter().min().unwrap()) {
//         let this_stats = &all_stats[i];
//         topn_length_contents.push_str(&format!(
//             "{}\t{}\t{}\t{:.2}{}\n",
//             i,
//             this_stats.name,
//             this_stats.length,
//             this_stats.qual,
//             if use_gc {
//                 format!("\t{:.2}", this_stats.gc.unwrap())
//             } else {
//                 "".to_string()
//             }
//         ))
//     }
//
//     // get lengths basic stats info
//     let mut lengths = Data::new(
//         all_stats
//             .iter()
//             .map(|x| x.length as f64)
//             .collect::<Vec<_>>(),
//     );
//     let reads_number = all_stats.len();
//     let total_length = lengths.iter().sum::<f64>();
//     let (n10, n50, n90) = get_n10_n50_n90(&all_stats, total_length);
//     let min_len = lengths.min() as u32;
//     let max_len = lengths.max() as u32;
//     let mean_len = lengths.mean().unwrap();
//     let std_len = lengths.iter().population_std_dev();
//     let len_quantile_25 = lengths.quantile(0.25);
//     let len_quantile_median = lengths.median();
//     let len_quantile_75 = lengths.quantile(0.75);
//
//     // get sub reads info depending on read length
//     let mut sub_reads_info = String::new();
//     if read_lengths.is_some() {
//         sub_reads_info.push_str(
//             "#ReadLength > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
//         );
//         let reads_infos = read_lengths
//             .unwrap()
//             .into_par_iter()
//             .map(|each_length| {
//                 let longer_reads = all_stats
//                     .iter()
//                     .take_while(|each_stats| each_stats.length > *each_length)
//                     .map(|each_stats| each_stats.length as usize)
//                     .collect::<Vec<_>>();
//                 let longer_reads_number = longer_reads.len();
//                 let longer_bases_number = longer_reads.iter().sum::<usize>();
//                 let mut longer_reads_info = String::from(&format!("ReadLength > {each_length}\t"));
//                 // let longer_reads_info =
//                 longer_reads_info.push_str(&format!(
//                     "{}({:.2}%); {:.6}Mb({:.2}%)\n",
//                     longer_reads_number,
//                     longer_bases_number as f64 / reads_number as f64 * 100.0,
//                     longer_bases_number as f64 / 1_000_000.0,
//                     longer_bases_number as f64 / total_length * 100.0
//                 ));
//                 longer_reads_info
//             })
//             .collect::<Vec<_>>();
//         for each_sub_reads_info in reads_infos {
//             sub_reads_info.push_str(&each_sub_reads_info)
//         }
//     }
//
//     // topn quality reads info
//     let mut topn_quality_contents = String::from(&format!(
//         "#Top {n} highest quality reads\nnth\tReadName\tReadLen\tReadQuality{}\n",
//         if use_gc { "\tGCContent" } else { "" }
//     ));
//     all_stats.par_sort_unstable_by(|x, y| y.qual.partial_cmp(&x.qual).unwrap());
//     for i in 0usize..(*[n as usize, all_stats.len()].iter().min().unwrap()) {
//         let this_stats = &all_stats[i];
//         topn_quality_contents.push_str(&format!(
//             "{}\t{}\t{}\t{:.2}{}\n",
//             i,
//             this_stats.name,
//             this_stats.length,
//             this_stats.qual,
//             if use_gc {
//                 format!("\t{:.2}", this_stats.gc.unwrap())
//             } else {
//                 "".to_string()
//             }
//         ))
//     }
//
//     // get sub reads info depending on read quality
//     sub_reads_info.push_str(
//         "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent); BasesNumber(BasesPercent)\n",
//     );
//     let reads_infos = read_qualities
//         .into_par_iter()
//         .map(|each_qual| {
//             let better_reads = all_stats
//                 .iter()
//                 .take_while(|each_stats| each_stats.qual as f64 > *each_qual)
//                 .map(|each_stats| each_stats.length as usize)
//                 .collect::<Vec<_>>();
//             let better_reads_number = better_reads.len();
//             let better_bases_number = better_reads.iter().sum::<usize>();
//             let mut longer_reads_info = String::from(&format!("ReadQuality > {each_qual}\t"));
//             longer_reads_info.push_str(&format!(
//                 "{}({:.2}%); {:.6}Mb({:.2}%)\n",
//                 better_reads_number,
//                 better_reads_number as f64 / reads_number as f64 * 100.0,
//                 better_bases_number as f64 / 1_000_000.0,
//                 better_bases_number as f64 / total_length as f64 * 100.0
//             ));
//             longer_reads_info
//         })
//         .collect::<Vec<_>>();
//     for each_sub_read_info in reads_infos {
//         sub_reads_info.push_str(&each_sub_read_info);
//     }
//     let mut qualities = Data::new(
//         all_stats
//             .iter()
//             .map(|x| x.qual as f64)
//             .collect::<Vec<f64>>(),
//     );
//
//     // get quality basic stats info
//     let min_qual = qualities.min();
//     let max_qual = qualities.max();
//     let mean_qual = qualities.mean().unwrap();
//     let std_qual = qualities.iter().population_std_dev();
//     let qual_quantile_25 = qualities.quantile(0.25);
//     let qual_quantile_median = qualities.median();
//     let qual_quantile_75 = qualities.quantile(0.75);
//
//     basic = BasicStatistics {
//         reads_number,
//         bases_number: total_length as usize,
//         n10,
//         n50,
//         n90,
//         min_len,
//         max_len,
//         mean_len: mean_len as f32,
//         std_len: std_len as f32,
//         quantile25_len: len_quantile_25 as f32,
//         median_len: len_quantile_median as f32,
//         quantile75_len: len_quantile_75 as f32,
//         min_qual: min_qual as f32,
//         max_qual: max_qual as f32,
//         mean_qual: mean_qual as f32,
//         std_qual: std_qual as f32,
//         quantile25_qual: qual_quantile_25 as f32,
//         median_qual: qual_quantile_median as f32,
//         quantile75_qual: qual_quantile_75 as f32,
//     };
//     let mut contents: String = String::default();
//     contents.push_str(&basic.basic_info());
//     contents.push_str(&sub_reads_info);
//     contents.push_str(&topn_length_contents);
//     contents.push_str(&topn_quality_contents);
//     if !basic_bam_statistics.is_empty() {
//         contents.push_str(&basic_bam_statistics.to_string())
//     }
//     (contents, basic)
// }
//
// pub fn write_summary(
//     all_stats: Vec<RecordEachStats>,
//     read_lengths: Option<&[u32]>,
//     read_qvalues: &[f64],
//     use_gc: bool,
//     n: usize,
//     basic_bam_stats: &BasicBamStatistics,
//     output: &str,
// ) -> BasicStatistics {
//     let (summary_info, basic_stats) = get_summary(
//         all_stats,
//         read_lengths,
//         read_qvalues,
//         use_gc,
//         n,
//         basic_bam_stats,
//     );
//     std::fs::write(output, &summary_info).expect(&format!(
//         "write summary info into {output}. The info is:\n{summary_info}"
//     ));
//     basic_stats
// }
