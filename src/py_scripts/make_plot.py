"""
This script is used for make plot for stats if --plot is set.
The content of this file will be written into /tmp/NanofqStatsPlot.py
And the default "python3" is called on the script using std::process::Command
"""
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
