import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import argparse


def get_arguments():
    parser = argparse.ArgumentParser("plot long reads length and quality")
    parser.add_argument("-i", "--input", help="a dataframe file separated by tab", required=True)
    parser.add_argument("-l", "--length_col", help="the length column index in dataframe, default 1", type=int, default=1)
    parser.add_argument("-q", "--quality_col", help="the quality column index in dataframe, default 2", type=int, default=2)
    parser.add_argument("-p", "--prefix", help="prefix of output figure name prefix", default="./plot_output")
    parser.add_argument("-f", "--format", help="what format do you need", action="append", choices=["png", "pdf", "jpg"])
    parser.add_argument("-m", "--plot_mean_length", help="whether plot the mean length", action="store_true")
    parser.add_argument("-M", "--mean_length", help="the mean of length, if not set, will calculate from --input file", type=float)
    parser.add_argument("-n", "--plot_n50", help="whether plot the N50 length", action="store_true")
    parser.add_argument("-N", "--n50", help="the N50 value, if not set, don't plot this line", type=float)
    parser.add_argument("-s", "--std", help="the standard deviation of length, if not set , will calculate from --input file. "
                                            "Values that outside the range of three std will be ignored", type=float)
    parser.add_argument("-Q", "--mean_quality", help="the mean of read Q value, if not set, don't plot this line", type=float)
    parser.add_argument("--header", help="if set, treat the first row as header", action="store_true")
    return parser.parse_args()


def plot(df: pd.DataFrame,plot_mean_length: bool,  mean_length: float, plot_n50: bool, n50: float, mean_quality: float, prefix: str, format: [str]):
    length_bins = 60
    quality_bins = 50
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout="constrained", sharex="col")
    ax00, ax01, ax10, ax11 = axes.flatten()
    ax00.hist(df['length'], bins=length_bins, weights=df['length'], rwidth=1)
    ax10.hist(df['length'], bins=length_bins, rwidth=1)
    ax01.hist(df['quality'], bins=quality_bins, weights=df['length'], rwidth=1)
    ax11.hist(df['quality'], bins=quality_bins, rwidth=1)
    for each_ax in [ax00, ax10, ax01, ax11]:
        ymin, ymax = each_ax.get_ylim()
        each_ax.set_ylim((ymin - ymax) * 0.02, ymax)

    for each_ax in [ax00, ax01]:
        labels = [f"{y / 1000000}Mb" for y in each_ax.get_yticks()]
        each_ax.set_yticklabels(labels)

    for each_ax in [ax10, ax11]:
        labels = [f"{y / 1000}k" for y in each_ax.get_yticks()]
        each_ax.set_yticklabels(labels)

    ax10.set_xticklabels([f"{x / 1000}k" for x in ax10.get_xticks()])
    ax00.set_ylabel("Number of bases")
    ax10.set_ylabel("Number of reads")
    ax10.set_xlabel("Read length")
    ax11.set_xlabel("Read quality")
    for ax in [ax00, ax10]:
        if plot_mean_length and mean_length > 0:
            ax.axvline(mean_length, 0, 1, color="black", linewidth=0.8, linestyle="dashed")  # always plot the of mean_length vline
            if ax == ax00:
                # just annotate the mean_length of ax00
                ax00.annotate(f"MeanLen={int(mean_length)}", xy=(mean_length, 1), xycoords=('data', 'axes fraction'),
                              ha='left', va='bottom', rotation=30, color='red', fontsize="xx-small")
        if plot_n50 and n50 is not None:
            ax.axvline(n50, 0, 1, color="red", linewidth=1, linestyle="dashed")  # if set the n50 value, so plot N50 vline
            if ax == ax00:
                ax00.annotate(f"N50={int(n50)}", xy=(n50, 1), xycoords=('data', 'axes fraction'), ha='left',
                              va='bottom', rotation=30, color="red", fontsize="xx-small")

    for ax in [ax01, ax11]:
        if mean_quality is not None:
            ax.axvline(mean_quality, 0, 1, color="black", linewidth=1, linestyle="dashed")  # if set the mean_quality, so plot mean quality vline
            if ax == ax01:
                ax01.annotate(f"MeanQ={round(mean_quality, 2)}", xy=(mean_quality, 1), xycoords=('data', 'axes fraction'),
                              ha='left', va='bottom', rotation=30, color="red", fontsize="xx-small")


    fig.suptitle("Read length and quality distribution")
    for f in format:
        fig.savefig(f"{prefix}__ReadLengthAndQualityDistribution.{f}")

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1.5), height_ratios=(1.5, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    ax = fig.add_subplot(gs[1, 0])
    ax_length = fig.add_subplot(gs[0, 0])
    ax_quality = fig.add_subplot(gs[1, 1])
    ax.scatter(df['length'], df['quality'], s=1)
    ax.set_xticklabels([f"{x / 1000}k" for x in ax.get_xticks()])
    ax.set_ylim(5, 30)
    ax_quality.set_ylim(5, 30)
    ax_length.hist(df['length'], bins=length_bins, rwidth=1)
    ymin, ymax = ax_length.get_ylim()
    ax_length.set_ylim((ymin - ymax) * 0.02, ymax)
    ax_quality.hist(df['quality'], bins=quality_bins, rwidth=1, orientation='horizontal')
    ax_quality.set_xticklabels([f"{x / 1000}k" for x in ax_quality.get_xticks()])
    ax_quality.set_yticklabels([])
    ax_length.set_xticklabels([])
    ax_length.set_yticklabels([f"{y / 1000}k" for y in ax_length.get_yticks()])
    ax.set_xlabel("Read length")
    ax.set_ylabel("Read quality")
    ax_length.set_ylabel("Number of reads")
    ax_quality.set_xlabel("Number of reads")
    fig.suptitle("Read length and quality distribution2")
    for f in format:
        fig.savefig(f"{prefix}__ReadLengthVsQualityMerged.{f}")


if __name__ == '__main__':
    args = get_arguments()
    df = pd.read_csv(args.input, sep="\t", header=0 if args.header else None, usecols=[args.length_col, args.quality_col])
    df.columns = ['length', 'quality']

    if args.mean_length is None:
        length_mean = np.mean(df['length'])
        # print("calculate mean length")
    else:
        length_mean = args.mean_length

    if args.plot_n50:
        if args.n50 is None:
            sorted_len = np.sort(df['length'])
            n50 = sorted_len[np.where(np.cumsum(sorted_len) >= 0.5 * np.sum(sorted_len))[0][0]]
            # print("calculate n50")
        else:
            n50 = args.n50
    else:
        n50 = None

    if args.std is None:
        std = np.std(df['length'])
        # print("calculate std")
    else:
        std = args.std


    fs = args.format if args.format else ["png"]
    min = length_mean - 3 * std
    max = length_mean + 3 * std
    df = df.query("@min <= length <= @max")
    mean_quality = args.mean_quality
    plot(df, args.plot_mean_length, length_mean,  args.plot_n50, n50, mean_quality, args.prefix, fs)
