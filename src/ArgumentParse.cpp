#include "ArgumentParse.h"
#include <fmt/core.h>

using std::cout;
using std::cerr;
using std::endl;


argparse::ArgumentParser& get_arguments(int argc, char* argv[]) {
    static argparse::ArgumentParser nanofq{"nanofq", "1.0"};
    nanofq.add_description("A tool for stats, filter, index, find, trim nanopore fastq reads");
    nanofq.add_epilog("Contact aadali@gmail.com");
    static argparse::ArgumentParser stats{"stats"};
    stats.add_description("stats nanopore fastq");

    stats.add_argument("-i", "--input")
         .help("the input fastq[.gz]")
         .required();
    stats.add_argument("-o", "--output")
         .help("the stats output file name, default print all results to stdout")
         .default_value("-");
    stats.add_argument("-s", "--summary")
        .help("output the summary into this file")
        .default_value("./summary.txt");
    stats.add_argument("-n", "--firstN")
        .help("get the top N longest and high quality reads info in summary, range (1, 1000)")
        .default_value(5)
        .scan<'i', int>();

    stats.add_argument("-q", "--quality")
        .help("count the reads number that whose quality is bigger than this value, can be set multi times, range (1, 50)")
        .append()
        .scan<'i', int>();

    stats.add_argument("-p", "--plot")
         .help("whether plot the stats result, if it's set, the value will be the figure file name prefix");
    stats.add_argument("-f", "--format")
         .help("what figure format you need, range {pdf, jpg, png}")
         .default_value<std::vector<std::string>>({"pdf"})
         .append();
    stats.add_argument("-g", "--gc")
         .help("whether the stats gc content, if not set, all gc content will be set 0.0")
         .flag();
    stats.add_argument("-t", "--threads")
         .help("threads number used, range (1, 16)").default_value(1)
         .scan<'i', int>();
    stats.add_argument("-c", "--chunk")
         .help("chunk number used, more chunk more memory, range (10000, 100000)")
         .default_value(20000)
         .scan<'i', int>();
    nanofq.add_subparser(stats);

    static argparse::ArgumentParser filter{"filter"};
    filter.add_description("filter the input fastq[.gz]");
    filter.add_argument("-i", "--input")
          .help("the input fastq[.gz]")
          .required();
    filter.add_argument("-o", "--output")
          .help("the stats output file name, default print all results to stdout")
          .default_value("-");
    filter.add_argument("-l", "--min_len")
          .help("read min length, range (" + std::to_string(MINL) + ", " + std::to_string(MAXL) + ")")
          .default_value(MINL)
          .scan<'i', int>();
    filter.add_argument("-L", "--max_len")
          .help("read max length, range (" + std::to_string(MINL) + ", " + std::to_string(MAXL) + ")")
          .default_value(MAXL)
          .scan<'i', int>();
    filter.add_argument("-q", "--min_quality, range (0.0, 100.0)")
          .help("read min quality")
          .default_value(8.0)
          .scan<'g', double>();
    filter.add_argument("", "--gc")
          .help("whether filter the gc content, used with --min_gc/--max_gc")
          .flag();
    filter.add_argument("-g", "--min_gc")
          .help("read min gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MIN_PERCENT)
          .scan<'g', double>();
    filter.add_argument("-G", "--max_gc")
          .help("read max gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MAX_PERCENT)
          .scan<'g', double>();
    filter.add_argument("-t", "--threads")
          .help("threads number used, range (1, 16)")
          .default_value(1)
          .scan<'i', int>();
    filter.add_argument("-c", "--chunk")
          .help("chunk number used, more chunk more memory, range (10000, 100000)")
          .default_value(20000)
          .scan<'i', int>();
    nanofq.add_subparser(filter);

    static argparse::ArgumentParser index{"index"};
    index.add_description("index the input fastq. Never do it for short reads");
    index.add_argument("-i", "--input")
         .help("the input fastq[.gz]")
         .required();
    index.add_argument("-k", "--key_len")
         .help(
             "the first N of read name used to make index, if this value bigger than readname, use readname as key, range >=8")
         .default_value(8)
         .scan<'i', int>();
    nanofq.add_subparser(index);

    static argparse::ArgumentParser find{"find"};
    find.add_description("find specified records from input fastq depends on the specified readnames");
    find.add_argument("-i", "--input")
        .help("the input fastq[.gz]")
        .required();
    find.add_argument("-o", "--output")
        .help("the output file name, default print all results to stdout")
        .default_value("-");
    find.add_argument("-r", "--reads")
        .help(
            "readnames of records expected in fastq, multi readnames separated by comma. Or file containing readnames, one read name per line")
        .required();
    find.add_argument("-u", "--use_index")
        .help(
            "whether use index to find reads, if true and no index file exists, it will make index firstly. Index fastq may need take a while, but once the index had finished, searching reads will be very fast")
        .flag();
    find.add_argument("-k", "--key_len")
        .help(
            "the first N of read name used to make index, if this value bigger than readname, use readname as key, range >=8")
        .default_value(8)
        .scan<'i', int>();
    nanofq.add_subparser(find);

    static argparse::ArgumentParser trim{"trim"};
    trim.add_description("Use local alignment to find possible adapter, barcode, primers in reads and trim them");
    trim.add_argument("-i", "--input").help("the input fastq[.gz]").required();
    trim.add_argument("-o", "--output").help("the output file").default_value("-");
    trim.add_argument("-l", "--log").help("the log output file").default_value("./log.txt");
    auto& group = trim.add_mutually_exclusive_group(true);
    group.add_argument("-k", "--kit")
         .help(R"(the sequence kit name, used with --barcode if the sequence kit is barcoded.
Each kit has it's own search parameter, but can be changed by [search parameter])")
         .choices("SQK-LSK114",
                  "SQK-RAD114",
                  "SQK-ULK114",
                  "SQK-PCS114",
                  "SQK-NBD114.24",
                  "SQK-NBD114.96",
                  "SQK-RBK114.24",
                  "SQK-RBK114.96",
                  "SQK-PCB114.24");

    group.add_argument("-p", "--primers")
         .help(
             R"(customer forward and reversed primers joined by comma or file containing two line which first line is forward primer and second is reversed.
if it's set, the following parameter will be set with default value
    5end_len = 180
    5end_align_percent = 0.8
    5end_align_identity = 0.8
    3end_len = 180
    3end_align_percent = 0.8
    3end_align_identity = 0.8
    5end_len_rc = 180
    5end_align_percent_rc = 0.8
    5end_align_identity_rc = 0.8
    3end_len_rc = 180
    3end_align_percent_rc = 0.8
    3end_align_identity_rc = 0.8
you can change this value by set specified parameter)");

    trim.add_argument("-b", "--barcode")
        .help("which barcode used, for 24 barcode kits, range (1, 24); for 96 barcode kits, range (1, 96)")
        .scan<'i', int>();


    trim.add_argument("-t", "--threads")
        .help("threads number used, range (1, 16)")
        .default_value(1)
        .scan<'i', int>();
    trim.add_argument("-c", "--chunk")
        .help("chunk number used, more chunk more memory, range (10000, 100000)")
        .default_value(20000)
        .scan<'i', int>();

    trim.add_argument("--match")
        .help("match score, positive int, range (-100, 0)")
        .default_value(3)
        .scan<'i', int>();
    trim.add_argument("--mismatch")
        .help("mismatch score, negative int or 0, range (-100, 0)")
        .default_value(-3)
        .scan<'i', int>();
    trim.add_argument("--gap_open")
        .help("gap opened score, negative int or 0, range (-100, 0)")
        .default_value(-7)
        .scan<'i', int>();
    trim.add_argument("--gap_extend")
        .help("gap extend score, negative int or 0, range (-100, 0)")
        .default_value(-1)
        .scan<'i', int>();

    trim.add_argument("--5end_len")
        .help("[search parameter]: check the first N bases from 5end of reads to find adapter, range (" +
            std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    trim.add_argument("--5end_align_percent")
        .help("[search parameter]: the align length between query and reads 5end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();
    trim.add_argument("--5end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();

    trim.add_argument("--3end_len")
        .help("[search parameter]: check the first N bases from 3end of reads to find adapter, range (" +
            std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    trim.add_argument("--3end_align_percent")
        .help("[search parameter]: the align length between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();
    trim.add_argument("--3end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();

    trim.add_argument("--5end_len_rc")
        .help(
            "[search parameter]: check the first N bases from 5end of reads to find adapter if this read is reverse complemented, range ("
            + std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    trim.add_argument("--5end_align_percent_rc")
        .help(
            "[search parameter]: the align length between query and reads 5end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();
    trim.add_argument("--5end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();

    trim.add_argument("--3end_len_rc")
        .help(
            "[search parameter]: check the first N bases from 3end of reads to find adapter if this read is reverse complemented, range ("
            + std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    trim.add_argument("--3end_align_percent_rc")
        .help(
            "[search parameter]: the align length between query and reads 3end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();
    trim.add_argument("--3end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', double>();
    nanofq.add_subparser(trim);

    try {
        nanofq.parse_args(argc, argv);
    }
    catch (const std::exception& e) {
        cerr << REDS << e.what() << COLOR_END << endl;
        if (nanofq.is_subcommand_used("stats")) {
            cerr << stats << endl;
        } else if (nanofq.is_subcommand_used("trim")) {
            cerr << trim << endl;
        } else if (nanofq.is_subcommand_used("find")) {
            cerr << find << endl;
        } else if (nanofq.is_subcommand_used("index")) {
            cerr << index << endl;
        } else if (nanofq.is_subcommand_used("filter")) {
            cerr << filter << endl;
        } else {
            cerr << nanofq << endl;
        }
        exit(1);
    }
    return nanofq;
}
