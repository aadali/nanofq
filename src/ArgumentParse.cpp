#include "ArgumentParse.h"

#include <fmt/core.h>

using std::cout;
using std::cerr;
using std::endl;


argparse::ArgumentParser& get_arguments(int argc, char* argv[]) {
    static argparse::ArgumentParser nanofq{"nanofq", "1.0"};
    nanofq.add_description("A tool for stats, filter, index, find, trim nanopore fastq reads");
    nanofq.add_epilog("Contact aadali@gmail.com");

    static argparse::ArgumentParser main{"main"};
    main.add_description("stats, filter and trim barcode in one command");
    main.add_argument("-i", "--input")
        .help("the input fastq[.gz] or directory containing some fastq[.gz]")
        .required();
    main.add_argument("-o", "--output")
        .help("the output fastq")
        .required();
    main.add_argument("-p", "--prefix")
        .help("the prefix of the output file names")
        .required();
    main.add_argument("--retain_failed")
        .help("whether to save the reads that cannot pass the filters")
        .flag();
    main.add_argument("-n", "--firstN")
        .help("write the top N longest and high quality reads info in summary file, range (1, 1000)")
        .default_value(5)
        .scan<'i', int>();
    main.add_argument("-q", "--quality")
        .help( "count the reads number that whose quality is bigger than this value, can be set multi times, default vector<int>{10,12,15,18,20,25}, range (1, 50)")
        .append()
        .scan<'i', int>();
    main.add_argument("-l", "--length")
        .help(
            "count the reads number that longer than this value, can be set multi times, no default value, range >= 1")
        .append()
        .scan<'i', int>();
    main.add_argument("--plot")
        .help("whether plot the stats result")
        .flag();
    main.add_argument("-f", "--format")
        .help("what figure format you need, range {pdf, jpg, png}, only used when --plot used")
        .default_value<std::vector<std::string>>({"pdf"})
        .append();
    main.add_argument("--plot_mean_length")
        .help("whether plot the filtered mean length in plot, only used when --plot used")
        .flag();
    main.add_argument("--plot_n50")
        .help("whether plot the filtered n50 in plot, only used when --plot used")
        .flag();
    main.add_argument("-g", "--gc")
        .help("whether the stats gc content, if not set, all gc content will be set 0.0")
        .flag();
    main.add_argument("--min_len")
        .help("read min length, range (" + std::to_string(MINL) + ", " + std::to_string(MAXL) + ")")
        .default_value(MINL)
        .scan<'i', int>();
    main.add_argument("--max_len")
        .help("read max length, range (" + std::to_string(MINL) + ", " + std::to_string(MAXL) + "0")
        .default_value(MAXL)
        .scan<'i', int>();
    main.add_argument("--min_quality")
        .help("read min quality, range (0.0, 100.0)")
        .default_value(8.0f)
        .scan<'g', float>();
    main.add_argument("--min_gc")
        .help("read min gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(
                static_cast<int>(MIN_PERCENT)) + ", " +
            std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
        .default_value(MIN_PERCENT)
        .scan<'g', float>();
    main.add_argument("--max_gc")
        .help("read max gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(
                static_cast<int>(MIN_PERCENT)) + ", " +
            std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
        .default_value(MAX_PERCENT)
        .scan<'g', float>();
    auto& group0 = main.add_mutually_exclusive_group(false);
    group0.add_argument("--kit")
          .help(
              R"(the sequence kit name, used with --barcode if the sequence kit is barcoded. Each kit has it's own search parameter, but can be changed by [search parameter])")
          .choices("SQK-LSK114",
                   "SQK-RAD114",
                   "SQK-ULK114",
                   "SQK-PCS114",
                   "SQK-NBD114.24",
                   "SQK-NBD114.96",
                   "SQK-RBK114.24",
                   "SQK-RBK114.96",
                   "SQK-PCB114.24");
    group0.add_argument("--primers")
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
    main.add_argument("-b", "--barcode")
        .help("which barcode used, for 24 barcode kits, range (1, 24); for 96 barcode kits, range (1, 96)")
        .scan<'i', int>();
    main.add_argument("--match")
        .help("match score, positive int, range (-100, 0)")
        .default_value(3)
        .scan<'i', int>();
    main.add_argument("--mismatch")
        .help("mismatch score, negative int or 0, range (-100, 0)")
        .default_value(-3)
        .scan<'i', int>();
    main.add_argument("--gap_open")
        .help("gap opened score, negative int or 0, range (-100, 0)")
        .default_value(-7)
        .scan<'i', int>();
    main.add_argument("--gap_extend")
        .help("gap extend score, negative int or 0, range (-100, 0)")
        .default_value(-1)
        .scan<'i', int>();

    main.add_argument("--5end_len")
        .help("[search parameter]: check the first N bases from 5end of reads to find adapter, range (" +
            std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    main.add_argument("--5end_align_percent")
        .help("[search parameter]: the align length between query and reads 5end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    main.add_argument("--5end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

    main.add_argument("--3end_len")
        .help("[search parameter]: check the first N bases from 3end of reads to find adapter, range (" +
            std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    main.add_argument("--3end_align_percent")
        .help("[search parameter]: the align length between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    main.add_argument("--3end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

    main.add_argument("--5end_len_rc")
        .help(
            "[search parameter]: check the first N bases from 5end of reads to find adapter if this read is reverse complemented, range ("
            + std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    main.add_argument("--5end_align_percent_rc")
        .help(
            "[search parameter]: the align length between query and reads 5end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    main.add_argument("--5end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

    main.add_argument("--3end_len_rc")
        .help(
            "[search parameter]: check the first N bases from 3end of reads to find adapter if this read is reverse complemented, range ("
            + std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    main.add_argument("--3end_align_percent_rc")
        .help(
            "[search parameter]: the align length between query and reads 3end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    main.add_argument("--3end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    main.add_argument("-t", "--threads")
        .help("threads number used, range (1, 16)")
        .default_value(1)
        .scan<'i', int>();
    main.add_argument("-c", "--chunk")
        .help(
            "chunk number used, more chunk more memory, range (10000, 100000), bigger chunk, more memory will be used")
        .default_value(20000)
        .scan<'i', int>();
    nanofq.add_subparser(main);

    static argparse::ArgumentParser stats{"stats"};
    stats.add_description("stats nanopore fastq");
    stats.add_argument("-i", "--input")
         .help("the input fastq[.gz] or directory containing some fastq[.gz]")
         .required();
    stats.add_argument("-o", "--output")
         .help("the stats output file path")
         .required();
    stats.add_argument("-s", "--summary")
         .help("output the stats summary into this file")
         .default_value("./summary.txt");
    stats.add_argument("-n", "--firstN")
         .help("write the top N longest and high quality reads info in summary file, range (1, 1000)")
         .default_value(5)
         .scan<'i', int>();

    stats.add_argument("-q", "--quality")
         .help(
             "count the reads number that whose quality is bigger than this value, can be set multi times, default vector<int>{10,12,15,18,20,25}, range (1, 50)")
         .append()
         .scan<'i', int>();
    stats.add_argument("-l", "--length")
         .help(
             "count the reads number that longer than this value, can be set multi times, no default value, range >= 1")
         .append()
         .scan<'i', int>();
    stats.add_argument("-p", "--plot")
         .help("whether plot the stats result, if it's set, the value will be the figure file name prefix");
    stats.add_argument("-f", "--format")
         .help("what figure format you need, range {pdf, jpg, png}, only used when --plot used")
         .default_value<std::vector<std::string>>({"pdf"})
         .append();
    stats.add_argument("--plot_mean_length")
         .help("whether plot the mean length in plot, only used when --plot used")
         .flag();
    stats.add_argument("--plot_n50")
         .help("whether plot the n50 in plot, only used when --plot used")
         .flag();

    stats.add_argument("-g", "--gc")
         .help("whether the stats gc content, if not set, all gc content will be set 0.0")
         .flag();
    stats.add_argument("-t", "--threads")
         .help("threads number used, range (1, 16)").default_value(1)
         .scan<'i', int>();
    stats.add_argument("-c", "--chunk")
         .help(
             "chunk number used, more chunk more memory, range (10000, 100000), bigger chunk, more memory will be used")
         .default_value(20000)
         .scan<'i', int>();
    nanofq.add_subparser(stats);

    /*
    static argparse::ArgumentParser sample{"sample"};
    sample.add_description("randomly sampling specified number bases from fastq");
    sample.add_argument("-i", "--input")
        .help("the input fastq[.gz]")
        .required();
    sample.add_argument("-o", "--output")
        .help("the output fastq")
        .required();
    sample.add_argument("-n", "--base_number")
        .help( "how many bases you need when sample from --input, format can be like 5G/1.2g/5m/4.21G/3.12M. If the total data that meet the filter condition is less than this value, then use all the passed data")
        .required();
    sample.add_argument("-s", "--stats")
        .help("the stats result file from subcommand stats, if not set, will use subcommand stats to generate");
    sample.add_argument("-q", "--quality")
        .help("min read quality you need when sample from --input")
        .default_value(8)
        .scan<'i', int>();
    sample.add_argument("-l", "--min_len")
        .help("min read length you need when sample from --input")
        .default_value(MINL)
        .scan<'i', int>();
    sample.add_argument("-L", "--max_len")
        .help("max read length you need when sample from --input")
        .default_value(MAXL)
        .scan<'i', int>();
    sample.add_argument("-S", "--seed")
        .help("the random seed")
        .default_value(1)
        .scan<'i', int>();
    nanofq.add_subparser(sample);
    */

    static argparse::ArgumentParser filter{"filter"};
    filter.add_description("filter the input fastq");
    filter.add_argument("-i", "--input")
          .help("the input fastq[.gz] or directory containing some fastq[.gz]")
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
    filter.add_argument("-q", "--min_quality")
          .help("read min quality, range (0.0, 100.0)")
          .default_value(8.0f)
          .scan<'g', float>();
    filter.add_argument("", "--gc")
          .help("whether filter the gc content, used with --min_gc/--max_gc")
          .flag();
    filter.add_argument("-g", "--min_gc")
          .help("read min gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(
                  static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MIN_PERCENT)
          .scan<'g', float>();
    filter.add_argument("-G", "--max_gc")
          .help("read max gc content, used with --gc, otherwise ignore this parameter, range (" + std::to_string(
                  static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MAX_PERCENT)
          .scan<'g', float>();
    filter.add_argument("-t", "--threads")
          .help("threads number used, range (1, 16)")
          .default_value(1)
          .scan<'i', int>();
    filter.add_argument("-c", "--chunk")
          .help(
              "chunk number used, more chunk more memory, range (10000, 100000), bigger chunk, more memory will be used")
          .default_value(20000)
          .scan<'i', int>();
    nanofq.add_subparser(filter);

    static argparse::ArgumentParser index{"index"};
    index.add_description("index the input fastq");
    index.add_argument("input")
         .help("the input fastq[.gz]")
         .required();
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
    nanofq.add_subparser(find);

    static argparse::ArgumentParser compress{"compress"};
    compress.add_description("convert text fastq to NanoBgzip fastq and build an index at same time");
    compress.add_argument("output")
            .help("the output NanoBgzip fastq, must endswith .gz")
            .required();
    compress.add_argument("input")
            .help("the input text fastq, default -, means from stdin")
            .default_value("-");
    compress.add_argument("-n", "--number")
            .help("compress how many reads into each block")
            .default_value(10)
            .scan<'i', int>();
    nanofq.add_subparser(compress);

    static argparse::ArgumentParser trim{"trim"};
    trim.add_description("use local alignment to find possible adapter, barcode, primers in reads and trim them");
    trim.add_argument("-i", "--input")
        .help("the input fastq[.gz] or directory containing some fastq[.gz]")
        .required();
    trim.add_argument("-o", "--output")
        .help("the output file")
        .default_value("-");
    trim.add_argument("-l", "--log")
        .help("the log output file")
        .default_value("./trim_log.txt");
    auto& group = trim.add_mutually_exclusive_group(true);
    group.add_argument("-k", "--kit")
         .help(
             R"(the sequence kit name, used with --barcode if the sequence kit is barcoded. Each kit has it's own search parameter, but can be changed by [search parameter])")
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
        .help(
            "chunk number used, more chunk more memory, range (10000, 100000), bigger chunk, more memory will be used")
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
        .scan<'g', float>();
    trim.add_argument("--5end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

    trim.add_argument("--3end_len")
        .help("[search parameter]: check the first N bases from 3end of reads to find adapter, range (" +
            std::to_string(MIN_TARGET) + ", " + std::to_string(MAX_TARGET) + ")")
        .scan<'i', int>();
    trim.add_argument("--3end_align_percent")
        .help("[search parameter]: the align length between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    trim.add_argument("--3end_align_identity")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

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
        .scan<'g', float>();
    trim.add_argument("--5end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 5end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();

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
        .scan<'g', float>();
    trim.add_argument("--3end_align_identity_rc")
        .help(
            "[search parameter]: the align identity between query and reads 3end should bigger than this value if this read is reverse complemented, range ("
            + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " + std::to_string(static_cast<int>(MAX_PERCENT)) +
            ")")
        .scan<'g', float>();
    nanofq.add_subparser(trim);

    try {
        nanofq.parse_args(argc, argv);
    }
    catch (const std::exception& e) {
        cerr << "\n" << REDS << e.what() << COLOR_END << "\n" << endl;
        // if (nanofq.is_subcommand_used("stats")) {
        //     cerr << stats << endl;
        // } else if (nanofq.is_subcommand_used("trim")) {
        //     // cerr << trim << endl;
        // } else if (nanofq.is_subcommand_used("find")) {
        //     cerr << find << endl;
        // } else if (nanofq.is_subcommand_used("index")) {
        //     cerr << index << endl;
        // } else if (nanofq.is_subcommand_used("filter")) {
        //     cerr << filter << endl;
        // } else {
        //     cerr << nanofq << endl;
        // }
        exit(1);
    }
    return nanofq;
}
