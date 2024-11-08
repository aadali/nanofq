#include "parse_arguments.h"
#include <fmt/core.h>

using std::cout;
using std::cerr;
using std::endl;


template <typename T>
void check_number_in_range(const std::string& parameter, const T& number, T min, T max,
                           argparse::ArgumentParser& parser, bool integer) {
    std::string type{integer ? "integer" : "float"};
    if (number < min || number > max) {
        cerr << fmt::format("{} should be a {}, and in range ({}, {})", parameter, type, min, max) << endl;
        cerr << parser << endl;
        exit(1);
    }
}

template <typename T>
void check_choices(const std::string& parameter, std::vector<T>& choices, std::vector<T>& allowed_choices,
                   argparse::ArgumentParser& parser) {
    for (T& candidate : choices) {
        if (std::find(allowed_choices.begin(), allowed_choices.end(), candidate) == allowed_choices.end()) {
            std::cerr << fmt::format("{} allowed choice should be in [{}]", parameter,
                                     join<T>(", ", allowed_choices)) << '\n';
            std::cerr << parser << endl;
            exit(1);
        }
    }
}

argparse::ArgumentParser& get_arguments(int argc, char* argv[]) {
    static argparse::ArgumentParser nanofq{"nanofq", "1.0"};
    nanofq.add_description("A tool for treat Nanopore fastq[.gz]");
    nanofq.add_epilog("aadali@gmail.com");
    argparse::ArgumentParser stats{"stats"};
    stats.add_description("stats the Nanopore fastq[.gz]");

    stats.add_argument("-i", "--input")
         .help("the input fastq[.gz]")
         .required();
    stats.add_argument("-o", "--output")
         .help("the stats output file name")
         .default_value("-");
    stats.add_argument("-p", "--plot")
         .help("whether plot the stats result, if it's set, the value will be the figure file name prefix");
    stats.add_argument("-f", "--format")
         .help("what figure format you need")
         .default_value<std::vector<std::string>>({"pdf"})
         .append();
    stats.add_argument("-g", "--gc")
         .help("whether the stats gc content, if not set, all gc content will be set 0.0")
         .flag();
    stats.add_argument("-t", "--threads")
         .help("threads number used").default_value(4)
         .scan<'i', int>();
    stats.add_argument("-c", "--chunk")
         .help("chunk number used, more chunk more memory")
         .default_value(20000)
         .scan<'i', int>();
    nanofq.add_subparser(stats);

    argparse::ArgumentParser filter{"filter"};
    filter.add_description("filter the input fastq[.gz]");
    filter.add_argument("-i", "--input")
          .help("the input fastq[.gz]")
          .required();
    filter.add_argument("-o", "--output")
          .help("the stats output file name")
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
          .help("read min gc content, used with --gc, range (" + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MIN_PERCENT)
          .scan<'g', double>();
    filter.add_argument("-G", "--max_gc")
          .help("read max gc content, used with --gc, range (" + std::to_string(static_cast<int>(MIN_PERCENT)) + ", " +
              std::to_string(static_cast<int>(MAX_PERCENT)) + ")")
          .default_value(MAX_PERCENT)
          .scan<'g', double>();
    filter.add_argument("-t", "--threads")
          .help("threads number used, range (1, 16)")
          .default_value(4)
          .scan<'i', int>();
    filter.add_argument("-c", "--chunk")
          .help("chunk number used, more chunk more memory, range (10000, 100000)")
          .default_value(20000)
          .scan<'i', int>();
    nanofq.add_subparser(filter);

    argparse::ArgumentParser index{"index"};
    index.add_description("index the input fastq");
    index.add_argument("-i", "--input")
         .help("the input fastq[.gz]")
         .required();
    index.add_argument("-k", "--key_len")
         .help(
             "the first N of read name used to make index, if this value bigger than readname, use readname as key, range >=8")
         .default_value(8)
         .scan<'i', int>();
    nanofq.add_subparser(index);

    argparse::ArgumentParser find{"find"};
    find.add_description("find specified record from input fastq depends on the specified readnames");
    find.add_argument("-i", "--input")
        .help("the input fastq[.gz]")
        .required();
    find.add_argument("-o", "--output")
        .help("the output file name, default, output to stdout")
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

    argparse::ArgumentParser trim{"trim"};
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
        .default_value(4)
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
        cerr << e.what() << endl;
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
    if (nanofq.is_subcommand_used("stats")) {
        cout << "=================stats================" << endl;
        auto input{stats.get("--input")};
        cout << "input: " << input << endl;;
        auto output{stats.get("--output")};
        cout << "output: " << output << endl;
        if (stats.is_used("--plot")) {
            auto plot{stats.get("--plot")};
            cout << "plot: " << plot << endl;
        } else {
            cout << "plot: " << "not set" << endl;
        }
        auto threads{stats.get<int>("--threads")};
        check_number_in_range("--threads", threads, MINT, MAXT, stats, true);
        cout << "threads: " << threads << endl;
        auto chunk{stats.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, stats, true);
        cout << "chunk: " << chunk << endl;
    } else if (nanofq.is_subcommand_used("filter")) {
        cout << "=================filter================" << endl;
        auto input{filter.get("--input")};
        cout << "input: " << input << endl;;
        auto output{filter.get("--output")};
        cout << "output: " << output << endl;
        auto min_len{filter.get<int>("--min_len")};
        auto max_len{filter.get<int>("--max_len")};
        check_number_in_range("--min_len", min_len, MINL, MAXL, filter, true);
        check_number_in_range("--max_len", max_len, MINL, MAXL, filter, true);
        cout << "min_len: " << min_len << endl;
        cout << "max_len: " << max_len << endl;
        if (filter.get<bool>("--gc")) {
            auto min_gc = filter.get<double>("--min_gc");
            auto max_gc = filter.get<double>("--max_gc");
            check_number_in_range("--min_gc", min_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
            check_number_in_range("--max_gc", max_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
            cout << "min_gc: " << min_gc << endl;
            cout << "max_gc: " << max_gc << endl;
        } else {
            cout << "gc: " << "not set" << endl;
        }
        auto threads{filter.get<int>("--threads")};
        check_number_in_range("--threads", threads, MINT, MAXT, filter, true);
        cout << "threads: " << threads << endl;
        auto chunk{filter.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, filter, true);
        cout << "chunk: " << chunk << endl;
    } else if (nanofq.is_subcommand_used("index")) {
        // TODO
    } else if (nanofq.is_subcommand_used("trim")) {
        std::string input = trim.get("--input");
        cout << "input: " << input << endl;
        std::string output = trim.get("--output");
        cout << "output: " << output << endl;
        auto threads{trim.get<int>("--threads")};
        check_number_in_range("--threads", threads, MINT, MAXT, trim, true);
        cout << "threads: " << threads << endl;
        auto chunk{trim.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, trim, true);
        cout << "chunk: " << chunk << endl;
        if (trim.is_used("--kit")) {
            auto kit = trim.get("--kit");
            cout << "kit: " << kit << endl;
            if (trim.is_used("--barcode")) {
                auto barcode = trim.get<int>("--barcode");
                if (kit.ends_with(".24")) {
                    if (barcode < MINB || barcode > MAX24B) {
                        cerr << "If kit with 24 barcodes used, --barcode should be a integer and  in range (1, 24)" <<
                            endl;
                        cerr << trim << endl;
                        exit(1);
                    }
                    cout << "barcode: " << barcode << endl;
                } else if (kit.ends_with(".96")) {
                    if (barcode < MINB || barcode > MAX96B) {
                        cerr << "If kit with 96 barcodes used, --barcode should be a integer and  in range (1, 96)" <<
                            endl;
                        cerr << trim << endl;
                        exit(1);
                    }
                    cout << "barcode: " << barcode << endl;
                } else {
                    cerr << "If kit with no barcode used, ignore --barcode" << endl;
                }
            }
        } else {
            auto primer = trim.get("--primers");
            cout << "primer: " << primer << endl;
            if (trim.is_used("--barcode")) {
                cerr << "If --kit no used, --barcode will be ignored" << endl;
            }
        }
        int match{trim.get<int>("--match")};
        check_number_in_range<int>("--match", match, 1, 100, trim, true);
        cout << "match: " << match << endl;
        int mismatch{trim.get<int>("--mismatch")};
        check_number_in_range<int>("--mismatch", mismatch, -100, 0, trim, true);
        cout << "mismatch: " << mismatch << endl;
        int gap_open{trim.get<int>("--gap_open")};
        check_number_in_range<int>("--gap_open", gap_open, -100, 0, trim, true);
        cout << "gap_open: " << gap_open << endl;
        int gap_extend{trim.get<int>("--gap_extend")};
        check_number_in_range<int>("--gap_extend", gap_extend, -100, 0, trim, true);
        cout << "gap_extend: " << gap_extend << endl;
        if (trim.is_used("--5end_len")) {
            int end5_len{trim.get<int>("--5end_len")};
            check_number_in_range("--5end_len", end5_len, MIN_TARGET, MAX_TARGET, trim, true);
            cout << "5end_len: " << end5_len << endl;
        }
        if (trim.is_used("--5end_align_percent")) {
            double end5_align_percent{trim.get<double>("--5end_align_percent")};
            check_number_in_range("--5end_align_percent", end5_align_percent, MIN_PERCENT, MAX_PERCENT, trim, false);
            cout << "5end_align_percent: " << end5_align_percent << endl;
        }
        if (trim.is_used("--5end_align_identity")) {
            double end5_align_identity{trim.get<double>("--5end_align_identity")};
            check_number_in_range("--5end_align_identity", end5_align_identity, MIN_PERCENT, MAX_PERCENT, trim, false);
            cout << "5end_align_identity: " << end5_align_identity << endl;
        }
        if (trim.is_used("--3end_len")) {
            int end3_len{trim.get<int>("--3end_len")};
            check_number_in_range("--3end_len", end3_len, MIN_TARGET, MAX_TARGET, trim, true);
            cout << "3end_len: " << end3_len << endl;
        }
        if (trim.is_used("--3end_align_percent")) {
            double end3_align_percent{trim.get<double>("--3end_align_percent")};
            check_number_in_range("--3end_align_percent", end3_align_percent, MIN_PERCENT, MAX_PERCENT, trim, false);
            cout << "3end_align_percent: " << end3_align_percent << endl;
        }
        if (trim.is_used("--3end_align_identity")) {
            double end3_align_identity{trim.get<double>("--3end_align_identity")};
            check_number_in_range("--3end_align_identity", end3_align_identity, MIN_PERCENT, MAX_PERCENT, trim, false);
            cout << "3end_align_identity: " << end3_align_identity << endl;
        }
        if (trim.is_used("--5end_len_rc")) {
            int end5_len_rc{trim.get<int>("--5end_len_rc")};
            check_number_in_range("--5end_len_rc", end5_len_rc, MIN_TARGET, MAX_TARGET, trim, true);
            cout << "5end_len_rc: " << end5_len_rc << endl;
        }
        if (trim.is_used("--5end_align_percent_rc")) {
            double end5_align_percent_rc{trim.get<double>("--5end_align_percent_rc")};
            check_number_in_range("--5end_align_percent_rc", end5_align_percent_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
            cout << "5end_align_percent_rc: " << end5_align_percent_rc << endl;
        }
        if (trim.is_used("--5end_align_identity_rc")) {
            double end5_align_identity_rc{trim.get<double>("--5end_align_identity_rc")};
            check_number_in_range("--5end_align_identity_rc", end5_align_identity_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
            cout << "5end_align_identity_rc: " << end5_align_identity_rc << endl;
        }
        if (trim.is_used("--3end_len_rc")) {
            int end3_len_rc{trim.get<int>("--3end_len_rc")};
            check_number_in_range("--3end_len_rc", end3_len_rc, MIN_TARGET, MAX_TARGET, trim, true);
            cout << "3end_len_rc: " << end3_len_rc << endl;
        }
        if (trim.is_used("--3end_align_percent_rc")) {
            double end3_align_percent_rc{trim.get<double>("--3end_align_percent_rc")};
            check_number_in_range("--3end_align_percent_rc", end3_align_percent_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
            cout << "3end_align_percent_rc: " << end3_align_percent_rc << endl;
        }
        if (trim.is_used("--3end_align_identity_rc")) {
            double end3_align_identity_rc{trim.get<double>("--3end_align_identity_rc")};
            check_number_in_range("--3end_align_identity_rc", end3_align_identity_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
            cout << "3end_align_identity_rc: " << end3_align_identity_rc << endl;
        }
    }
    return nanofq;
}
