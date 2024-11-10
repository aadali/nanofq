#ifndef SUBMAIN_H
#define SUBMAIN_H
#include <thread>
#include <filesystem>

#include "FastqReader.h"
#include "Work.h"
#include "SequenceInfo.h"
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
#include "Adapter.h"
#include "ArgumentParse.h"

int sub_main(int argc, char* argv[]) {
    argparse::ArgumentParser& nanofq{get_arguments(argc, argv)};
    if (nanofq.is_subcommand_used("stats")) {
        argparse::ArgumentParser& stats{nanofq.at<argparse::ArgumentParser>("stats")};
        std::string input{stats.get("--input")};
        std::string output{stats.get("--output")};
        std::string summary{stats.get("--summary")};
        int n{stats.get<int>("--firstN")};
        check_number_in_range("--firstN", n, 1, 1000, stats, true);
        std::vector<int> quals;
        if (!stats.is_used("--quality")) {
            // quals = {9, 12, 15, 18, 20, 25};
            quals = {25, 20, 18, 15, 12, 9};
        } else {
            quals = {stats.get<std::vector<int>>("--quality")};
            for (int i : quals) {
                check_number_in_range("--quality", i, 1, 50, stats, true);
            }
            std::ranges::sort(quals, greater<>());
        }
        bool make_plot{false};
        std::string plot_prefix;
        if (stats.is_used("--plot")) {
            make_plot = true;
            plot_prefix = stats.get("--plot");
        }
        int threads{stats.get<int>("--threads")};
        check_number_in_range("--threads", threads, 1, 16, stats, true);
        int chunk{stats.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, stats, true);
        bool gc{stats.get<bool>("--gc")};
        std::vector<std::string> format{stats.get<std::vector<std::string>>("--format")};
        std::vector<std::string> allowed_choices{"pdf", "jpg", "png"};
        check_choices<std::string>("--format", format, allowed_choices, stats);
        FastqReader fq{input, static_cast<unsigned>(chunk)};
        std::ofstream out;
        if (output != "-") {
            out.open(output.data(), std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        Work work{fq, static_cast<unsigned>(threads), gc, output == "-" ? std::cout : out};
        if (threads == 1) {
            work.run_stats();
        } else {
            std::thread t1{&FastqReader::read_chunk_fastq, &fq};
            std::thread t2{&Work::run_stats, &work};
            t1.join();
            t2.join();
        }
        work.save_summary(n, quals, summary);
        if (out.is_open()) { out.close(); }
    } else if (nanofq.is_subcommand_used("filter")) {
        argparse::ArgumentParser& filter{nanofq.at<argparse::ArgumentParser>("filter")};
        std::string input{filter.get("--input")};
        std::string output{filter.get("--output")};
        int min_length{filter.get<int>("--min_length")};
        check_number_in_range("--min_length", min_length, MINL, MAXL, filter, true);
        int max_length{filter.get<int>("--max_length")};
        check_number_in_range("--max_length", max_length, MINL, MAXL, filter, true);
        double min_quality{filter.get<double>("--min_quality")};
        check_number_in_range("--min_double", min_quality, 0.0, 100.0, filter, false);
        bool gc{filter.get<bool>("--gc")};
        double min_gc{filter.get<double>("--min_gc")};
        check_number_in_range("--min_gc", min_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
        double max_gc{filter.get<double>("--max_gc")};
        check_number_in_range("--max_gc", max_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
        int threads{filter.get<int>("--threads")};
        check_number_in_range("--threads", threads, 1, 16, filter, true);
        int chunk{filter.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, filter, true);
        FastqReader fq{input, static_cast<unsigned>(chunk)};
        std::ofstream out;
        if (output == "-") {
            out = std::ofstream{output.data(), std::ios::out};
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        Work work{fq, static_cast<unsigned>(threads), gc, output == "-" ? std::cout : out};
        if (threads == 1) {
            work.run_filter(min_length, max_length, min_quality, min_gc, max_gc);
        } else {
            std::thread t1{&FastqReader::read_chunk_fastq, &fq};
            std::thread t2{&Work::run_filter, &work, min_length, max_length, min_quality, min_gc, max_gc};
            t1.join();
            t2.join();
        }
        if (out.is_open()) out.close();
    } else if (nanofq.is_subcommand_used("index")) {
        argparse::ArgumentParser& index{nanofq.at<argparse::ArgumentParser>("index")};
        std::string input{index.get("--input")};
        int key_len{index.get<int>("--key_len")};
        check_number_in_range("--key_len", key_len, 8, 100, index, true);
        FastqReader fq{input, 5000};
        Work work{fq};
        work.run_index(key_len);
    } else if (nanofq.is_subcommand_used("find")) {
        argparse::ArgumentParser& find{nanofq.at<argparse::ArgumentParser>("find")};
        std::string input{find.get("--input")};
        std::string output{find.get("--output")};
        std::string reads{find.get("--reads")};
        bool use_index{find.get<bool>("--use_index")};
        int key_len;
        if (!use_index) {
            if (find.is_used("--key_len")) {
                cerr << WARNS + "if --use_index is not set, ignore --key_len" + COLOR_END << endl;
            }
        } else {
            key_len = find.get<int>("--key_len");
            check_number_in_range("--key_len", key_len, 8, 100, find, true);
        }
        FastqReader fq{input, 5000};
        std::ofstream out;
        if (output != "-") {
            out.open(output, std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        Work work{fq, 1, true, output == "-" ? std::cout : out};
        work.run_find(reads, use_index, key_len);
        if (out.is_open()) out.close();
    } else if (nanofq.is_subcommand_used("trim")) {
        argparse::ArgumentParser& trim{nanofq.at<argparse::ArgumentParser>("trim")};
        std::string input{trim.get("--input")};
        std::string output{trim.get("--output")};
        std::string log{trim.get("--log")};
        auto threads{trim.get<int>("--threads")};
        check_number_in_range("--threads", threads, MINT, MAXT, trim, true);
        auto chunk{trim.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, trim, true);
        std::string kit;
        std::string forward;
        std::string reversed;
        int barcode;
        if (trim.is_used("--kit")) {
            kit = trim.get("--kit");
            if (kit.ends_with(".24") || kit.ends_with(".96")) {
                if (!trim.is_used("--barcode")) {
                    cerr << REDS + "If kit with barcodes used, --barcode must be set" + COLOR_END << endl;
                    cerr << trim << endl;
                    exit(1);
                }
            }
            if (trim.is_used("--barcode")) {
                barcode = trim.get<int>("--barcode");
                if (kit.ends_with(".24")) {
                    if (barcode < MINB || barcode > MAX24B) {
                        cerr << REDS +
                            "If kit with 24 barcodes used, --barcode should be a integer and  in range (1, 24)" +
                            COLOR_END <<
                            endl;
                        cerr << trim << endl;
                        exit(1);
                    }
                    kit = fmt::format("{}-{}", kit, barcode);
                } else if (kit.ends_with(".96")) {
                    if (barcode < MINB || barcode > MAX96B) {
                        cerr << REDS +
                            "If kit with 96 barcodes used, --barcode should be a integer and  in range (1, 96)" +
                            COLOR_END <<
                            endl;
                        cerr << trim << endl;
                        exit(1);
                    }
                    kit = fmt::format("{}-{}", kit, barcode);
                } else {
                    cerr << WARNS + "If kit with no barcode used, ignore --barcode" + COLOR_END << endl;
                }
            }
        } else {
            std::string primers = trim.get("--primers");
            if (exists(std::filesystem::path{primers.data()})) {
                std::ifstream primers_file{primers, std::ios::in};
                std::getline(primers_file, forward);
                std::getline(primers_file, reversed);
                primers_file.close();
            } else {
                auto primers_vec = myUtility::split(primers, ",");
                if (primers_vec.size() != 2) {
                    cerr << REDS + "if --primer is not file, it should be a pair of primers separated by one comma" +
                        COLOR_END << endl;
                    cerr << trim << endl;
                    exit(1);
                }
                forward = std::string{primers_vec[0]};
                reversed = std::string{primers_vec[1]};
            }
        }
        int match{trim.get<int>("--match")};
        check_number_in_range<int>("--match", match, 1, 100, trim, true);
        int mismatch{trim.get<int>("--mismatch")};
        check_number_in_range<int>("--mismatch", mismatch, -100, 0, trim, true);
        int gap_open{trim.get<int>("--gap_open")};
        check_number_in_range<int>("--gap_open", gap_open, -100, 0, trim, true);
        int gap_extend{trim.get<int>("--gap_extend")};
        check_number_in_range<int>("--gap_extend", gap_extend, -100, 0, trim, true);
        int end5_len = DEFAULT_INT, end3_len = DEFAULT_INT, end5_len_rc = DEFAULT_INT, end3_len_rc = DEFAULT_INT;
        double end5_align_percent = DEFAULT_FLOAT, end5_align_identity = DEFAULT_FLOAT;
        double end3_align_percent = DEFAULT_FLOAT, end3_align_identity = DEFAULT_FLOAT;
        double end5_align_percent_rc = DEFAULT_FLOAT, end5_align_identity_rc = DEFAULT_FLOAT;
        double end3_align_percent_rc = DEFAULT_FLOAT, end3_align_identity_rc = DEFAULT_FLOAT;
        if (trim.is_used("--5end_len")) {
            end5_len = trim.get<int>("--5end_len");
            check_number_in_range("--end5_len", end5_len, MIN_TARGET, MAX_TARGET, trim, true);
        }
        if (trim.is_used("--5end_align_percent")) {
            end5_align_percent = trim.get<double>("--5end_align_percent");
            check_number_in_range("--5end_align_percent", end5_align_percent, MIN_PERCENT, MAX_PERCENT, trim, false);
        }
        if (trim.is_used("--5end_align_identity")) {
            end5_align_identity = trim.get<double>("--5end_align_identity");
            check_number_in_range("--5end_align_identity", end5_align_identity, MIN_PERCENT, MAX_PERCENT, trim, false);
        }
        if (trim.is_used("--3end_len")) {
            end3_len = trim.get<int>("--3end_len");
            check_number_in_range("--3end_len", end3_len, MIN_TARGET, MAX_TARGET, trim, true);
        }
        if (trim.is_used("--3end_align_percent")) {
            end3_align_percent = trim.get<double>("--3end_align_percent");
            check_number_in_range("--3end_align_percent", end3_align_percent, MIN_PERCENT, MAX_PERCENT, trim, false);
        }
        if (trim.is_used("--3end_align_identity")) {
            end3_align_identity = trim.get<double>("--3end_align_identity");
            check_number_in_range("--3end_align_identity", end3_align_identity, MIN_PERCENT, MAX_PERCENT, trim, false);
        }
        if (trim.is_used("--5end_len_rc")) {
            end5_len_rc = trim.get<int>("--5end_len_rc");
            check_number_in_range("--5end_len_rc", end5_len_rc, MIN_TARGET, MAX_TARGET, trim, true);
        }
        if (trim.is_used("--5end_align_percent_rc")) {
            end5_align_percent_rc = trim.get<double>("--5end_align_percent_rc");
            check_number_in_range("--5end_align_percent_rc", end5_align_percent_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
        }
        if (trim.is_used("--5end_align_identity_rc")) {
            end5_align_identity_rc = trim.get<double>("--5end_align_identity_rc");
            check_number_in_range("--5end_align_identity_rc", end5_align_identity_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
        }
        if (trim.is_used("--3end_len_rc")) {
            end3_len_rc = trim.get<int>("--3end_len_rc");
            check_number_in_range("--3end_len_rc", end3_len_rc, MIN_TARGET, MAX_TARGET, trim, true);
        }
        if (trim.is_used("--3end_align_percent_rc")) {
            end3_align_percent_rc = trim.get<double>("--3end_align_percent_rc");
            check_number_in_range("--3end_align_percent_rc", end3_align_percent_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
        }
        if (trim.is_used("--3end_align_identity_rc")) {
            end3_align_identity_rc = trim.get<double>("--3end_align_identity_rc");
            check_number_in_range("--3end_align_identity_rc", end3_align_identity_rc, MIN_PERCENT, MAX_PERCENT, trim,
                                  false);
        }
        // INIT WORKFLOW
        auto trim_info{barcode_info::get_trim_info()};
        FastqReader fq{input, static_cast<unsigned>(chunk)};
        SequenceInfo sequence_info = !kit.empty() ? trim_info.find(kit)->second : SequenceInfo{forward, reversed};
        sequence_info.update_sequence_info(
            end5_len,
            end5_align_percent,
            end5_align_percent,
            end3_len,
            end3_align_percent,
            end3_align_identity,
            end5_len_rc,
            end5_align_percent_rc,
            end5_align_percent_rc,
            end3_len_rc,
            end3_align_percent_rc,
            end3_align_identity_rc
        );
        std::ofstream out;
        if (output != "-") {
            out.open(output, std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        Work work{fq, static_cast<unsigned>(threads), true, output == "-" ? std::cout : out};
        trim_direction td{myUtility::how_trim(sequence_info)};
        AlignmentConfig align_config{match, mismatch, gap_open, gap_extend};
        std::vector<AlignmentConfig> align_configs;
        for (int i{0}; i < threads; i++) {
            align_configs.push_back(align_config);
        }
        std::fstream logfile{log, std::ios::out};
        if (!logfile) {
            cerr << REDS + "Failed opening log" + COLOR_END << endl;
            exit(1);
        }
        logfile << sequence_info.seq_info() << '\n';
        if (threads == 1) {
            work.run_trim(sequence_info, td, align_configs, logfile);
        } else {
            std::thread t1{&FastqReader::read_chunk_fastq, &fq};
            std::thread t2{
                &Work::run_trim,
                &work,
                std::ref(sequence_info),
                std::ref(td),
                std::ref(align_configs),
                std::ref(logfile)
            };
            t1.join();
            t2.join();
        }
        if (out.is_open()) out.close();
    }
    return 0;
}

#endif //SUBMAIN_H
