#ifndef SUBMAIN_H
#define SUBMAIN_H
#include <thread>
#include <filesystem>
#include <regex>

#include "FastqReader.h"
#include "Work.h"
#include "SequenceInfo.h"
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
#include "Adapter.h"
#include "ArgumentParse.h"

std::pair<std::shared_ptr<SequenceInfo>, std::tuple<int, int, int, int>> parse_trim_arguments(
    argparse::ArgumentParser& parser, bool require_trim)
{
    std::string kit;
    std::string forward;
    std::string reversed;
    int barcode;
    bool kit_used{parser.is_used("--kit")};
    bool barcode_used{parser.is_used("--barcode")};
    bool primers_used{parser.is_used("--primers")};
    if (kit_used && primers_used) {
        std::cerr << "Argument '-p/--primers VAR' not allowed with '-k/--kit VAR'" << std::endl;
        exit(1);
    }
    if (require_trim) {
        if (!kit_used && !primers_used) {
            std::cerr << "One of the arguments '-k/--kit VAR' or '-p/--primers VAR' is required" << std::endl;
            exit(1);
        }
    }
    if (kit_used) {
        kit = parser.get("--kit");
        if (kit.ends_with(".24") || kit.ends_with(".96")) {
            if (!parser.is_used("--barcode")) {
                std::cerr << REDS + "If kit with barcodes used, --barcode must be set" + COLOR_END << std::endl;
                exit(1);
            }
        }
        if (barcode_used) {
            barcode = parser.get<int>("--barcode");
            if (kit.ends_with(".24")) {
                if (barcode < MINB || barcode > MAX24B) {
                    std::cerr << REDS +
                        "If kit with 24 barcodes used, --barcode should be a integer and  in range (1, 24)" +
                        COLOR_END << std::endl;
                    exit(1);
                }
                kit = fmt::format("{}-{}", kit, barcode);
            } else if (kit.ends_with(".96")) {
                if (barcode < MINB || barcode > MAX96B) {
                    std::cerr << REDS +
                        "If kit with 96 barcodes used, --barcode should be a integer and  in range (1, 96)" +
                        COLOR_END << std::endl;
                    exit(1);
                }
                kit = fmt::format("{}-{}", kit, barcode);
            } else {
                std::cerr << WARNS + "If kit with no barcode used, ignore --barcode" + COLOR_END << std::endl;
            }
        }
    } else if (primers_used) {
        std::string primers{parser.get("--primers")};
        if (std::filesystem::exists(std::filesystem::path{primers.data()})) {
            std::ifstream primers_file{primers, std::ios::in};
            std::getline(primers_file, forward);
            std::getline(primers_file, reversed);
            primers_file.close();
        } else {
            auto primers_vec{myutility::split(primers, ",")};
            if (primers_vec.size() != 2) {
                std::cerr << REDS + "if --primer is not file, it should be a pair of primers separated by one comma"
                    + COLOR_END << std::endl;
                exit(1);
            }
            forward = std::string{primers_vec[0]};
            reversed = std::string{primers[1]};
        }
    } else {
        // NO TRIM
    }
    int match{parser.get<int>("--match")};
    check_number_in_range<int>("--match", match, 1, 100, parser, true);
    int mismatch{parser.get<int>("--mismatch")};
    check_number_in_range<int>("--mismatch", mismatch, -100, 0, parser, true);
    int gap_open{parser.get<int>("--gap_open")};
    check_number_in_range<int>("--gap_open", gap_open, -100, 0, parser, true);
    int gap_extend{parser.get<int>("--gap_extend")};
    check_number_in_range<int>("--gap_extend", gap_extend, -100, 0, parser, true);
    int end5_len = DEFAULT_INT, end3_len = DEFAULT_INT, end5_len_rc = DEFAULT_INT, end3_len_rc = DEFAULT_INT;
    float end5_align_percent = DEFAULT_FLOAT, end5_align_identity = DEFAULT_FLOAT;
    float end3_align_percent = DEFAULT_FLOAT, end3_align_identity = DEFAULT_FLOAT;
    float end5_align_percent_rc = DEFAULT_FLOAT, end5_align_identity_rc = DEFAULT_FLOAT;
    float end3_align_percent_rc = DEFAULT_FLOAT, end3_align_identity_rc = DEFAULT_FLOAT;
    if (parser.is_used("--5end_len")) {
        end5_len = parser.get<int>("--5end_len");
        check_number_in_range("--end5_len", end5_len, MIN_TARGET, MAX_TARGET, parser, true);
    }
    if (parser.is_used("--5end_align_percent")) {
        end5_align_percent = parser.get<float>("--5end_align_percent");
        check_number_in_range("--5end_align_percent", end5_align_percent, MIN_PERCENT, MAX_PERCENT, parser, false);
    }
    if (parser.is_used("--5end_align_identity")) {
        end5_align_identity = parser.get<float>("--5end_align_identity");
        check_number_in_range("--5end_align_identity", end5_align_identity, MIN_PERCENT, MAX_PERCENT, parser, false);
    }
    if (parser.is_used("--3end_len")) {
        end3_len = parser.get<int>("--3end_len");
        check_number_in_range("--3end_len", end3_len, MIN_TARGET, MAX_TARGET, parser, true);
    }
    if (parser.is_used("--3end_align_percent")) {
        end3_align_percent = parser.get<float>("--3end_align_percent");
        check_number_in_range("--3end_align_percent", end3_align_percent, MIN_PERCENT, MAX_PERCENT, parser, false);
    }
    if (parser.is_used("--3end_align_identity")) {
        end3_align_identity = parser.get<float>("--3end_align_identity");
        check_number_in_range("--3end_align_identity", end3_align_identity, MIN_PERCENT, MAX_PERCENT, parser, false);
    }
    if (parser.is_used("--5end_len_rc")) {
        end5_len_rc = parser.get<int>("--5end_len_rc");
        check_number_in_range("--5end_len_rc", end5_len_rc, MIN_TARGET, MAX_TARGET, parser, true);
    }
    if (parser.is_used("--5end_align_percent_rc")) {
        end5_align_percent_rc = parser.get<float>("--5end_align_percent_rc");
        check_number_in_range("--5end_align_percent_rc", end5_align_percent_rc, MIN_PERCENT, MAX_PERCENT, parser,
                              false);
    }
    if (parser.is_used("--5end_align_identity_rc")) {
        end5_align_identity_rc = parser.get<float>("--5end_align_identity_rc");
        check_number_in_range("--5end_align_identity_rc", end5_align_identity_rc, MIN_PERCENT, MAX_PERCENT, parser,
                              false);
    }
    if (parser.is_used("--3end_len_rc")) {
        end3_len_rc = parser.get<int>("--3end_len_rc");
        check_number_in_range("--3end_len_rc", end3_len_rc, MIN_TARGET, MAX_TARGET, parser, true);
    }
    if (parser.is_used("--3end_align_percent_rc")) {
        end3_align_percent_rc = parser.get<float>("--3end_align_percent_rc");
        check_number_in_range("--3end_align_percent_rc", end3_align_percent_rc, MIN_PERCENT, MAX_PERCENT, parser,
                              false);
    }
    if (parser.is_used("--3end_align_identity_rc")) {
        end3_align_identity_rc = parser.get<float>("--3end_align_identity_rc");
        check_number_in_range("--3end_align_identity_rc", end3_align_identity_rc, MIN_PERCENT, MAX_PERCENT, parser,
                              false);
    }
    auto parser_info{barcode_info::get_trim_info()};
    SequenceInfo sequence_info{
        !kit_used && !primers_used
            ? parser_info.find("SQK-LSK114")->second
            : !kit.empty()
            ? parser_info.find(kit)->second
            : SequenceInfo{forward, reversed}
    };
    sequence_info.update_sequence_info(
        end5_len,
        end5_align_percent,
        end5_align_identity,
        end3_len,
        end3_align_percent,
        end3_align_identity,
        end5_len_rc,
        end5_align_percent_rc,
        end5_align_identity_rc,
        end3_len_rc,
        end3_align_percent_rc,
        end3_align_identity_rc
    );
    return std::make_pair<std::shared_ptr<SequenceInfo>, std::tuple<int, int, int, int>>(std::make_shared<SequenceInfo>(sequence_info),
    std::make_tuple(match, mismatch, gap_open, gap_extend));
}

int sub_main(int argc, char* argv[])
{
    argparse::ArgumentParser& nanofq{get_arguments(argc, argv)};
    if (nanofq.is_subcommand_used("main")) {
        argparse::ArgumentParser& main{nanofq.at<argparse::ArgumentParser>("main")};
        std::string input{main.get("--input")};
        std::string output{main.get("--output")};
        std::string prefix{main.get("--prefix")};
        bool retain_failed{main.get<bool>("--retain_failed")};
        int n{main.get<int>("--firstN")};
        check_number_in_range("--firstN", n, 1, 1000, main, true);
        std::vector<int> quals;
        if (!main.is_used("--quality")) {
            quals = {25, 20, 18, 15, 12, 10};
        } else {
            quals = {main.get<std::vector<int>>("--quality")};
            for (int i : quals) {
                check_number_in_range<size_t>("--quality", i, 1, 50, main, true);
            }
            std::ranges::sort(quals, std::greater<>());
        }
        std::vector<int> lengths;
        if (main.is_used("--length")) {
            lengths = main.get<std::vector<int>>("--length");
            for (int i : lengths) {
                check_number_in_range<size_t>("--length", i, MINL, MAXL, main, true);
            }
        }
        bool make_plot{main.get<bool>("--plot")};
        std::vector<std::string> format{main.get<std::vector<std::string>>("--format")};
        std::vector<std::string> allowed_choices{"pdf", "jpg", "png"};
        check_choices("--format", format, allowed_choices, main);
        bool plot_mean_length{main.get<bool>("--plot_mean_length")};
        bool plot_n50{main.get<bool>("--plot_n50")};
        bool gc{main.get<bool>("--gc")};
        int min_len{main.get<int>("--min_len")};
        check_number_in_range<size_t>("--min_len", static_cast<size_t>(min_len), MINL, MAXL, main, true);
        int max_len{main.get<int>("--max_len")};
        check_number_in_range<size_t>("--max_len", max_len, MINL, MAXL, main, true);
        float min_quality{main.get<float>("--min_quality")};
        check_number_in_range("--min_quality", min_quality, 0.0f, 100.0f, main, false);
        float min_gc{main.get<float>("--min_gc")};
        float max_gc{main.get<float>("--max_gc")};
        check_number_in_range("--min_gc", min_gc, MIN_PERCENT, MAX_PERCENT, main, false);
        check_number_in_range("--max_gc", max_gc, MIN_PERCENT, MAX_PERCENT, main, false);
        int threads{main.get<int>("--threads")};
        check_number_in_range("--threads", threads, 1, 16, main, true);
        int chunk{main.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, main, true);
        bool do_trim{false};
        if (main.is_used("--kit") || main.is_used("--primers")) {
            do_trim = true;
        }
        // auto [sequence_info, align_config]{parse_trim_arguments(main, do_trim)};
        auto [sequence_info, align_arguments]{parse_trim_arguments(main, do_trim)};
        int max_target_len {std::ranges::max(std::vector<int>{
            std::get<0>(sequence_info->m_top5end),
            std::get<0>(sequence_info->m_top3end),
            std::get<0>(sequence_info->m_bot5end),
            std::get<0>(sequence_info->m_bot3end)
        })};
        int max_query_len { std::ranges::max(std::vector<int>{
            static_cast<int>(sequence_info->m_top5end_query.size()),
            static_cast<int>(sequence_info->m_top5end_query.size()),
            static_cast<int>(sequence_info->m_bot5end_query.size()),
            static_cast<int>(sequence_info->m_bot3end_query.size()),
        })};
        const trim_direction td{myutility::how_trim(*sequence_info)};
        auto [match, mismatch, gap_open, gap_extend] {align_arguments};
        ThreadPool tp{threads};
        std::barrier<> bar{threads};
        auto fqs{myutility::get_fastqs(input)};
        std::vector<read_stats_result> all_stats_result{};
        std::vector<read_stats_result> passed_stats_result{};
        std::vector<AlignmentConfig> align_configs;
        align_configs.reserve(threads);
        for (int i{0}; i < threads; ++i) {
            align_configs.emplace_back(max_target_len, max_query_len, match, mismatch, gap_open, gap_extend);
            // align_configs.push_back(*align_config);
        }
        for (auto& x : align_configs) {
            cout << &x << endl;
        }
        std::filesystem::path prefix_path{prefix};
        std::filesystem::path out_path{output};
        std::filesystem::path all_stats_path{prefix_path.string() + ".raw.stats.tsv"};
        std::filesystem::path passed_stats_path{prefix_path.string() + ".passed.stats.tsv"};
        std::filesystem::path failed_stats_path{prefix_path.string() + ".failed.stats.tsv"};
        std::filesystem::path summary_all_path{prefix_path.string() + ".summary.txt"};
        // std::filesystem::path summary_passed_path{prefix_path.string + ".passed.summary.txt"};
        std::filesystem::path trim_log_path{prefix_path.string() + ".trim_log.txt"};
        std::filesystem::path failed_path{prefix_path.string() + ".failed.fastq"};
        std::ofstream out_ofstream{output, std::ios::out};
        if (!out_ofstream) {
            std::cerr << REDS + "Error happened when opening " + output + COLOR_END << std::endl;
            exit(1);
        }
        std::ofstream all_stats_ofstream{all_stats_path, std::ios::out};
        if (!all_stats_ofstream) {
            std::cerr << REDS + "Error happened when opening " << all_stats_path << COLOR_END << std::endl;
            exit(1);
        }
        std::ofstream passed_stats_ofstream{passed_stats_path, std::ios::out};
        if (!passed_stats_ofstream) {
            std::cerr << REDS + "Error happened when opening " << passed_stats_path << COLOR_END << std::endl;
            exit(1);
        }
        std::ofstream failed_stats_ofstream{failed_stats_path, std::ios::out};
        if (!failed_stats_ofstream) {
            std::cerr << REDS + "Error happened when opening " << failed_stats_path << COLOR_END << std::endl;
            exit(1);
        }
        std::ofstream summary_all_ofstream{summary_all_path, std::ios::out};
        if (!summary_all_ofstream) {
            std::cerr << REDS + "Error happened when opening " << summary_all_path << COLOR_END << std::endl;
            exit(1);
        }
        std::ofstream trim_log_ofstream;
        if (main.is_used("--kit") || main.is_used("--primers")) {
            trim_log_ofstream.open(trim_log_path, std::ios::out);
            if (!trim_log_ofstream) {
                std::cerr << REDS + "Error happened when opening " << summary_all_path << COLOR_END << std::endl;
                exit(1);
            }
            trim_log_ofstream << sequence_info->seq_info() << '\n';
        }
        std::ofstream failed_ofstream;
        if (retain_failed) {
            failed_ofstream.open(failed_path, std::ios::out);
            if (!failed_ofstream) {
                std::cerr << REDS + "Error happened when opening " << summary_all_path << COLOR_END << std::endl;
                exit(1);
            }
        }
        bool is_directory{fqs.has_value()};
        FastqReader fq{input, chunk, is_directory};
        Work work{fq, tp};
        std::mutex all_mtx;
        std::mutex passed_mtx;
        if (is_directory) {
            work.run_main_multi_fqs_in_multi_threads(fqs.value(),
                                                     all_stats_result,
                                                     passed_stats_result,
                                                     gc,
                                                     min_len,
                                                     max_len,
                                                     min_quality,
                                                     min_gc,
                                                     max_gc,
                                                     do_trim,
                                                     *sequence_info,
                                                     td,
                                                     align_configs,
                                                     out_ofstream,
                                                     all_stats_ofstream,
                                                     passed_stats_ofstream,
                                                     failed_stats_ofstream,
                                                     trim_log_ofstream,
                                                     retain_failed,
                                                     failed_ofstream,
                                                     all_mtx,
                                                     passed_mtx);
        } else {
            work.run_main(all_stats_result,
                          passed_stats_result,
                          gc,
                          min_len,
                          max_len,
                          min_quality,
                          min_gc,
                          max_gc,
                          do_trim,
                          *sequence_info,
                          td,
                          align_configs,
                          out_ofstream,
                          all_stats_ofstream,
                          passed_stats_ofstream,
                          failed_stats_ofstream,
                          trim_log_ofstream,
                          retain_failed,
                          failed_ofstream,
                          all_mtx,
                          passed_mtx,
                          bar);
        }
        cout << "all_stats_result: " << all_stats_result.size() << endl;
        cout << "passed_stats_result: " << passed_stats_result.size() << endl;
        auto all_stats_info{work.save_summary(n, quals, lengths, all_stats_result, summary_all_path.c_str(), false)};
        auto passed_stats_info{
            work.save_summary(n, quals, lengths, passed_stats_result, summary_all_path.c_str(), true)
        };

        auto x = all_stats_info;
        auto y = passed_stats_info;
        cout << "raw_mean_length: " << std::get<0>(x) << endl;
        cout << "raw_n50 " << std::get<1>(x) << endl;
        cout << "raw_mean_quality " << std::get<2>(x) << endl;
        cout << "raw_len_std " << std::get<3>(x) << endl;
        cout << "passed_mean_length " << std::get<0>(y) << endl;
        cout << "passed_n50 " << std::get<1>(y) << endl;
        cout << "passed_mean_quality " << std::get<2>(y) << endl;
        cout << "passed_len_std " << std::get<3>(y) << endl;
        // TODO rebuild plot.py
        /*
        auto [all_reads_stats, passed_reads_stats]{get_all_and_passed_read_stats_result(main_stats_result)};
            auto all_summary_info_tuple{work.save_summary(n, quals, lengths, all_reads_stats, summary_all_path.c_str())};
            auto passed_summary_info_tuple{work.save_summary(n, quals, lengths, passed_reads_stats, summary_passed_path.c_str())};
            if (make_plot) {
                auto [all_mean_len, all_n50, all_mean_quality, all_std]{all_summary_info_tuple};
                auto [passed_mean_len, passed_n50, passed_mean_quality, passed_std]{passed_summary_info_tuple};
                std::vector<std::thread> plot_threads;
                auto bin_path{std::string{argv[0]}};
                plot_threads.emplace_back([&all_mean_quality, &format, &work, plot_mean_length, prefix, stats_path, bin_path, all_mean_len, plot_n50, all_n50, all_std]{
                    work.plot(bin_path,
                              stats_path.c_str(),
                              prefix + ".all",
                              plot_mean_length,
                              all_mean_len,
                              plot_n50,
                              all_n50,
                              all_std,
                              format,
                              all_mean_quality);
                });
                plot_threads.emplace_back([&passed_mean_quality, &format, &work, plot_mean_length, prefix, stats_path, bin_path, all_mean_len, plot_n50, all_n50, all_std]{
                    work.plot(bin_path,
                              stats_path.c_str(),
                              prefix + ".passed",
                              plot_mean_length,
                              all_mean_len,
                              plot_n50,
                              all_n50,
                              all_std,
                              format,
                              passed_mean_quality);
                });
                for (auto& t : plot_threads) t.join();
            }
            */
        if (out_ofstream.is_open()) out_ofstream.close();
        if (all_stats_ofstream.is_open()) all_stats_ofstream.close();
        if (passed_stats_ofstream.is_open()) passed_stats_ofstream.close();
        if (summary_all_ofstream.is_open()) summary_all_ofstream.close();
        // if (summary_passed_ofstream.is_open()) summary_passed_ofstream.close();
        if (trim_log_ofstream.is_open()) trim_log_ofstream.close();
        if (failed_ofstream.is_open()) failed_ofstream.close();
    } else if (nanofq.is_subcommand_used("stats")) {
        argparse::ArgumentParser& stats{nanofq.at<argparse::ArgumentParser>("stats")};
        std::string input{stats.get("--input")};
        std::string output{stats.get("--output")}; // output should not be stdout
        std::string summary{stats.get("--summary")};
        int n{stats.get<int>("--firstN")};
        check_number_in_range("--firstN", n, 1, 1000, stats, true);
        std::vector<int> quals;
        if (!stats.is_used("--quality")) {
            quals = {25, 20, 18, 15, 12, 10};
        } else {
            quals = {stats.get<std::vector<int>>("--quality")};
            for (int i : quals) {
                check_number_in_range("--quality", i, 1, 50, stats, true);
            }
            std::ranges::sort(quals, std::greater<>());
        }
        std::vector<int> lengths;
        if (stats.is_used("--length")) {
            lengths = stats.get<std::vector<int>>("--length");
            for (int i : lengths) {
                check_number_in_range("--length", i, MINL, MAXL, stats, true);
            }
            std::ranges::sort(lengths, std::greater<>());
        }
        bool make_plot{false};
        bool plot_mean_length{false};
        bool plot_n50{false};
        std::string plot_prefix;
        if (stats.is_used("--plot")) {
            make_plot = true;
            plot_prefix = stats.get("--plot");
            plot_mean_length = stats.get<bool>("--plot_mean_length");
            plot_n50 = stats.get<bool>("--plot_n50");
        }
        int threads{stats.get<int>("--threads")};
        check_number_in_range("--threads", threads, 1, 16, stats, true);
        int chunk{stats.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, stats, true);
        bool gc{stats.get<bool>("--gc")};
        std::vector<std::string> format{stats.get<std::vector<std::string>>("--format")};
        std::vector<std::string> allowed_choices{"pdf", "jpg", "png"};
        check_choices<std::string>("--format", format, allowed_choices, stats);
        std::ofstream out;
        if (output != "-") {
            out.open(output.data(), std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        ThreadPool tp{threads};
        std::barrier<> bar{threads};
        std::vector<read_stats_result> stats_result{};
        // std::vector<main_read_stats_result> main_stats_result{};
        auto fqs{myutility::get_fastqs(input)};
        FastqReader fq{input, chunk, fqs.has_value()};
        Work work{fq, tp};
        if (fqs.has_value()) {
            work.run_stats_multi_fqs_in_multi_threads(fqs.value(), stats_result, out, gc);
        } else {
            work.run_stats(stats_result, output != "-" ? out : std::cout, gc );
        }
        auto summary_info_tuple{
            work.save_summary(
                n, quals, lengths, stats_result, summary, false)
        };
        if (make_plot) {
            auto [mean_len, n50, mean_quality, std] = summary_info_tuple;
            work.plot(std::string{argv[0]},
                      output,
                      plot_prefix,
                      plot_mean_length,
                      mean_len,
                      plot_n50,
                      n50,
                      std,
                      format,
                      mean_quality);
        }
        if (out.is_open()) { out.close(); }
    } else if (nanofq.is_subcommand_used("filter")) {
        argparse::ArgumentParser& filter{nanofq.at<argparse::ArgumentParser>("filter")};
        std::string input{filter.get("--input")};
        std::string output{filter.get("--output")};
        int min_length{filter.get<int>("--min_len")};
        check_number_in_range("--min_len", min_length, MINL, MAXL, filter, true);
        int max_length{filter.get<int>("--max_len")};
        check_number_in_range("--max_len", max_length, MINL, MAXL, filter, true);
        float min_quality{filter.get<float>("--min_quality")};
        check_number_in_range("--min_quality", min_quality, 0.0f, 100.0f, filter, false);
        bool gc{filter.get<bool>("--gc")};
        float min_gc{filter.get<float>("--min_gc")};
        check_number_in_range("--min_gc", min_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
        float max_gc{filter.get<float>("--max_gc")};
        check_number_in_range("--max_gc", max_gc, MIN_PERCENT, MAX_PERCENT, filter, false);
        int threads{filter.get<int>("--threads")};
        check_number_in_range("--threads", threads, 1, 16, filter, true);
        int chunk{filter.get<int>("--chunk")};
        check_number_in_range("--chunk", chunk, MINC, MAXC, filter, true);
        std::ofstream out;
        if (output != "-") {
            out = std::ofstream{output.data(), std::ios::out};
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        ThreadPool tp{threads};
        std::barrier<> bar{threads};
        auto fqs{myutility::get_fastqs(input)};
        FastqReader fq{input, chunk, fqs.has_value()};
        Work work{fq, tp};
        if (fqs.has_value()) {
            work.run_filter_multi_fqs_in_multi_threads(
                fqs.value(),
                gc,
                min_length,
                max_length,
                min_quality,
                min_gc,
                max_gc,
                output != "-" ? out : std::cout);
        } else {
            std::atomic<size_t> counter{0};
            work.run_filter(counter,
                            gc,
                            min_length,
                            max_length,
                            min_quality,
                            min_gc,
                            max_gc,
                            output != "-" ? out : std::cout);
        }
        if (out.is_open()) out.close();
    } else if (nanofq.is_subcommand_used("index")) {
        argparse::ArgumentParser& index{nanofq.at<argparse::ArgumentParser>("index")};
        std::string input{index.get("input")};
        FastqReader fq{input, 20000, false};
        ThreadPool tp{1};
        Work work{fq, tp};
        work.run_index(true);
    } else if (nanofq.is_subcommand_used("find")) {
        argparse::ArgumentParser& find{nanofq.at<argparse::ArgumentParser>("find")};
        std::string input{find.get("--input")};
        std::string output{find.get("--output")};
        std::string reads{find.get("--reads")};
        bool use_index{find.get<bool>("--use_index")};
        FastqReader fq{input, 5000, false};
        std::ofstream out;
        if (output != "-") {
            out.open(output, std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        ThreadPool tp{1};
        Work work{fq, tp};
        work.run_find(reads, output != "-" ? out : std::cout, use_index);
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

        std::ofstream out;
        if (output != "-") {
            out.open(output, std::ios::out);
            if (!out) {
                std::cerr << REDS + "Failed when opened " + output << COLOR_END << std::endl;
                exit(1);
            }
        }
        auto [sequence_info, align_arguments]{parse_trim_arguments(trim, true)};
        int max_target_len {std::ranges::max(std::vector<int>{
            std::get<0>(sequence_info->m_top5end),
            std::get<0>(sequence_info->m_top3end),
            std::get<0>(sequence_info->m_bot5end),
            std::get<0>(sequence_info->m_bot3end)
        })};
        int max_query_len { std::ranges::max(std::vector<int>{
            static_cast<int>(sequence_info->m_top5end_query.size()),
            static_cast<int>(sequence_info->m_top5end_query.size()),
            static_cast<int>(sequence_info->m_bot5end_query.size()),
            static_cast<int>(sequence_info->m_bot3end_query.size()),
        })};
        ThreadPool tp{threads};
        trim_direction td{myutility::how_trim(*sequence_info)};
        auto [match, mismatch, gap_open, gap_extend] {align_arguments};
        std::vector<AlignmentConfig> align_configs;
        align_configs.reserve(threads);
        for (int i{0}; i < threads; i++) {
            align_configs.emplace_back(max_target_len, max_query_len, match, mismatch, gap_open, gap_extend);
        }
        std::fstream logfile{log, std::ios::out};
        if (!logfile) {
            cerr << REDS + "Failed opening log" + COLOR_END << endl;
            exit(1);
        }
        logfile << sequence_info->seq_info() << '\n';
        auto fqs = myutility::get_fastqs(input);
        FastqReader fq{input, chunk, fqs.has_value()};
        Work work{fq, tp};
        std::barrier<> bar{threads};
        if (fqs.has_value()) {
            work.run_trim_multi_fqs_in_multi_threads(
                fqs.value(),
                *sequence_info,
                td,
                align_configs,
                logfile,
                output != "-" ? out : std::cout);
        } else {
            std::atomic<size_t> counter;
            work.run_trim(counter, *sequence_info, td, align_configs, logfile, output != "-" ? out : std::cout, bar);
        }
        if (out.is_open()) out.close();
    } else if (nanofq.is_subcommand_used("compress")) {
        argparse::ArgumentParser& compress{nanofq.at<argparse::ArgumentParser>("compress")};
        std::string input{compress.get("input")};
        std::string output{compress.get("output")};
        int reads_number{compress.get<int>("--number")};
        check_number_in_range("--number", reads_number, 5, 50, compress, true);
        nanobgzip::nano_compress(input, output, fmt::format("{}.index", output), reads_number);
    }
    return 0;
}

#endif //SUBMAIN_H
