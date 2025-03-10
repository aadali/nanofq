#include "Work.h"

#include <limits>
#include <thread>
#include <tuple>
#include <vector>
#include <fmt/core.h>

using std::cout;
using std::endl;

Work::Work(
    FastqReader& fq,
    ThreadPool& threads_pool) :
    m_fq(fq),
    m_threads_pool(threads_pool) {}

std::vector<std::pair<unsigned, unsigned>> Work::get_edges(int size) const
{
    std::vector<std::pair<unsigned, unsigned>> idx_ranges;
    const int step{size / m_threads_pool.threads_number()};
    int start{0}, stop{0};
    for (int i{0}; i < size; i += step) {
        start = i;
        stop = step + start;
        if (stop + step > size) {
            stop = size;
            idx_ranges.emplace_back(start, stop);
            break;
        }
        idx_ranges.emplace_back(start, stop);
    }
    return idx_ranges;
}

void Work::run_stats(
    std::vector<read_stats_result>& stats_result_vec,
    std::ostream& out,
    bool gc)
{
    if (m_threads_pool.threads_number() == 1) {
        if (!stats_result_vec.empty()) {
            cerr << "When used one thread, the parameter stats_result must be empty" << endl;
            exit(1);
        }
        while (true) {
            Read read{m_fq.read_one_fastq()}; // get one fastq record from FastqReader
            if (*read.get_id() == finished_read_name) {
                // if get the last read, return.
                // The readname of last record is "FINISHED FINISHED FINISHED"
                return;
            }
            stats_one_thread(read, stats_result_vec, out, gc);
        }
    }
    size_t total_reads{0};
    while (true) {
        auto reads_ptr{m_fq.read_chunk_fastq()};
        total_reads += reads_ptr->size();
        auto bins = get_edges(reads_ptr->size());
        for (auto [start, end] : bins) {
            m_threads_pool.enqueue([this, start, end, gc, &stats_result_vec, &out, reads_ptr ](){
                stats(start, end, reads_ptr, stats_result_vec, out, gc);
            });
        }
        if (m_fq.read_finish())
            break;
    }
    // wait for all reads statsed finished
    while (total_reads != stats_result_vec.size()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

void Work::run_main(
    std::vector<read_stats_result>& all_reads_stats_result,
    std::vector<read_stats_result>& passed_reads_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    std::vector<AlignmentConfig>& align_configs,
    std::ofstream& out_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    std::ofstream& trim_log_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    std::mutex& all_mtx,
    std::mutex& passed_mtx,
    std::barrier<>& bar
)
{
    if (m_threads_pool.threads_number() == 1) {
        if (!all_reads_stats_result.empty() || !passed_reads_stats_result.empty()) {
            std::cerr << REDS <<
                "When used one thread, the parameter all_reads_stats_result/passed_reads_stats_result must be empty" <<
                COLOR_END <<
                std::endl;
            exit(1);
        }
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (*read.get_id() == finished_read_name) {
                return;
            }
            main_one_thread(
                read,
                all_reads_stats_result,
                passed_reads_stats_result,
                gc,
                min_len,
                max_len,
                min_quality,
                min_gc,
                max_gc,
                do_trim,
                seq_info,
                td,
                align_configs[0],
                out_ofstream,
                all_stats_ofstream,
                passed_stats_ofstream,
                failed_stats_ofstream,
                trim_log_ofstream,
                retain_failed,
                failed_ofstream,
                all_mtx,
                passed_mtx
            );
        }
    }
    size_t total_reads{0};
    while (true) {
        auto reads_ptr{m_fq.read_chunk_fastq()};
        total_reads += reads_ptr->size();
        auto bins{get_edges(reads_ptr->size())};
        for (int i{0}; i < bins.size(); ++i) {
            auto [start, end] = bins[i];
            m_threads_pool.enqueue(
                [i, this, start, end, reads_ptr, &all_reads_stats_result, &passed_reads_stats_result, gc, min_len,
                    max_len, min_quality, min_gc,
                    max_gc, do_trim, &seq_info, &td, &align_configs, &out_ofstream, &trim_log_ofstream, &
                    all_stats_ofstream, &passed_stats_ofstream, &failed_stats_ofstream, retain_failed, &failed_ofstream,
                    &all_mtx, &passed_mtx, &bar]{
                    this->main(start,
                               end,
                               reads_ptr,
                               all_reads_stats_result,
                               passed_reads_stats_result,
                               gc,
                               min_len,
                               max_len,
                               min_quality,
                               min_gc,
                               max_gc,
                               do_trim,
                               seq_info,
                               td,
                               align_configs[i],
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
            );
        }
        if (m_fq.read_finish()) break;
    }
    while (total_reads != all_reads_stats_result.size()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

void Work::run_main_multi_fqs_in_multi_threads(
    const std::vector<std::filesystem::path>& paths,
    std::vector<read_stats_result>& all_reads_stats_result,
    std::vector<read_stats_result>& passed_reads_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    std::vector<AlignmentConfig>& align_configs,
    std::ofstream& out_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    std::ofstream& trim_log_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    std::mutex& all_mtx,
    std::mutex& passed_mtx
)
{
    auto bins{get_edges(paths.size())};
    cout << "fastq numbers: " << paths.size() << endl;
    std::vector<std::thread> threads;
    for (int i{0}; i < bins.size(); ++i) {
        auto [start, end]{bins[i]};
        cout << "start: " << start << "; end: " << end << std::endl;
        threads.emplace_back(
            [this, i, &paths, start, end, &all_reads_stats_result, &passed_reads_stats_result,gc, min_len, max_len,
                min_quality, min_gc, max_gc, do_trim,
                &seq_info, &td, &align_configs, &out_ofstream, &trim_log_ofstream, &all_stats_ofstream, &
                passed_stats_ofstream, &failed_stats_ofstream, retain_failed, &
                failed_ofstream, &all_mtx, &passed_mtx
            ]{
                main_multi_fqs_in_one_thread(
                    paths, start, end, all_reads_stats_result, passed_reads_stats_result, gc, min_len, max_len,
                    min_quality,
                    min_gc, max_gc, do_trim,
                    seq_info, td, align_configs[i], out_ofstream, trim_log_ofstream,
                    all_stats_ofstream, passed_stats_ofstream, failed_stats_ofstream, retain_failed, failed_ofstream,
                    all_mtx, passed_mtx);
            }
        );
    }
    for (std::thread& t : threads) {
        t.join();
    }
}

void Work::run_stats_multi_fqs_in_multi_threads(
    const std::vector<std::filesystem::path>& paths,
    std::vector<read_stats_result>& stats_result,
    std::ostream& out,
    bool gc)
{
    auto bins{get_edges(paths.size())};
    std::vector<std::thread> threads;
    for (auto [start, end] : bins) {
        threads.emplace_back([start, end, gc, &stats_result, &paths, &out, this]{
            stats_multi_fqs_in_one_thread(paths, start, end, &stats_result, out, gc);
        });
    }
    for (std::thread& t : threads) {
        t.join();
    }
}


void Work::run_filter(
    std::atomic<size_t>& counter,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    std::ostream& out) const
{
    if (m_threads_pool.threads_number() == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (*read.get_id() == finished_read_name) {
                return;
            }
            filter_one_thread(read, gc, min_len, max_len, min_quality, min_gc, max_gc, out);
        }
    }

    size_t total_reads{0};
    while (true) {
        auto reads_ptr{m_fq.read_chunk_fastq()};
        total_reads += reads_ptr->size();
        auto bins = get_edges(reads_ptr->size());
        for (auto [start, end] : bins) {
            m_threads_pool.enqueue(
                [this, start, end, reads_ptr, gc, min_len, max_len, min_quality, min_gc, max_gc, &out, &counter](){
                    filter(start, end, reads_ptr, counter, gc, min_len, max_len, min_quality, min_gc, max_gc, out);
                });
        }
        if (m_fq.read_finish()) {
            break;
        }
    }
    while (total_reads != counter) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

void Work::run_filter_multi_fqs_in_multi_threads(
    const std::vector<std::filesystem::path>& paths,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    std::ostream& out) const
{
    auto bins{get_edges(paths.size())};
    std::vector<std::thread> threads;
    for (auto [start, end] : bins) {
        threads.emplace_back([this, &paths, start, end, gc, min_len, max_len, min_quality, min_gc, max_gc, &out]{
            filter_multi_fqs_in_one_thread(paths, start, end, gc, min_len, max_len, min_quality, min_gc, max_gc, out);
        });
    }
    for (std::thread& t : threads) {
        t.join();
    }
}

void Work::run_find(
    const std::string& input_reads,
    std::ostream& out,
    bool use_index) const
{
    m_fq.find(input_reads, out, use_index);
}


std::tuple<float, int, float, float> Work::save_summary(
    int n,
    const std::vector<int>& read_quals,
    const std::vector<int>& read_lengths,
    std::vector<read_stats_result>& reads_stats_result,
    const std::string& summary_file_path,
    bool is_passed
)
{
    std::stringstream summary_content;
    auto summary_info{summary_stats_result(n, read_quals, read_lengths, reads_stats_result)};
    summary_content << std::get<0>(summary_info);
    std::ofstream summary_file;
    if (!is_passed) {
        summary_file.open(summary_file_path, std::ios::out);
    } else {
        summary_file.open(summary_file_path, std::ios::app);
    }
    if (summary_file) {
        if (is_passed) summary_file << "### The following is information about passed reads\n";
        summary_file << summary_content.str();
        summary_file.close();
    } else {
        std::cerr << WARNS + "Failed when opening " + summary_file_path + ". Try write the summary into ~/SumMarY.txt" +
            COLOR_END
            << std::endl;
        std::ofstream try_summary_file{"./SumMarY.txt", std::ios::out};
        try {
            try_summary_file << summary_content.str();
            try_summary_file.close();
        }
        catch (const std::exception& e) {
            std::cerr << REDS << "Failed when trying write summary into ~/SumMarY.txt" + COLOR_END << std::endl;
            exit(1);
        }
    }
    auto res{
        std::make_tuple(
            std::get<1>(summary_info),
            std::get<2>(summary_info),
            std::get<3>(summary_info),
            std::get<4>(summary_info)
        )
    };
    return res; // mean_read_len, read_len_n50, mean_read_quality, read_len_std
}

void Work::plot(
    const std::string& argv0,
    const std::string& input,
    const std::string& prefix,
    bool plot_mean_length,
    float mean_length,
    bool plot_n50,
    int n50,
    float std,
    const std::vector<std::string>& fmt,
    float mean_quality)
{
    auto bin_dir = std::filesystem::path(argv0.data()).parent_path();
    auto script = std::string{bin_dir.append("plot.py")};
    std::string cmd{fmt::format("python3 {} -i {} -p {} ", script, input, prefix)};
    if (plot_mean_length) {
        cmd.append(fmt::format("--plot_mean_length -M {} ", mean_length));
    }

    if (plot_n50) {
        cmd.append(fmt::format("--plot_n50 -N {} ", n50));
    }

    for (const std::string& f : fmt) {
        cmd.append(fmt::format("-f {} ", f));
    }
    cmd.append(fmt::format("-Q {} -s {} ", mean_quality, std));
    // cout << cmd << endl;
    // int res = std::system(cmd.data());
    if (int res{std::system(cmd.data())}; res != 0) {
        std::cerr << WARNS << "Error found when making plot, you can plot using " << script << " manually" << COLOR_END
            << std::endl;
    }
}

void Work::run_index(bool force_index) const
{
    m_fq.index(force_index);
}

void Work::run_trim(
    std::atomic<size_t>& counter,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    std::vector<AlignmentConfig>& align_configs,
    std::ostream& log_fstream,
    std::ostream& out,
    std::barrier<>& bar) const
{
    if (m_threads_pool.threads_number() == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (*read.get_id() == finished_read_name) {
                return;
            }
            trim_one_thread(read, seq_info, td, align_configs[0], log_fstream, out);
        }
    }
    size_t total_reads{0};
    int a{0};
    while (true) {
        auto reads_ptr{m_fq.read_chunk_fastq()};
        total_reads += reads_ptr->size();
        auto bins = get_edges(reads_ptr->size());
        for (int i{0}; i < bins.size(); ++i) {
            ++a;
            auto [start, end] = bins[i];
            m_threads_pool.enqueue(
                [i, this, start, end, reads_ptr, &counter, &seq_info, &td, &align_configs,&log_fstream, &out, &bar]{
                    trim(start, end, reads_ptr, counter, seq_info, td, align_configs[i], log_fstream, out, bar);
                });
        }
        if (m_fq.read_finish()) {
            break;
        }
    }
    while (total_reads != counter) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

void Work::run_trim_multi_fqs_in_multi_threads(
    const std::vector<std::filesystem::path>& paths,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    std::vector<AlignmentConfig>& align_configs,
    std::ostream& log_fstream,
    std::ostream& out) const
{
    auto bins{get_edges(paths.size())};
    std::vector<std::thread> threads;
    for (int i{0}; i < bins.size(); ++i) {
        auto [start, end] = bins[i];
        threads.emplace_back([this, start, end, &paths, &seq_info, &td, &align_configs, &log_fstream, &out, i]{
            trim_multi_fqs_in_one_thread(paths, start, end, seq_info, td, align_configs[i], log_fstream, out);
        });
    }
    for (std::thread& t : threads) {
        t.join();
    }
}

void Work::stats(
    int start,
    int end,
    std::shared_ptr<std::vector<Read>> reads_ptr,
    std::vector<read_stats_result>& read_stats,
    std::ostream& out,
    bool gc)
{
    for (int idx{start}; idx < end; idx++) {
        const Read& read = (*reads_ptr)[idx];
        auto name{read.get_id()};
        unsigned length{read.get_length()};
        float quality{read.calculate_read_quality()};
        float gc_content{gc ? read.get_gc_content() : 0.0f};
        std::osyncstream{out} << fmt::format("{}\t{}\t{}\t{}\n", *name, length, quality, gc_content);
        {
            std::unique_lock lock{m_mtx};
            read_stats.emplace_back(name, length, quality, gc_content);
        }
    }
}

std::tuple<std::string, float, int, float, float> Work::summary_stats_result(
    int n,
    const std::vector<int>& read_quals,
    const std::vector<int>& read_lengths,
    std::vector<read_stats_result>& stats_result)
// return std::tuple<read_name, mean_read_len, read_len_n50, mean_read_quality, read_len_std>
{
    std::tuple<std::string, float, int, float, float> summary_tuple;
    std::stringstream summary_stream;
    ulong total_bases_number{
        std::accumulate(stats_result.begin(),
                        stats_result.end(),
                        0ul,
                        [](ulong x, const read_stats_result& y){ return x + static_cast<ulong>(std::get<1>(y)); })
    };
    ulong total_reads_number{stats_result.size()};
    summary_stream << fmt::format("BasesNumber\t{:.6f}Mb\n", static_cast<double>(total_bases_number) / 1000000);
    summary_stream << fmt::format("ReadsNumber\t{}\n", total_reads_number);
    // decreased by read length
    std::ranges::sort(stats_result, [](auto& x, auto& y){ return std::get<1>(x) > std::get<1>(y); });
    size_t n_idx{0};
    double n_length;
    auto get_n_percent{
        [&n_idx, &n_length, &stats_result, total_bases_number](double percent){
            for (; n_idx < stats_result.size(); n_idx++) {
                ulong this_len{get<1>(stats_result[n_idx])};
                n_length += this_len;
                if (n_length / total_bases_number > percent) {
                    return this_len;
                }
            }
            exit(1);
        }
    };
    ulong n10{get_n_percent(0.1)};
    ulong n50{get_n_percent(0.5)};
    ulong n90{get_n_percent(0.9)};
    summary_stream << fmt::format("N10\t{}\nN50\t{}\nN90\t{}\n", n10, n50, n90);

    auto get_len_quantile{
        [&](double quantile){
            size_t quantile_idx{static_cast<size_t>(total_reads_number * quantile)};
            return std::get<1>(stats_result[total_reads_number - quantile_idx + 1]);
        }
    };
    auto get_quantile_quality{
        [&](double quantile){
            size_t quantile_idx{static_cast<size_t>(total_reads_number * quantile)};
            return std::get<2>(stats_result[total_reads_number - quantile_idx + 1]);
        }
    };

    unsigned read_len_quantile25{get_len_quantile(0.25)};
    unsigned read_len_quantile50{get_len_quantile(0.5)};
    unsigned read_len_quantile75{get_len_quantile(0.75)};
    double mean_read_len{static_cast<double>(total_bases_number) / static_cast<double>(total_reads_number)};
    summary_stream << fmt::format("ReadLenQuantile25\t{}\nReadLenQuantile50\t{}\nReadLenQuantile75\t{}\n",
                                  read_len_quantile25, read_len_quantile50, read_len_quantile75);
    summary_stream << fmt::format("ReadMeanLen\t{:.2f}\n", mean_read_len);
    double sum_std = std::accumulate(stats_result.begin(),
                                     stats_result.end(),
                                     0.0,
                                     [&](double x, const read_stats_result& y){
                                         return x + std::pow(static_cast<double>(std::get<1>(y)) - mean_read_len, 2);
                                     });
    double std{std::sqrt(sum_std / total_reads_number)};
    summary_stream << fmt::format("ReadLenStd\t{:.2f}\n", std);
    std::stringstream longest_reads;
    longest_reads << fmt::format("#Top {} longest  reads\nnth\tReadName\tReadLen\tReadQuality\tGC\n", n);
    for (int i{0}; i < n; i++) {
        if (i >= total_reads_number) {
            break;
        }
        longest_reads << fmt::format(
            "{}\t{}\t{}\t{}\t{}\n",
            i + 1,
            *std::get<0>(stats_result[i]),
            std::get<1>(stats_result[i]),
            fmt::format("{:.2f}", std::get<2>(stats_result[i])),
            fmt::format("{:.2f}", std::get<3>(stats_result[i])));
    }
    std::stringstream lengths_info;
    int read_idx{0};
    ulong reads_count = 0;
    ulong bases_count = 0;
    auto stats_depend_length{
        [&](int length){
            for (; read_idx < total_reads_number; read_idx++) {
                if (std::get<1>(stats_result[read_idx]) < length || read_idx == total_reads_number - 1) {
                    lengths_info << fmt::format(
                        "ReadLength > {}\t{}({});{}({})\n",
                        length,
                        reads_count,
                        fmt::format(
                            "{:.2f}%",
                            100 * static_cast<double>(reads_count) / total_reads_number),
                        fmt::format("{:.6f}Mb", static_cast<double>(bases_count) / 1000000),
                        fmt::format(
                            "{:.2f}%",
                            100 * static_cast<double>(bases_count) / total_bases_number));
                    break;
                }
                reads_count += 1;
                bases_count += std::get<1>(stats_result[read_idx]);
            }
        }
    };
    if (!read_lengths.empty()) {
        for (const int& length : read_lengths) {
            stats_depend_length(length);
        }
    }
    // decreased by read quality
    std::ranges::sort(stats_result, [](auto& x, auto& y){ return std::get<2>(x) > std::get<2>(y); });
    float total_error_probability{
        std::accumulate(stats_result.begin(),
                        stats_result.end(),
                        0.0f,
                        [](float x, const read_stats_result& y){ return x + std::pow(10.0f, std::get<2>(y) / -10.0f); })
    };
    double mean_quality{std::log10(total_error_probability / total_reads_number) * -10};
    double read_quality_quantile25{get_quantile_quality(0.25f)};
    double read_quality_quantile50{get_quantile_quality(0.5f)};
    double read_quality_quantile75{get_quantile_quality(0.75f)};
    summary_stream << fmt::format("ReadQualityQuantile25\t{}\nReadQualityQuantile50\t{}\nReadQualityQuantile75\t{}\n",
                                  fmt::format("{:.2f}", read_quality_quantile25),
                                  fmt::format("{:.2f}", read_quality_quantile50),
                                  fmt::format("{:.2f}", read_quality_quantile75));
    summary_stream << fmt::format("ReadMeanQuality\t{}", fmt::format("{:.2f}\n", mean_quality));

    summary_stream << "#ReadQuality > SpecifiedValue\tReadsNumber(ReadsPercent);BasesNumber(BasesPercent)\n";
    read_idx = 0;
    reads_count = 0;
    bases_count = 0;

    auto stats_depend_quality{
        [&](int quality){
            if (read_idx >= total_reads_number - 1) {
                summary_stream << fmt::format(
                    "ReadQuality > {}\t{}({});{}({})\n",
                    quality,
                    total_reads_number,
                    fmt::format(
                        "{:.2f}%",
                        100 * static_cast<double>(total_reads_number) / total_reads_number),
                    fmt::format("{:.6f}Mb", static_cast<double>(total_bases_number) / 1000000),
                    fmt::format(
                        "{:.2f}%",
                        100 * static_cast<double>(total_bases_number) / total_bases_number));
            }
            for (; read_idx < total_reads_number; read_idx++) {
                double this_read_quality{std::get<2>(stats_result[read_idx])};
                // stats_result must be decreased by read quality
                if (this_read_quality < quality) {
                    summary_stream << fmt::format(
                        "ReadQuality > {}\t{}({});{}({})\n",
                        quality,
                        reads_count,
                        fmt::format(
                            "{:.2f}%",
                            100 * static_cast<double>(reads_count) / total_reads_number),
                        fmt::format("{:.6f}Mb", static_cast<double>(bases_count) / 1000000),
                        fmt::format(
                            "{:.2f}%",
                            100 * static_cast<double>(bases_count) / total_bases_number));
                    break;
                }
                ++reads_count;
                bases_count += std::get<1>(stats_result[read_idx]);
                if (read_idx == total_reads_number - 1){
                    summary_stream << fmt::format(
                         "ReadQuality > {}\t{}({});{}({})\n",
                         quality,
                         reads_count,
                         fmt::format(
                             "{:.2f}%",
                             100 * static_cast<double>(reads_count) / total_reads_number),
                         fmt::format("{:.6f}Mb", static_cast<double>(bases_count) / 1000000),
                         fmt::format(
                             "{:.2f}%",
                             100 * static_cast<double>(bases_count) / total_bases_number));
                }
            }
        }
    };
    for (const int& quality : read_quals) {
        std::cout << "calculate: "<< quality << std::endl;
        stats_depend_quality(quality);
    }
    summary_stream << lengths_info.str();
    summary_stream << longest_reads.str();
    summary_stream << fmt::format("#Top {} high quality reads\nnth\tReadName\tReadLen\tReadQuality\tGC\n", n);
    for (int i{0}; i < n; i++) {
        if (i >= total_reads_number) {
            break;
        }
        summary_stream << fmt::format("{}\t{}\t{}\t{}\t{}\n",
                                      i + 1,
                                      *std::get<0>(stats_result[i]),
                                      std::get<1>(stats_result[i]),
                                      fmt::format("{:.2f}", std::get<2>(stats_result[i])),
                                      fmt::format("{:.2f}", std::get<3>(stats_result[i])));
    }
    summary_tuple = std::make_tuple(summary_stream.str(), mean_read_len, n50, mean_quality, std);
    return summary_tuple;
}

// std::tuple<std::string, float, int, float, float> Work::main_summary_stats_result(
//     int n,
//     const std::vector<int>& read_quals,
//     const std::vector<int>& read_lengths,
//     std::vector<read_stats_result>& stats_result_vec
// )
// {
//     std::vector<read_stats_result> stats_result;
//     for (auto& element: stats_result_vec){
//         stats_result.emplace_back(
//             std::get<0>(element),
//             std::get<1>(element),
//             std::get<2>(element),
//             std::get<3>(element)
//         );
//     }
//     auto summary_tuple{summary_stats_result(n, reads,)}
//     if (get_all) {
//         for (auto& element : stats_result_vec) {
//             stats_result.emplace_back(std::get<0>(element),
//                                       std::get<1>(element),
//                                       std::get<2>(element),
//                                       std::get<3>(element));
//         }
//         auto summary_tuple{summary_stats_result(n, read_quals, read_lengths, stats_result)};
//         return summary_tuple;
//     }
//     for (auto& element : stats_result_vec) {
//         if (std::get<7>(element)) {
//             stats_result.emplace_back(
//                 std::get<0>(element),
//                 std::get<4>(element),
//                 std::get<5>(element),
//                 std::get<6>(element)
//             );
//         }
//     }
//     auto summary_tuple{summary_stats_result(n, read_quals, read_lengths, stats_result)};
//     return summary_tuple;
// }

/**
 * @brief stats a Read object, store the stats result into a vector and output stats information for further using
 * @param read a Rread object
 * @param read_stats_vec vector to store stats result of each read
 * @param out ostream or ofstream to output stats result for further use
 * @return void
 */
void Work::stats_one_thread(
    const Read& read,
    std::vector<read_stats_result>& read_stats_vec,
    std::ostream& out,
    bool gc)
{
    unsigned len{read.get_length()};
    float quality{read.calculate_read_quality()};
    float gc_content{gc ? read.get_gc_content() : 0.0f};
    read_stats_vec.emplace_back(read.get_id(), len, quality, gc_content);
    auto line = fmt::format("{}\t{}\t{}\t{}\n",
                            *read.get_id(),
                            len,
                            fmt::format("{:.{}f}", quality, 2),
                            fmt::format("{:.{}f}", gc_content, 2));
    out << line;
}

/**
 * @brief Processes and stats multiple FASTQ files in a single thread, storing the results and outputting them for further use
 * read and stats FASTQ file one by one
 * @param paths A vector of FASTQ file paths to be processed
 * @param start the start index of paths
 * @param end  the end index of paths
 * @param stats_result pointer to a vector to store the stats result of each read from the FASTQ files
 * @param out An ostream or ofstream to output the stats result for further use
 * @return void
 */
void Work::stats_multi_fqs_in_one_thread(
    const std::vector<std::filesystem::path>& paths,
    size_t start,
    size_t end,
    std::vector<read_stats_result>* stats_result,
    std::ostream& out,
    bool gc)
{
    for (size_t i{start}; i < end; i++) {
        int l;
        gzFile file = gzopen(paths[i].c_str(), "rb");
        kseq_t* seq = kseq_init(file);
        while (true) {
            l = kseq_read(seq);
            Read read{FastqReader::fastq_record_ok(l, seq, paths[i].c_str())};
            if (*read.get_id() == finished_read_name) {
                break;
            }
            auto shared_name{read.get_id()};
            unsigned length{read.get_length()};
            float quality{read.calculate_read_quality()};
            float gc_content{gc ? read.get_gc_content() : 0.0f};
            std::osyncstream{out} << fmt::format("{}\t{}\t{}\t{}\n", *shared_name, length, quality, gc_content);
            {
                std::unique_lock lock{m_mtx};
                stats_result->emplace_back(shared_name, length, quality, gc_content);
            }
        }
        kseq_destroy(seq);
        gzclose(file);
    }
}

void Work::filter(
    int start,
    int end,
    std::shared_ptr<std::vector<Read>> reads_ptr,
    std::atomic<size_t>& counter,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    std::ostream& out)
{
    for (int idx{start}; idx < end; idx++) {
        const Read& read = (*reads_ptr)[idx];
        unsigned len{read.get_length()};
        float quality{read.calculate_read_quality()};
        if (gc) {
            if (float gc_content{read.get_gc_content()};
                len >= min_len && len <= max_len && quality > min_quality &&
                gc_content > min_gc && gc_content < max_gc) {
                std::osyncstream{out} << read.get_record();
            }
        } else {
            if (len >= min_len && len <= max_len && quality > min_quality) {
                std::osyncstream{out} << read.get_record();
            }
        }
        ++counter;
    }
}

void Work::filter_one_thread(
    const Read& read,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    std::ostream& out)
{
    unsigned len{read.get_length()};
    float quality{read.calculate_read_quality()};
    if (gc) {
        if (float gc_content{read.get_gc_content()};
            len >= min_len &&
            len <= max_len &&
            quality > min_quality &&
            gc_content > min_gc &&
            gc_content < max_gc) {
            out << read.get_record();
        }
    } else {
        if (len >= min_len && len <= max_len && quality > min_quality) {
            out << read.get_record();
        }
    }
}

void Work::filter_multi_fqs_in_one_thread(
    const std::vector<std::filesystem::path>& paths,
    size_t start,
    size_t end,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    std::ostream& out)
{
    for (size_t i{start}; i < end; ++i) {
        int l;
        gzFile file = gzopen(paths[i].c_str(), "rb");
        kseq_t* seq = kseq_init(file);
        while (true) {
            l = kseq_read(seq);
            Read read{FastqReader::fastq_record_ok(l, seq, paths[i].c_str())};
            if (*read.get_id() == finished_read_name) {
                break;
            }
            unsigned len{read.get_length()};
            float quality{read.calculate_read_quality()};
            if (gc) {
                if (float gc_content{read.get_gc_content()};
                    len >= min_len &&
                    len <= max_len &&
                    quality > min_quality &&
                    gc_content > min_gc &&
                    gc_content < max_gc) {
                    std::osyncstream{out} << read.get_record();
                }
            } else {
                if (len >= min_len && len <= max_len && quality > min_quality) {
                    std::osyncstream{out} << read.get_record();
                }
            }
        }
        kseq_destroy(seq);
        gzclose(file);
    }
}

void Work::trim(
    int start,
    int end,
    std::shared_ptr<std::vector<Read>> reads_ptr,
    std::atomic<size_t>& counter,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& align_config,
    std::ostream& log_fstream,
    std::ostream& out,
    std::barrier<>& bar
)
{
    for (int idx{start}; idx < end; idx++) {
        (*reads_ptr)[idx].trim(seq_info, td, align_config, log_fstream);
        std::osyncstream{out} << (*reads_ptr)[idx];
        ++counter;
    }
    bar.arrive_and_wait();
}

void Work::trim_one_thread(
    Read& read,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ostream& log_fstream,
    std::ostream& out)
{
    read.trim(seq_info, td, alignment_config, log_fstream);
    out << read;
}

void Work::trim_multi_fqs_in_one_thread(
    const std::vector<std::filesystem::path>& paths,
    size_t start,
    size_t end,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ostream& log_fstream,
    std::ostream& out)
{
    for (size_t i{start}; i < end; i++) {
        int l;
        gzFile file = gzopen(paths[i].c_str(), "rb");
        kseq_t* seq = kseq_init(file);
        while (true) {
            l = kseq_read(seq);
            Read read{FastqReader::fastq_record_ok(l, seq, paths[i].c_str())};
            if (*read.get_id() == finished_read_name) {
                break;
            }
            read.trim(seq_info, td, alignment_config, log_fstream);
            std::osyncstream{out} << read;
        }
        kseq_destroy(seq);
        gzclose(file);
    }
}

void Work::main(
    int start,
    int end,
    std::shared_ptr<std::vector<Read>> reads_ptr,
    std::vector<read_stats_result>& all_reads_stats_result,
    std::vector<read_stats_result>& passed_reads_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ofstream& out_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    std::ofstream& trim_log_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    std::mutex& all_mtx,
    std::mutex& passed_mtx,
    std::barrier<>& bar)
{
    for (int idx{start}; idx < end; ++idx) {
        Read& read{(*reads_ptr)[idx]};
        main_core(
            read,
            all_reads_stats_result,
            passed_reads_stats_result,
            gc,
            min_len,
            max_len,
            min_quality,
            min_gc,
            max_gc,
            do_trim,
            seq_info,
            td,
            alignment_config,
            out_ofstream,
            all_stats_ofstream,
            passed_stats_ofstream,
            failed_stats_ofstream,
            trim_log_ofstream,
            retain_failed,
            failed_ofstream,
            true,
            all_mtx,
            passed_mtx
        );
    }
    bar.arrive_and_wait();
}


void Work::main_one_thread(
    Read& read,
    std::vector<read_stats_result>& all_reads_stats_result,
    std::vector<read_stats_result>& passed_reads_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ofstream& out_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    std::ofstream& trim_log_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    std::mutex& all_mtx,
    std::mutex& passed_mtx
)
{
    main_core(read,
              all_reads_stats_result,
              passed_reads_stats_result,
              gc,
              min_len,
              max_len,
              min_quality,
              min_gc,
              max_gc,
              do_trim,
              seq_info,
              td,
              alignment_config,
              out_ofstream,
              all_stats_ofstream,
              passed_stats_ofstream,
              failed_stats_ofstream,
              trim_log_ofstream,
              retain_failed,
              failed_ofstream,
              false,
              all_mtx,
              passed_mtx
    );
}

void Work::main_multi_fqs_in_one_thread(
    const std::vector<std::filesystem::path>& paths,
    size_t start,
    size_t end,
    std::vector<read_stats_result>& all_stats_result,
    std::vector<read_stats_result>& passed_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ofstream& out_ofstream,
    std::ofstream& trim_log_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    std::mutex& all_mtx,
    std::mutex& passed_mtx
)
{
    for (size_t i{start}; i < end; ++i) {
        int l;
        gzFile file = gzopen(paths[i].c_str(), "rb");
        kseq_t* seq = kseq_init(file);
        while (true) {
            l = kseq_read(seq);
            Read read{FastqReader::fastq_record_ok(l, seq, paths[i].c_str())};
            if (*read.get_id() == finished_read_name) {
                break;
            }
            main_core(read,
                      all_stats_result,
                      passed_stats_result,
                      gc,
                      min_len,
                      max_len,
                      min_quality,
                      min_gc,
                      max_gc,
                      do_trim,
                      seq_info,
                      td,
                      alignment_config,
                      out_ofstream,
                      all_stats_ofstream,
                      passed_stats_ofstream,
                      failed_stats_ofstream,
                      trim_log_ofstream,
                      retain_failed,
                      failed_ofstream,
                      true,
                      all_mtx,
                      passed_mtx
            );
        }
    }
}

void Work::main_core(
    Read& read,
    std::vector<read_stats_result>& all_reads_stats_result,
    std::vector<read_stats_result>& passed_reads_stats_result,
    bool gc,
    unsigned min_len,
    unsigned max_len,
    float min_quality,
    float min_gc,
    float max_gc,
    bool do_trim,
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& alignment_config,
    std::ofstream& out_ofstream,
    std::ofstream& all_stats_ofstream,
    std::ofstream& passed_stats_ofstream,
    std::ofstream& failed_stats_ofstream,
    std::ofstream& trim_log_ofstream,
    bool retain_failed,
    std::ofstream& failed_ofstream,
    bool sync,
    std::mutex& all_mtx,
    std::mutex& passed_mtx
)
{
    auto name{read.get_id()};
    unsigned before_len{read.get_length()};
    float before_quality{read.calculate_read_quality()};
    float before_gc_content{gc ? read.get_gc_content() : 0.001f};
    // using trimmed read to compare the threshold provided by user. if no need trim, using the raw read
    unsigned after_len{before_len};
    float after_quality{before_quality};
    float after_gc_content{before_gc_content};
    if (do_trim) {
        // std::osyncstream{cout} << std::this_thread::get_id() << ": " << &alignment_config << endl;
        std::string read_name{*read.get_id()};
        if (read_name == "0116d6a2-48f5-4cb5-b6a3-3c2058658c68" ||
            read_name == "00880d01-b743-4fd3-8b1d-e7f11153e1ef"
            ||
            read_name == "0354fa18-c3e4-431c-b632-484e4479ae34" ||
            read_name == "049e55d2-2b4d-47ba-97d3-6354fce7a723" ||
            read_name == "06ce0f44-8850-41b5-95c1-ad8eea0ad2aa" ||
            read_name == "0841876a-1c57-41ec-a6fe-8f0ff7c8e164" ||
            read_name == "084ed772-573d-4f1f-abe7-6387829b292a" ||
            read_name == "08f57dda-29cc-48f8-a7b8-94a58a4f0c20" ||
            read_name == "099852b4-6c57-4d5c-9eb5-66ecb6566650" ||
            read_name == "1102c8b6-e058-4792-b34d-405c872d4291" ||
            read_name == "292a885c-2923-47a7-9080-b443f1e23c0c" ||
            read_name == "86dbc28a-b26c-497c-b730-e568d6a7d0bd" ||
            read_name == "c9196921-429a-4d1a-9956-f7044b91d344" ||
            read_name == "fefc2f2f-2a2e-45a7-8ed7-157ad8496b22"
        ) {
            read.trim(seq_info, td, alignment_config, trim_log_ofstream);
            after_len = read.get_length();
            after_quality = read.calculate_read_quality();
            after_gc_content = gc ? read.get_gc_content() : 0.001f;
        } else {
            read.trim(seq_info, td, alignment_config, trim_log_ofstream);
            after_len = read.get_length();
            after_quality = read.calculate_read_quality();
            after_gc_content = gc ? read.get_gc_content() : 0.001f;
        }
        // read.trim(seq_info, td, alignment_config, trim_log_ofstream);
        // after_len = read.get_length();
        // after_quality = read.calculate_read_quality();
        // after_gc_content = gc ? read.get_gc_content() : 0.001f;
    }
    // {
    //     std::unique_lock lock{m_mtx};
    //     cout << *read.get_id() << "\t" << before_len << "\t" << after_len << '\n';
    // }
    bool passed{
        after_len >= min_len &&
        after_len <= max_len &&
        after_quality > min_quality &&
        after_gc_content > min_gc &&
        after_gc_content < max_gc
    };
    if (!gc) {
        after_gc_content = 0.0f;
        before_gc_content = 0.0f;
    }
    if (sync) {
        if (passed) {
            {
                // output passed reads stats info
                std::unique_lock lock{passed_mtx};
                out_ofstream << read; // output passed reads
                passed_stats_ofstream << fmt::format("{}\t{}\t{}\t{}\n",
                                                     *name,
                                                     after_len,
                                                     after_quality,
                                                     after_gc_content); // output passed reads stats result
                passed_reads_stats_result.emplace_back(
                    name,
                    after_len,
                    after_quality,
                    after_gc_content); // store passed reads stats result
            }
        } else {
            std::osyncstream{failed_stats_ofstream} << fmt::format("{}\t{}\t{}\t{}\n",
                                                                   *name, after_len, after_quality, after_gc_content);
            if (retain_failed) { std::osyncstream{failed_ofstream} << read; }
        }
        {
            std::unique_lock lock{all_mtx};
            all_stats_ofstream << fmt::format(
                "{}\t{}\t{}\t{}\n",
                *name,
                before_len,
                before_quality,
                before_gc_content); // output all reads stats result
            all_reads_stats_result.emplace_back(
                name,
                before_len,
                before_quality,
                before_gc_content); // store all reads stats result
        }
    } else {
        if (passed) {
            out_ofstream << read;
            passed_stats_ofstream << fmt::format(
                "{}\t{}\t{}\t{}\n",
                *name, after_len, after_quality, after_gc_content
            );
            passed_reads_stats_result.emplace_back(name, after_len, after_quality, after_gc_content);
        } else {
            failed_stats_ofstream << fmt::format("{}\t{}\t{}\t{}\n",
                                                 *name, after_len, after_quality, after_gc_content);
            if (retain_failed) { failed_ofstream << read; }
        }
        all_stats_ofstream << fmt::format(
            "{}\t{}\t{}\t{}\n",
            *name,
            before_len,
            before_quality,
            before_gc_content); // output all reads stats result
        all_reads_stats_result.emplace_back(
            name,
            before_len,
            before_quality,
            before_gc_content); // store all reads stats result
    }
}
