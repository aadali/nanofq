#include "Work.h"

#include <thread>
#include <fmt/core.h>

#include "fmt/xchar.h"

using std::cout;
using std::endl;

Work::Work(FastqReader& fq, ThreadPool& threads_pool) : m_fq(fq), m_threads_pool(threads_pool) {}

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

void Work::run_stats(std::vector<read_stats_result>& stats_result, std::osyncstream& out, bool gc)
{
    if (m_threads_pool.threads_number() == 1) {
        if (!stats_result.empty()) {
            cerr << "When used one thread, the parameter stats_result must be empty" << endl;
            exit(1);
        }
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (read.get_id() == finished_read_name)
                return;
            stats_one_thread(read, stats_result, out, gc);
        }
    }
    size_t total_reads{0};
    while (true) {
        auto reads_ptr{m_fq.read_chunk_fastq()};
        total_reads += reads_ptr->size();
        auto bins = get_edges(reads_ptr->size());
        for (auto [start, end] : bins) {
            m_threads_pool.enqueue([this, start, end, gc, &stats_result, &out, reads_ptr](){
                stats(start, end, reads_ptr, stats_result, out, gc);
            });
        }
        if (m_fq.read_finish())
            break;
    }
    while (total_reads != stats_result.size()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
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
    std::osyncstream& out) const
{
    if (m_threads_pool.threads_number() == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (read.get_id() == finished_read_name) return;
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
        if (m_fq.read_finish()) break;
    }
    while (total_reads != counter) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
}

//
// void Work::run_trim(const SequenceInfo& seq_info,
//                     const trim_direction& td,
//                     std::vector<AlignmentConfig>& align_configs,
//                     std::ostream& log_fstream)
// {
//     if (m_thread == 1) {
//         while (true) {
//             Read read{m_fq.read_one_fastq()};
//             if (read.get_id() == finished_read_name) return;
//             trim_one_thread(read, seq_info, td, align_configs[0], log_fstream);
//         }
//     }
//     while (true) {
//         if (std::optional<shared_vec_reads> reads = m_fq.get_reads(); reads.has_value()) {
//             std::vector<std::jthread> threads;
//             threads.reserve(m_thread);
//             auto bins{get_bins(reads.value()->size())};
//             for (int i{0}; i < m_thread; i++) {
//                 threads.emplace_back(&Work::trim,
//                                      this,
//                                      bins[i].first,
//                                      bins[i].second,
//                                      reads.value(),
//                                      std::ref(seq_info),
//                                      std::ref(td),
//                                      std::ref(align_configs[i]), std::ref(log_fstream));
//             }
//             // trim(0, reads.value()->size(), reads.value(), seq_info, td, align_config, log_fstream);
//         }
//         if (m_fq.read_finish() && m_fq.is_empty()) break;
//     }
// }

void Work::save_summary(int n, const std::vector<int>& read_quals, const std::vector<int>& read_lengths,
                        std::vector<read_stats_result>& stats_result,
                        const std::string& summary_file_path)
{
    std::string summary_info{summary_stats_result(n, read_quals, read_lengths, stats_result)};
    std::ofstream summary_file{summary_file_path, std::ios::out};
    if (summary_file) {
        summary_file << summary_info;
        summary_info.clear();
    } else {
        std::cerr << WARNS + "Failed when opening " + summary_file_path + ". Try write the summary into ~/SumMarY.txt" +
            COLOR_END
            << std::endl;
        std::ofstream try_summary_file{summary_file_path, std::ios::out};
        try {
            try_summary_file << summary_info;
            try_summary_file.close();
        }
        catch (const std::exception& e) {
            std::cerr << REDS << "Failed when trying write summary into ~/SumMarY.txt" + COLOR_END << std::endl;
            exit(1);
        }
    }
}

// void Work::run_find(const std::string& input_reads, bool use_index, unsigned key_length) const
// {
//     m_fq.find_reads(input_reads, m_out, use_index, key_length);
// }
//
void Work::run_index(unsigned key_length) const
{
    m_fq.index(key_length);
}

void Work::stats(int start,
                 int end,
                 std::shared_ptr<std::vector<Read>> reads_ptr,
                 std::vector<read_stats_result>& stats_results,
                 std::osyncstream& out,
                 bool gc)
{
    for (int idx{start}; idx < end; idx++) {
        const Read& read = (*reads_ptr)[idx];
        std::string name{read.get_id()};
        unsigned length{read.get_length()};
        float quality{read.calculate_read_quality()};
        float gc_content{gc ? read.get_gc_content() : 0.0f};
        out << fmt::format("{}\t{}\t{}\t{}\n", name, length, quality, gc_content);
        {
            std::unique_lock lock{m_mtx};
            stats_results.emplace_back(name, length, quality, gc_content);
        }
    }
}

std::string Work::summary_stats_result(int n, const std::vector<int>& read_quals,
                                       const std::vector<int>& read_lengths,
                                       std::vector<read_stats_result>& stats_result)
{
    std::stringstream summary_stream;
    // for (int i{1}; i < m_stats_result.size(); i++) {
    //     for (read_stats_result& item : m_stats_result[i]) {
    //         m_stats_result[0].push_back(std::move(item));
    //     }
    // }
    // std::vector<read_stats_result>& stats_result = m_stats_result[0];
    // std::span<read_stats_result> stats_result_span{m_stats_result[0]};
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
        longest_reads << fmt::format("{}\t{}\t{}\t{}\t{}\n",
                                     i + 1,
                                     std::get<0>(stats_result[i]),
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
                    lengths_info << fmt::format("ReadLength > {}\t{}({});{}({})\n",
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
            for (; read_idx < total_reads_number; read_idx++) {
                if (std::get<2>(stats_result[read_idx]) < quality || read_idx == total_reads_number - 1) {
                    summary_stream << fmt::format("ReadQuality > {}\t{}({});{}({})\n",
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
                reads_count += 1;
                bases_count += std::get<1>(stats_result[read_idx]);
            }
        }
    };
    for (const int& quality : read_quals) {
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
                                      std::get<0>(stats_result[i]),
                                      std::get<1>(stats_result[i]),
                                      fmt::format("{:.2f}", std::get<2>(stats_result[i])),
                                      fmt::format("{:.2f}", std::get<3>(stats_result[i])));
    }
    return summary_stream.str();
}

void Work::stats_one_thread(const Read& read, std::vector<read_stats_result>& stats_result, std::ostream& out, bool gc)
{
    unsigned len{read.get_length()};
    float quality{read.calculate_read_quality()};
    float gc_content{gc ? read.get_gc_content() : 0.0f};
    stats_result.emplace_back(read.get_id(), len, quality, gc_content);
    auto line = fmt::format("{}\t{}\t{}\t{}\n",
                            read.get_id(),
                            len,
                            fmt::format("{:.{}f}", quality, 2),
                            fmt::format("{:.{}f}", gc_content, 2));
    out << line;
}

void Work::filter(int start,
                  int end,
                  std::shared_ptr<std::vector<Read>> reads_ptr,
                  std::atomic<size_t>& counter,
                  bool gc,
                  unsigned min_len,
                  unsigned max_len,
                  float min_quality,
                  float min_gc,
                  float max_gc,
                  std::osyncstream& out)
{
    for (int idx{start}; idx < end; idx++) {
        const Read& read = (*reads_ptr)[idx];
        unsigned len{read.get_length()};
        float quality{read.calculate_read_quality()};
        if (gc) {
            if (float gc_content{read.get_gc_content()};
                len >= min_len && len <= max_len && quality > min_quality &&
                gc_content > min_gc && gc_content < max_gc) {
                out << read.get_record();
            }
        } else {
            if (len >= min_len && len <= max_len && quality > min_quality) {
                out<< read.get_record();
            }
        }
        ++counter;
    }
}

void Work::filter_one_thread(const Read& read,
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

/*
void Work::trim(unsigned start, unsigned end, const shared_vec_reads& reads, const SequenceInfo& seq_info,
                const trim_direction& td, AlignmentConfig& align_config, std::ostream& log_fstream)
{
    for (unsigned idx{start}; idx < end; idx++) {
        (*reads)[idx]->trim(seq_info, td, align_config, log_fstream);
        std::osyncstream{m_out} << *(*reads)[idx];
    }
    m_bar.arrive_and_wait();
}

void Work::trim_one_thread(Read& read, const SequenceInfo& seq_info, const trim_direction& td,
                           AlignmentConfig& alignment_config, std::ostream& log_fstream) const
{
    read.trim(seq_info, td, alignment_config, log_fstream);
    m_out << read;
}

*/
