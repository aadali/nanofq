#include "Work.h"

#include <thread>
#include <fmt/core.h>

using std::cout;
using std::endl;

Work::Work(FastqReader &fq, const unsigned thread, const bool gc, std::string_view outfile_path) :
        m_fq(fq), m_thread{thread}, m_gc(gc), m_outfile_path(outfile_path),
        m_bar(static_cast<int>(m_thread)) {
    m_sub_stats_result.reserve(m_thread);
    for (int i{0}; i < m_thread; i++) {
        m_sub_stats_result.emplace_back();
    }
    m_outfile_stream = std::ofstream(m_outfile_path.data(), std::ios::out);
    if (!m_outfile_stream) {
        throw std::runtime_error(fmt::format("Cannot open file: {}", m_outfile_path));
    }
}

std::vector<std::pair<unsigned , unsigned >> Work::get_bins(const unsigned length) const {
    std::vector<std::pair<unsigned , unsigned >> idx_ranges;
    const unsigned step{length / m_thread};
    unsigned start{0}, stop{0};
    for (unsigned i{0}; i < length; i += step) {
        start = i;
        stop = step + start;
        if (stop + step > length) {
            stop = length;
            idx_ranges.emplace_back(start, stop);
            break;
        }
        idx_ranges.emplace_back(start, stop);
    }
    return idx_ranges;
}

void Work::run_stats() {
    while (true) {
        if (std::optional<shared_vec_reads> reads = m_fq.get_reads(); reads.has_value()) {
            std::vector<std::jthread> threads;
            threads.reserve(m_thread);
            auto bins{get_bins(reads.value()->size())};
            for (int i{0}; i < m_thread; i++) {
                threads.emplace_back(&Work::stats,
                                     this,
                                     bins[i].first,
                                     bins[i].second,
                                     reads.value(),
                                     std::ref(m_sub_stats_result[i]));
            }
        }
        if (m_fq.read_finish() && m_fq.is_empty()) break;
    }
    for (std::vector<read_stats_result> &item: m_sub_stats_result) {
        for (read_stats_result &x: item) {
            auto line = fmt::format("{}\t{}\t{}\t{}\n", std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x));
            m_outfile_stream << line;
        }
    }
//    cout << m_stats_result.size() << endl;
}


void Work::run_filter(const unsigned min_len,
                      const unsigned max_len,
                      const float min_quality,
                      const float min_gc,
                      const float max_gc) {
    while (true) {
        if (std::optional<shared_vec_reads> reads = m_fq.get_reads(); reads.has_value()) {
            std::vector<std::jthread> threads;
            threads.reserve(m_thread);
            auto bins{get_bins(reads.value()->size())};;
            for (int i{0}; i < m_thread; i++) {
                threads.emplace_back(&Work::filter,
                                     this,
                                     bins[i].first,
                                     bins[i].second,
                                     min_len,
                                     max_len,
                                     min_quality,
                                     min_gc,
                                     max_gc,
                                     reads.value());
            }
        }
        if (m_fq.read_finish()) break;
    }
}

void Work::run_find(const std::string& input_reads, bool use_index, unsigned key_length) {
    m_fq.find_reads(input_reads, m_outfile_stream, use_index, key_length);
}

void Work::run_index(unsigned key_length) {
    m_fq.index(key_length);
}

void Work::stats(const unsigned start,
                 const unsigned end,
                 const shared_vec_reads &reads,
                 std::vector<read_stats_result> &sub_stats_result) {
    for (unsigned idx{start}; idx < end; idx++) {
        unsigned len{(*reads)[idx]->get_length()};
        float quality{(*reads)[idx]->calculate_read_quality()};
        float gc_content{m_gc ? (*reads)[idx]->get_gc_content() : 0.0f};
        sub_stats_result.emplace_back((*reads)[idx]->get_id(), len, quality, gc_content);
    }
    m_bar.arrive_and_wait();
}


void Work::filter(const unsigned start,
                  const unsigned end,
                  const unsigned min_len,
                  const unsigned max_len,
                  const float min_quality,
                  const float min_gc,
                  const float max_gc,
                  const shared_vec_reads &reads) {
    for (unsigned idx{start}; idx < end; idx++) {
        unsigned len{(*reads)[idx]->get_length()};
        float quality{(*reads)[idx]->calculate_read_quality()};
        if (m_gc) {
            if (float gc_content{(*reads)[idx]->get_gc_content()};
                    len >= min_len && len <= max_len && quality > min_quality &&
                    gc_content > min_gc && gc_content < max_gc) {
                std::osyncstream{m_outfile_stream} << (*reads)[idx]->get_record();
            }
        } else {
            if (len >= min_len && len <= max_len && quality > min_quality) {
                std::osyncstream{m_outfile_stream} << (*reads)[idx]->get_record();
            }
        }
    }
    m_bar.arrive_and_wait();
}

Work::~Work() {
    if (m_outfile_stream.is_open()) m_outfile_stream.close();
}

