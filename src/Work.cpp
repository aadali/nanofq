#include "Work.h"

#include <thread>
#include <fmt/core.h>

using std::cout;
using std::endl;

Work::Work(FastqReader& fq, unsigned thread, const bool gc, std::ostream& out) :
    m_fq(fq),
    m_thread{thread},
    m_gc(gc),
    m_out(out),
    m_bar(static_cast<int>(m_thread)) {
    m_stats_result.reserve(m_thread);
    for (int i{0}; i < m_thread; i++) {
        m_stats_result.emplace_back();
    }
}

Work::Work(FastqReader& fq): m_fq{fq}, m_bar(1) {}


std::vector<std::pair<unsigned, unsigned>> Work::get_bins(unsigned length) const {
    std::vector<std::pair<unsigned, unsigned>> idx_ranges;
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
    if (m_thread == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (read.get_id() == finished_read_name) return;
            stats_one_thread(read);
        }
    }
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
                                     std::ref(m_stats_result[i]));
            }
        }
        if (m_fq.read_finish() && m_fq.is_empty()) break;
    }
    for (std::vector<read_stats_result>& item : m_stats_result) {
        for (read_stats_result& x : item) {
            auto line = fmt::format("{}\t{}\t{}\t{}\n",
                                    std::get<0>(x),
                                    std::get<1>(x),
                                    fmt::format("{:.{}f}", std::get<2>(x), 2),
                                    fmt::format("{:.{}f}", std::get<3>(x), 2));
            m_out << line;
        }
    }
}


void Work::run_filter(unsigned min_len,
                      unsigned max_len,
                      float min_quality,
                      float min_gc,
                      float max_gc) {
    if (m_thread == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (read.get_id() == finished_read_name) return;
            filter_one_thread(read, min_len, max_len, min_quality, min_gc, max_gc);
        }
    }
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
        if (m_fq.read_finish() && m_fq.is_empty()) break;
    }
}

void Work::run_trim(const SequenceInfo& seq_info,
                    const trim_direction& td,
                    std::vector<AlignmentConfig>& align_configs,
                    std::ostream& log_fstream) {
    if (m_thread == 1) {
        while (true) {
            Read read{m_fq.read_one_fastq()};
            if (read.get_id() == finished_read_name) return;
            trim_one_thread(read, seq_info, td, align_configs[0], log_fstream);
        }
    }
    while (true) {
        if (std::optional<shared_vec_reads> reads = m_fq.get_reads(); reads.has_value()) {
            std::vector<std::jthread> threads;
            threads.reserve(m_thread);
            auto bins{get_bins(reads.value()->size())};
            for (int i{0}; i < m_thread; i++) {
                threads.emplace_back(&Work::trim,
                                     this,
                                     bins[i].first,
                                     bins[i].second,
                                     reads.value(),
                                     std::ref(seq_info),
                                     std::ref(td),
                                     std::ref(align_configs[i]), std::ref(log_fstream));
            }
            // trim(0, reads.value()->size(), reads.value(), seq_info, td, align_config, log_fstream);
        }
        if (m_fq.read_finish() && m_fq.is_empty()) break;
    }
}

void Work::run_find(const std::string& input_reads, bool use_index, unsigned key_length) const {
    m_fq.find_reads(input_reads, m_out, use_index, key_length);
}

void Work::run_index(unsigned key_length) const {
    m_fq.index(key_length);
}


void Work::stats(unsigned start,
                 unsigned end,
                 const shared_vec_reads& reads,
                 std::vector<read_stats_result>& sub_stats_result) {
    for (unsigned idx{start}; idx < end; idx++) {
        unsigned len{(*reads)[idx]->get_length()};
        float quality{(*reads)[idx]->calculate_read_quality()};
        float gc_content{m_gc ? (*reads)[idx]->get_gc_content() : 0.0f};
        sub_stats_result.emplace_back((*reads)[idx]->get_id(), len, quality, gc_content);
    }
    m_bar.arrive_and_wait();
}

void Work::stats_one_thread(const Read& read) {
    unsigned len{read.get_length()};
    float quality{read.calculate_read_quality()};
    float gc_content{m_gc ? read.get_gc_content() : 0.0f};
    m_stats_result[0].emplace_back(read.get_id(), len, quality, gc_content);
    auto line = fmt::format("{}\t{}\t{}\t{}\n",
                            read.get_id(),
                            len,
                            fmt::format("{:.{}f}", quality, 2),
                            fmt::format("{:.{}f}", gc_content, 2));
    m_out << line;
}


void Work::filter(unsigned start,
                  unsigned end,
                  unsigned min_len,
                  unsigned max_len,
                  float min_quality,
                  float min_gc,
                  float max_gc,
                  const shared_vec_reads& reads) {
    for (unsigned idx{start}; idx < end; idx++) {
        unsigned len{(*reads)[idx]->get_length()};
        float quality{(*reads)[idx]->calculate_read_quality()};
        if (m_gc) {
            if (float gc_content{(*reads)[idx]->get_gc_content()};
                len >= min_len && len <= max_len && quality > min_quality &&
                gc_content > min_gc && gc_content < max_gc) {
                std::osyncstream{m_out} << (*reads)[idx]->get_record();
            }
        } else {
            if (len >= min_len && len <= max_len && quality > min_quality) {
                std::osyncstream{m_out} << (*reads)[idx]->get_record();
            }
        }
    }
    m_bar.arrive_and_wait();
}

void Work::filter_one_thread(const Read& read, unsigned min_len, unsigned max_len, float min_quality, float min_gc,
                             float max_gc) const {
    unsigned len{read.get_length()};
    float quality{read.calculate_read_quality()};
    if (m_gc) {
        if (float gc_content{read.get_gc_content()};
            len >= min_len &&
            len <= max_len &&
            quality > min_quality &&
            gc_content > min_gc &&
            gc_content < max_gc) {
            m_out << read.get_record();
        }
    } else {
        if (len >= min_len && len <= max_len && quality > min_quality) {
            m_out << read.get_record();
        }
    }
}


void Work::trim(unsigned start, unsigned end, const shared_vec_reads& reads, const SequenceInfo& seq_info,
                const trim_direction& td, AlignmentConfig& align_config, std::ostream& log_fstream) {
    for (unsigned idx{start}; idx < end; idx++) {
        (*reads)[idx]->trim(seq_info, td, align_config, log_fstream);
        std::osyncstream{m_out} << *(*reads)[idx];
    }
    m_bar.arrive_and_wait();
}

void Work::trim_one_thread(Read& read, const SequenceInfo& seq_info, const trim_direction& td,
                           AlignmentConfig& alignment_config, std::ostream& log_fstream) const {
    read.trim(seq_info, td, alignment_config, log_fstream);
    m_out << read;
}
