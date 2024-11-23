#pragma once

#include <algorithm>
#include <ranges>
#include <string>
#include <iostream>
#include <array>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <fmt/core.h>
#include "SequenceInfo.h"
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
#include "myUtility.h"

static const std::array<double, 256> s_char_to_score_table = []{
    std::array<double, 256> a{};
    for (int q{33}; q < 127; q++) {
        a[q] = std::pow(10.0, static_cast<double>(33 - q) / 10.0);
    }
    return a;
}();

class Read {
public:
    std::string m_id;
    std::string m_desc;
    std::string m_sequence;
    std::string m_quality;

public:
    Read() = delete;

    Read(std::string& id, std::string& desc, std::string& sequence, std::string& quality);

    Read(char* id, char* desc, char* sequence, char* quality);

    Read(const Read& read) = delete;

    Read& operator=(const Read& read) = delete;

    Read(Read&& read) noexcept ;

    Read& operator=(Read&& read) noexcept;

    ~Read() = default;

public:
    [[nodiscard]] const std::string& get_id() const { return this->m_id; }
    [[nodiscard]] const std::string& get_sequence() const { return this->m_sequence; }
    [[nodiscard]] const std::string& get_quality() const { return this->m_quality; }


    [[nodiscard]] unsigned get_length() const { return m_sequence.size(); }

    [[nodiscard]] float get_gc_content() const;

    void rev_com();

    [[nodiscard]] float calculate_read_quality() const;
    [[nodiscard]] bool is_passed(const unsigned min_length, const unsigned max_length, const float quality) const;


    [[nodiscard]] bool is_passed(const unsigned min_length, const unsigned max_length, float quality, float min_gc,
                                 float max_gc) const;

    [[nodiscard]] std::string get_record() const;
    void trim(const SequenceInfo& seq_info, const trim_direction& td, AlignmentConfig& align_config,
               std::ostream& log);

private:
    size_t trim_positive_strand_left(std::string_view top5end_query, const trim_end& top5end, AlignmentConfig& align_config,
                           AlignmentResult& align_5end_result) const;
    size_t trim_positive_strand_right(std::string_view& left_trimmed_seq_view, std::string_view top3end_query, const trim_end& top3end, AlignmentConfig& align_config,
                            AlignmentResult& align_3end_result) const;
    size_t trim_negative_strand_left(std::string_view bot5end_query, const trim_end& bot5end, AlignmentConfig& align_config,
                           AlignmentResult& align_5end_result) const;
    size_t trim_negative_strand_right(std::string_view& left_trimmed_seq_view, std::string_view bot3end_query, const trim_end& bot3end, AlignmentConfig& align_config,
                            AlignmentResult& align_3end_result) const;
};

inline std::ostream& operator<<(std::ostream& c, const Read& read) {
    if (read.m_desc.empty()) {
        c << fmt::format("@{}\n{}\n+\n{}\n", read.m_id, read.m_sequence, read.m_quality);
    } else {
        c << fmt::format("@{} {}\n{}\n+\n{}\n", read.m_id, read.m_desc, read.m_sequence, read.m_quality);
    }
    return c;
}
