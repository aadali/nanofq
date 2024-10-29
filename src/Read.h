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

static const std::array<double, 256> s_char_to_score_table = [] {
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

    Read(std::string &id, std::string &desc, std::string &sequence, std::string &quality) :
            m_id(std::move(id)), m_desc(std::move(desc)), m_sequence(std::move(sequence)),
            m_quality(std::move(quality)) {}

    Read(char *id, char *desc, char *sequence, char *quality) :
            m_id(id), m_desc(desc ? desc : ""), m_sequence(sequence), m_quality(quality) {}

    Read(const Read &read) = delete;

    Read &operator=(const Read &read) = delete;

    Read(Read &&read) = default;

    Read &operator=(Read &&read) = delete;

    ~Read() = default;

public:
    [[nodiscard]] const std::string &get_id() const { return this->m_id; }

    [[nodiscard]] unsigned get_length() const { return m_sequence.size(); }

    [[nodiscard]] float get_gc_content() const {
        auto gc_number{
                std::ranges::count_if(m_sequence,
                                      [](const char &c) { return c == 'G' || c == 'C' || c == 'g' || c == 'c'; })
        };
        return static_cast<float >(gc_number) / static_cast<float >(m_sequence.size());
    }

    void rev_com() {
        std::ranges::transform(m_sequence,
                               std::begin(m_sequence),
                               [](const char &c) {
                                   switch (c) {
                                       case 'A':
                                           return 'T';
                                       case 'T':
                                           return 'A';
                                       case 'G':
                                           return 'C';
                                       case 'C':
                                           return 'G';
                                       case 'a':
                                           return 't';
                                       case 't':
                                           return 'a';
                                       case 'g':
                                           return 'c';
                                       case 'c':
                                           return 'a';
                                       default:
                                           return 'N';
                                   }
                               });
        std::ranges::reverse(m_sequence);
        std::ranges::reverse(m_quality);
    }

    [[nodiscard]] float calculate_read_quality() const {
        auto error_p{m_quality | std::views::transform([](const char &c) { return s_char_to_score_table[c]; })};
        double total_error_p{
                std::accumulate(std::cbegin(error_p), std::cend(error_p), 0.0) / static_cast<float >(m_quality.size())
        };
        return std::log10(total_error_p) * -10.0;
    }

    [[nodiscard]] bool is_passed(const unsigned min_length, const unsigned max_length, const float quality) const {
        return m_sequence.size() >= min_length
               && m_sequence.size() <= max_length
               && calculate_read_quality() > quality;
    }

    [[nodiscard]] bool is_passed(const unsigned min_length, const unsigned max_length, float quality, float min_gc,
                                 float max_gc) const {
        float gc{get_gc_content()};
        return m_sequence.size() >= min_length
               && m_sequence.size() <= max_length
               && calculate_read_quality() > quality
               && gc > min_gc
               && gc < max_gc;
    }

    [[nodiscard]] std::string get_record() const {
        if (m_desc.empty()) {
            return fmt::format("@{}\n{}\n+\n{}\n", m_id, m_sequence, m_quality);
        }
        return fmt::format("@{} {}\n{}\n+\n{}\n", m_id, m_desc, m_sequence, m_quality);
    }
};

inline std::ostream &operator<<(std::ostream &c, const Read &read) {
    c << read.m_id + read.m_desc << '\n';
    c << read.m_sequence + '\n' << read.m_quality;
    return c;
}
