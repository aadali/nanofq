#include <syncstream>
#include "Read.h"

Read::Read(std::string& id, std::string& desc, std::string& sequence, std::string& quality) :
    m_id(std::move(id)), m_desc(std::move(desc)), m_sequence(std::move(sequence)),
    m_quality(std::move(quality)) {}

Read::Read(char* id, char* desc, char* sequence, char* quality) :
    m_id(id), m_desc(desc ? desc : ""), m_sequence(sequence), m_quality(quality) {}

float Read::get_gc_content() const {
    auto gc_number{
        std::ranges::count_if(m_sequence,
                              [](const char& c){ return c == 'G' || c == 'C' || c == 'g' || c == 'c'; })
    };
    return static_cast<float>(gc_number) / static_cast<float>(m_sequence.size());
}

void Read::rev_com() {
    std::ranges::transform(m_sequence,
                           std::begin(m_sequence),
                           [](const char& c){
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

float Read::calculate_read_quality() const {
    auto error_p{m_quality | std::views::transform([](const char& c){ return s_char_to_score_table[c]; })};
    double total_error_p{
        std::accumulate(std::cbegin(error_p), std::cend(error_p), 0.0) / static_cast<float>(m_quality.size())
    };
    return std::log10(total_error_p) * -10.0;
}

bool Read::is_passed(const unsigned min_length, const unsigned max_length, const float quality) const {
    return m_sequence.size() >= min_length
        && m_sequence.size() <= max_length
        && calculate_read_quality() > quality;
}

bool Read::is_passed(const unsigned min_length, const unsigned max_length, float quality, float min_gc,
                     float max_gc) const {
    float gc{get_gc_content()};
    return m_sequence.size() >= min_length
        && m_sequence.size() <= max_length
        && calculate_read_quality() > quality
        && gc > min_gc
        && gc < max_gc;
}

std::string Read::get_record() const {
    if (m_desc.empty()) {
        return fmt::format("@{}\n{}\n+\n{}\n", m_id, m_sequence, m_quality);
    }
    return fmt::format("@{} {}\n{}\n+\n{}\n", m_id, m_desc, m_sequence, m_quality);
}

Read& Read::trim(const SequenceInfo& seq_info, const trim_direction& td,
                 AlignmentConfig& align_config,
                 std::ostream& log) {
    std::string_view sequence_view{m_sequence};
    unsigned trim_start_idx{0};
    unsigned trim_end_idx{static_cast<unsigned>(m_sequence.size())};
    bool already_trimmed{false};

    if (td.trim_top5end) {
        std::string_view target_5end_view{
            sequence_view.size() > get<0>(seq_info.m_top5end)
                ? sequence_view.substr(0, get<0>(seq_info.m_top5end))
                : sequence_view
        };
        AlignmentResult align_result_top5end;
        myUtility::smith_waterman(target_5end_view, seq_info.m_top5end_query, align_config, align_result_top5end);
        if (align_result_top5end.get_percent(seq_info.m_top5end_query) > get<1>(seq_info.m_top5end) &&
            align_result_top5end.get_identity() > get<2>(seq_info.m_top5end)) {
            already_trimmed = true;
            auto [query_5start, target_5start] = align_result_top5end.get_start_idx();
            auto [query_5end, target_5end] = align_result_top5end.get_stop_idx();
            std::osyncstream{log} << fmt::format("{} Left\t{}\t{}\t{}\t{}\n{}\n",
                                                 m_id,
                                                 target_5start,
                                                 target_5end,
                                                 query_5start,
                                                 query_5end,
                                                 align_result_top5end.to_string()
            );
            if (target_5end >= sequence_view.size()) { // target_5end will never is bigger than sequence_view.size();
                m_sequence = "";
                m_quality = "";
                return *this;
            }
            trim_start_idx = target_5end;
            sequence_view.remove_prefix(trim_start_idx);
        }

        if (td.trim_top3end) {
            std::string_view target_3end_view{
                sequence_view.size() > get<0>(seq_info.m_top3end)
                    ? sequence_view.substr(sequence_view.size() - get<0>(seq_info.m_top3end))
                    : sequence_view
            };
            AlignmentResult align_result_top3end;
            myUtility::smith_waterman(target_3end_view, seq_info.m_top3end_query, align_config, align_result_top3end);
            if (align_result_top3end.get_percent(seq_info.m_top3end_query) > get<1>(seq_info.m_top3end) &&
                align_result_top3end.get_identity() > get<2>(seq_info.m_top3end)) {
                already_trimmed = true;
                auto [query_3start, target_3start] = align_result_top3end.get_start_idx();
                auto [query_3end, target_3end] = align_result_top3end.get_stop_idx();
                std::osyncstream{log} << fmt::format("{} Right\t-{}\t-{}\t{}\t{}\n{}\n",
                                                     m_id,
                                                     target_3end_view.size() - target_3start,
                                                     target_3end_view.size() - target_3end,
                                                     query_3start,
                                                     query_3end,
                                                     align_result_top3end.to_string()
                );
                if (target_3end >= sequence_view.size()) {
                    m_sequence = "";
                    m_quality = "";
                    return *this;
                }
                trim_end_idx = m_sequence.size() - target_3start; // the bases number should be trimmed from read right
            }
        }
    }


    if (td.trim_bot5end && !already_trimmed) {
        std::string_view target_5end_view{
            sequence_view.size() > get<0>(seq_info.m_bot5end)
                ? sequence_view.substr(0, get<0>(seq_info.m_bot5end))
                : sequence_view
        };
        AlignmentResult align_result_top5end;
        myUtility::smith_waterman(target_5end_view, seq_info.m_bot5end_query, align_config, align_result_top5end);
        if (align_result_top5end.get_percent(seq_info.m_bot5end_query) > get<1>(seq_info.m_bot5end) &&
            align_result_top5end.get_identity() > get<2>(seq_info.m_bot5end)) {
            auto [query_5start, target_5start] = align_result_top5end.get_start_idx();
            auto [query_5end, target_5end] = align_result_top5end.get_stop_idx();
            std::osyncstream{log} << fmt::format("{} Left\t{}\t{}\t{}\t{}\n{}\n",
                                                 m_id,
                                                 target_5start,
                                                 target_5end,
                                                 query_5start,
                                                 query_5end,
                                                 align_result_top5end.to_string()
            );
            if (target_5end >= sequence_view.size()) {
                m_sequence = "";
                m_quality = "";
                return *this;
            }
            trim_start_idx = target_5end;
            sequence_view.remove_prefix(trim_start_idx);
        }

        if (td.trim_bot3end) {
            std::string_view target_3end_view{
                sequence_view.size() > get<0>(seq_info.m_bot3end)
                    ? sequence_view.substr(sequence_view.size() - get<0>(seq_info.m_bot3end))
                    : sequence_view
            };
            AlignmentResult align_result_bot3end;
            myUtility::smith_waterman(target_3end_view, seq_info.m_bot3end_query, align_config, align_result_bot3end);
            if (align_result_bot3end.get_percent(seq_info.m_bot3end_query) > get<1>(seq_info.m_bot3end) &&
                align_result_bot3end.get_identity() > get<2>(seq_info.m_bot3end)) {
                auto [query_3start, target_3start] = align_result_bot3end.get_start_idx();
                auto [query_3end, target_3end] = align_result_bot3end.get_stop_idx();
                std::osyncstream{log} << fmt::format("{} Right\t-{}\t-{}\t{}\t{}\n{}\n",
                                                     m_id,
                                                     target_3end_view.size() - target_3start,
                                                     target_3end_view.size() - target_3end,
                                                     query_3start,
                                                     query_3end,
                                                     align_result_bot3end.to_string());
                if (target_3end >= sequence_view.size()) {
                    m_sequence = "";
                    m_quality = "";
                    return *this;
                }
                trim_end_idx = m_sequence.size() - target_3end;
            }
        }
    }

    m_sequence.erase(trim_start_idx, trim_end_idx - trim_start_idx);
    m_quality.erase(trim_start_idx, trim_end_idx - trim_start_idx);
    return *this;
}
