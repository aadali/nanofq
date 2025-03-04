#include <syncstream>
#include "Read.h"

#include <thread>

Read::Read(
    std::string& id,
    std::string& desc,
    std::string& sequence,
    std::string& quality):
    m_id(std::make_shared<std::string>(id)),
    m_desc(std::make_shared<std::string>(desc)),
    m_sequence(std::make_shared<std::string>(sequence)),
    m_quality(std::make_shared<std::string>(quality)) {}

Read::Read(
    char* id,
    char* desc,
    char* sequence,
    char* quality):
    m_id(std::make_shared<std::string>(id)),
    m_desc(desc ? std::make_shared<std::string>(desc): std::make_shared<std::string>("")),
    m_sequence(std::make_shared<std::string>(sequence)),
    m_quality(std::make_shared<std::string>(quality)) {}

Read& Read::operator=(Read&& read) noexcept
{
    if (&read != this) {
        m_id = std::move(read.m_id);
        m_desc = std::move(read.m_desc);
        m_sequence = std::move(read.m_sequence);
        m_quality = std::move(read.m_quality);
    }
    return *this;
}

Read::Read(Read&& read) noexcept
{
    m_id = std::move(read.m_id);
    m_desc = std::move(read.m_desc);
    m_sequence = std::move(read.m_sequence);
    m_quality = std::move(read.m_quality);
}

float Read::get_gc_content() const
{
    auto gc_number{
        std::ranges::count_if(*m_sequence,
                              [](const char& c){ return c == 'G' || c == 'C' || c == 'g' || c == 'c'; })
    };
    return static_cast<float>(gc_number) / static_cast<float>(m_sequence->size());
}

void Read::rev_com()
{
    std::ranges::transform(*m_sequence,
                           std::begin(*m_sequence),
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
    std::ranges::reverse(*m_sequence);
    std::ranges::reverse(*m_quality);
}

float Read::calculate_read_quality() const
{
    auto error_probability{*m_quality | std::views::transform([](const char& c){ return s_char_to_score_table[c]; })};
    double total_error_probability{
        std::accumulate(std::cbegin(error_probability), std::cend(error_probability), 0.0) / static_cast<float>(
            m_quality->size())
    };
    return std::log10(total_error_probability) * -10.0;
}

bool Read::is_passed(
    const unsigned min_length,
    const unsigned max_length,
    const float quality) const
{
    return m_sequence->size() >= min_length
        && m_sequence->size() <= max_length
        && calculate_read_quality() > quality;
}

bool Read::is_passed(
    const unsigned min_length,
    const unsigned max_length,
    float quality,
    float min_gc,
    float max_gc) const
{
    float gc{get_gc_content()};
    return m_sequence->size() >= min_length
        && m_sequence->size() <= max_length
        && calculate_read_quality() > quality
        && gc > min_gc
        && gc < max_gc;
}

std::string Read::get_record() const
{
    if (m_desc->empty()) {
        return fmt::format("@{}\n{}\n+\n{}\n", *m_id, *m_sequence, *m_quality);
    }
    return fmt::format("@{} {}\n{}\n+\n{}\n", *m_id, *m_desc, *m_sequence, *m_quality);
}


size_t Read::trim_positive_strand_left(
    const std::string& top5end_query,
    const trim_end& top5end,
    AlignmentConfig& align_config,
    AlignmentResult& align_5end_result) const
{
    myutility::smith_waterman(*m_sequence, top5end_query, align_config, align_5end_result);
    if (align_5end_result.get_percent(top5end_query) > get<1>(top5end) &&
        align_5end_result.get_identity() > get<2>(top5end)
    ) {
        auto [query_align_5end, target_align_5end] = align_5end_result.get_stop_idx();
        return target_align_5end;
    }
    return 0;
}

size_t Read::trim_positive_strand_right(
    // std::string_view left_trimmed_seq_view,
    const std::string& top3end_query,
    const trim_end& top3end,
    AlignmentConfig& align_config,
    AlignmentResult& align_3end_result) const
{
    std::string left_target_sequence;
    if (m_sequence->size() > align_config.m_max_target_len) {
        left_target_sequence = m_sequence->substr(m_sequence->size()- align_config.m_max_target_len);
    } else {
        left_target_sequence = *m_sequence;
    }
    myutility::smith_waterman(left_target_sequence, top3end_query, align_config, align_3end_result);
    if (align_3end_result.get_percent(top3end_query) > get<1>(top3end) &&
        align_3end_result.get_identity() > get<2>(top3end)
    ) {
        auto [query_align_3start, target_align_3start] = align_3end_result.get_start_idx();
        return m_sequence->size() - (left_target_sequence.size() - target_align_3start);
    }
    return m_sequence->size();
}

size_t Read::trim_negative_strand_left(
    const std::string& bot5end_query,
    const trim_end& bot5end,
    AlignmentConfig& align_config,
    AlignmentResult& align_5end_result) const
{
    myutility::smith_waterman(*m_sequence, bot5end_query, align_config, align_5end_result);
    if (align_5end_result.get_percent(bot5end_query) > get<1>(bot5end) &&
        align_5end_result.get_identity() > get<2>(bot5end)
    ) {
        auto [query_align_5end, target_align_5end] = align_5end_result.get_stop_idx();
        return target_align_5end;
    }
    return 0;
}

size_t Read::trim_negative_strand_right(
    const std::string& bot3end_query,
    const trim_end& bot3end,
    AlignmentConfig& align_config,
    AlignmentResult& align_3end_result) const
{
    std::string left_target_sequence;
    if (m_sequence->size() > align_config.m_max_target_len){
        left_target_sequence = m_sequence->substr(m_sequence->size() - align_config.m_max_target_len);
    } else {
        left_target_sequence = *m_sequence;
    }
    myutility::smith_waterman(left_target_sequence, bot3end_query, align_config, align_3end_result);
    if (align_3end_result.get_percent(bot3end_query) > get<1>(bot3end) &&
        align_3end_result.get_identity() > get<2>(bot3end)
    ) {
        auto [query_align_3start, target_align_3start] = align_3end_result.get_start_idx();
        return m_sequence->size() - (left_target_sequence.size() - target_align_3start);
    }
    return m_sequence->size();
}

void Read::trim(
    const SequenceInfo& seq_info,
    const trim_direction& td,
    AlignmentConfig& align_config,
    std::ostream& log)
{
    AlignmentResult align_5end_result{true};
    AlignmentResult align_3end_result{false};
    size_t trim_start_idx{0};
    size_t trim_stop_idx{m_sequence->size()};
    std::stringstream align_string;
    if (td.trim_top5end) {
        trim_start_idx = trim_positive_strand_left(seq_info.m_top5end_query, seq_info.m_top5end, align_config,
                                                   align_5end_result);
        // align_config.to_zero();
        if (trim_start_idx != 0) {
            align_string << align_5end_result.to_string(20);
        }
        // if (m_sequence->substr(trim_start_idx, 20) != m_reads_seq->at(*m_id).substr(0, 20)){
        //     std::cout << '\n' << std::endl;
        // }
        if (td.trim_top3end) {
            trim_stop_idx = trim_positive_strand_right(seq_info.m_top3end_query, seq_info.m_top3end, align_config, align_3end_result);
            // std::string left_trimmed_seq{m_sequence->substr(trim_start_idx)};
            // trim_stop_idx = trim_positive_strand_right(left_trimmed_seq, seq_info.m_top3end_query,
            //                                            seq_info.m_top3end, align_config, align_3end_result);
        }
        if (trim_stop_idx != m_sequence->size()) {
            align_string << align_3end_result.to_string(std::get<0>(seq_info.m_top3end));
        }
    }

    if (td.trim_bot5end && trim_start_idx == 0 && trim_stop_idx == m_sequence->size()) {
        align_5end_result.to_empty();
        trim_start_idx = trim_negative_strand_left(seq_info.m_bot5end_query, seq_info.m_bot5end, align_config,
                                                   align_5end_result);
        // align_config.to_zero();
        if (trim_start_idx != 0) {
            align_string << align_5end_result.to_string(20);
        }
        if (td.trim_bot3end) {
            align_3end_result.to_empty();
            // std::string left_trimmed_seq_view{sequence_view.substr(trim_start_idx)};
            trim_stop_idx = trim_negative_strand_right(seq_info.m_bot3end_query,
                                                       seq_info.m_bot3end, align_config, align_3end_result);
            // align_config.to_zero();
        }
        if (trim_stop_idx != m_sequence->size()) {
            align_string << align_3end_result.to_string(std::get<0>(seq_info.m_top3end));
        }
    }

    if (trim_start_idx != 0 || trim_stop_idx != m_sequence->size()) {
        std::osyncstream{log} << fmt::format("{}\t{}\t{}\t{}\t{}\n{}\n",
                                             *m_id,
                                             trim_start_idx == 0 ? "None" : m_sequence->substr(0, trim_start_idx),
                                             trim_start_idx == 0 ? 0 : trim_start_idx,
                                             trim_stop_idx == m_sequence->size()
                                                 ? "None"
                                                 : m_sequence->substr(
                                                     trim_stop_idx, m_sequence->size() - trim_stop_idx),
                                             trim_stop_idx == m_sequence->size() ? 0 : m_sequence->size() - trim_stop_idx,
                                             align_string.str());
    }
    *m_sequence = m_sequence->substr(trim_start_idx, trim_stop_idx - trim_start_idx);
    *m_quality = m_quality->substr(trim_start_idx, trim_stop_idx - trim_start_idx);
}
