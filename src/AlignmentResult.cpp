#include "AlignmentResult.h"
#include <algorithm>
#include <fmt/core.h>

std::string AlignmentResult::to_string() {
    std::string align_str;
    align_str.append(fmt::format("target: {:<4} {} {:>4}\n",
                                 m_align_start_idx.second,
                                 m_target_align_seq,
                                 m_align_stop_idx.second));
    align_str.append(fmt::format("             {}\n", m_line));
    align_str.append(fmt::format(" query: {:<4} {} {:>4}\n",
                                 m_align_start_idx.first,
                                 m_query_align_seq,
                                 m_align_stop_idx.first));
    return align_str;
}

void AlignmentResult::reverse_align() {
    std::reverse(m_target_align_seq.begin(), m_target_align_seq.end());
    std::reverse(m_line.begin(), m_line.end());
    std::reverse(m_query_align_seq.begin(), m_query_align_seq.end());
}

bool AlignmentResult::is_empty() const {
    return m_query_align_seq.empty() &&
           m_target_align_seq.empty() &&
           m_line.empty() &&
           m_max_score == 0 &&
           m_align_start_idx.first == 0 &&
           m_align_start_idx.second == 0 &&
           m_align_stop_idx.first == 0 &&
           m_align_stop_idx.second == 0;
}