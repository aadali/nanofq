#include "AlignmentResult.h"
#include <algorithm>
#include <iostream>
#include <fmt/core.h>

std::string AlignmentResult::to_string(size_t target_3end_len) {
    std::string align_str;
    if (!m_left) {
        align_str.append(fmt::format("target: {:<4} {} {:>4}\n",
                                     static_cast<int>(m_align_start_idx.second) - static_cast<int>(target_3end_len),
                                     m_target_align_seq,
                                     static_cast<int>(m_align_stop_idx.second) - static_cast<int>(target_3end_len)));
    } else {
        align_str.append(fmt::format("target: {:<4} {} {:>4}\n",
                                     m_align_start_idx.second,
                                     m_target_align_seq,
                                     m_align_stop_idx.second));
    }
    align_str.append(fmt::format("             {}\n", m_line));
    align_str.append(fmt::format(" query: {:<4} {} {:>4}\n",
                                 m_align_start_idx.first,
                                 m_query_align_seq,
                                 m_align_stop_idx.first));
    // std::cout << m_target_align_seq << '\n';
    // std::cout << m_query_align_seq << '\n';
    // std::cout << align_str<<'\n';
    return align_str;
}

void AlignmentResult::reverse_align() {
    std::ranges::reverse(m_target_align_seq);
    std::ranges::reverse(m_line);
    std::ranges::reverse(m_query_align_seq);
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
