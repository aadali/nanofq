#include "Alignment.h"
#include <algorithm>
#include <ranges>
#include <sstream>
#include <fmt/core.h>

Alignment::Alignment(const std::string &target_seq) {
    m_target_seq = target_seq;
    if (m_target_seq.size() > MAX_TARGET_LEN) {
        throw std::runtime_error("the length of query sequence should not be larger than 200");
    }
}

void Alignment::align() {
    int diagonal_score{0}, up_score{0}, left_score{0};
    for (int row{1}; row < ms_query_seq5.size(); row++) {
        for (int col{1}; col < m_target_seq.size(); col++) {
            /* calculate score from diagonal */
            int score{m_target_seq[col - 1] == ms_query_seq5[row - 1] ? ms_match : ms_mismatch};
            diagonal_score = ms_score_matrix[row - 1][col - 1] + score;

            /* calculate score from left
             * 1 in direction matrix mean from left */
            score = ms_direction_matrix[row][col - 1] == 1 ? ms_gap_extend : ms_gap_open;
            left_score = ms_score_matrix[row][col - 1] + score;

            /* calculate score from up
             * 2 in direction matrix mean from up*/
            score = ms_direction_matrix[row - 1][col] == 2 ? ms_gap_extend : ms_gap_open;
            up_score = ms_score_matrix[row - 1][col] + score;

            // update the value of this cell
            ms_score_matrix[row][col] = std::max({diagonal_score, left_score, up_score});

            // record which direction the value of this cell is from
            if (ms_score_matrix[row][col] == diagonal_score) {
                ms_direction_matrix[row][col] = 0;
            } else if (ms_score_matrix[row][col] == left_score) {
                ms_direction_matrix[row][col] = 1;
            } else {
                ms_direction_matrix[row][col] = 2;
            }

            // update the max score and the row, column index of max score
            if (ms_score_matrix[row][col] > m_align_result.max_score) {
                m_align_result.max_score = ms_score_matrix[row][col];
                m_align_result.align_stop_idx = {row, col};
            }
        }
    }

    // trace back depend on the direction matrix
    auto [row, col] = m_align_result.align_stop_idx;
    while (row > 0 && col > 0 && ms_score_matrix[row][col] > 0) {
        if (ms_direction_matrix[row][col] == 0) {
            m_align_result.target_align_seq.push_back(m_target_seq[col - 1]);
            m_align_result.query_align_seq.push_back(ms_query_seq5[row - 1]);
            // '|' mean match, ':' mean mismatch
            m_align_result.line += m_target_seq[col - 1] == ms_query_seq5[row - 1] ? '|' : ':';
            --row;
            --col;
        } else if (ms_direction_matrix[row][col] == 1) {
            m_align_result.target_align_seq.push_back(m_target_seq[col - 1]);
            m_align_result.query_align_seq.push_back('-');
            //  ' ' mean gap
            m_align_result.line += ' ';
            --col;
        } else {
            m_align_result.target_align_seq.push_back('-');
            m_align_result.query_align_seq.push_back(ms_query_seq5[row - 1]);
            m_align_result.line += ' ';
            --row;
        }
    }
    std::reverse(m_align_result.target_align_seq.begin(), m_align_result.target_align_seq.end());
    std::reverse(m_align_result.query_align_seq.begin(), m_align_result.query_align_seq.end());
    std::reverse(m_align_result.line.begin(), m_align_result.line.end());
    m_align_result.align_start_idx = {row, col};
    m_align_finished = true;
    return;
}


std::string Alignment::alignment2string() const {
    if (!m_align_finished) {
        throw std::runtime_error("sequences had not been aligned");
    }
    std::string align_str;
    align_str.append(fmt::format("target: {:<4} {} {:>4}\n",
                                 m_align_result.align_start_idx.second,
                                 m_align_result.target_align_seq,
                                 m_align_result.align_stop_idx.second));
    align_str.append(fmt::format("             {}\n", m_align_result.line));
    align_str.append(fmt::format(" query: {:<4} {} {:>4}\n",
                                 m_align_result.align_start_idx.first,
                                 m_align_result.query_align_seq,
                                 m_align_result.align_stop_idx.first));
    return align_str;
}

// the following three function must be called before Alignment class used;
void Alignment::init_query(const std::string &query_seq5, const std::string& query_seq3) {
    if (!ms_init_query_seq) {
        if (query_seq5.size() > MAX_QUERY_LEN || query_seq3.size() > MAX_QUERY_LEN)
            throw std::runtime_error("the length of target sequence should not be larger than 1000");
        ms_query_seq5 = query_seq5;
        ms_query_seq3 = query_seq3;
        ms_init_query_seq = true;
    }
}

void Alignment::init_matrix() {
    if (!ms_init_matrix_finished) {
        Alignment::ms_score_matrix[0][0] = 0;
        Alignment::ms_direction_matrix[0][0] = 0;
        for (int col{1}; col < ms_score_matrix[0].size(); col++) {
            ms_score_matrix[0][col] = 0; // init the first row by zero
            ms_direction_matrix[0][col] = 1; // the first row of direction matrix inited by 1, one means from left
        }
        for (int row{0}; row < ms_score_matrix.size(); row++) {
            ms_score_matrix[row][0] = 0; // init the first column by zero
            ms_direction_matrix[row][0] = 2; // the first column of direction matrix inited by 2, two mean from up
        }
        ms_init_matrix_finished = true;
    }
}

void Alignment::init_penalty(int match, int mismatch, int gap_open, int gap_extend) {
    if (!ms_init_penalty_finished) {
        ms_match = match;
        ms_mismatch = mismatch;
        ms_gap_open = gap_open;
        ms_gap_extend = gap_extend;
        ms_init_penalty_finished = true;
    }
}

