#include "AlignmentConfig.h"

AlignmentConfig::AlignmentConfig(int max_target_len, int max_query_len, int match, int mismatch, int gap_open,
                                 int gap_extend)

    : m_max_target_len(max_target_len),
    m_max_query_len(max_query_len),
    m_match(match),
      m_mismatch(mismatch),
      m_gap_open(gap_open),
      m_gap_extend(gap_extend),
      m_score_matrix(max_query_len + 1),
      m_direction_matrix(max_query_len + 1)
{
    for (int i{0}; i < max_query_len + 1; ++i) {
        m_score_matrix[i].resize(m_max_target_len + 1, 0);
        m_direction_matrix[i].resize(m_max_target_len + 1, Direction::Diag);
    }
    // m_score_matrix[0][0] = 0;
    // m_direction_matrix[0][0] = Direction::Diag;
    for (int col{1}; col < m_score_matrix[0].size(); col++) {
        m_score_matrix[0][col] = 0; // init the first row by zero
        m_direction_matrix[0][col] = Direction::Left;
        //  the first row of direction matrix inited by 1, one means from left
    }
    for (int row{0}; row < m_score_matrix.size(); row++) {
        m_score_matrix[row][0] = 0; // init the first column by zero
        m_direction_matrix[row][0] = Direction::Up; // the first col of direction matrix inited by 2, two means from up
    }
}

AlignmentConfig::AlignmentConfig(const AlignmentConfig& src)
{
    m_max_target_len = src.m_max_target_len;
    m_max_query_len = src.m_max_query_len;
    m_match = src.m_match;
    m_mismatch = src.m_mismatch;
    m_gap_open = src.m_gap_open;
    m_gap_extend = src.m_gap_extend;
    m_score_matrix = src.m_score_matrix;
    m_direction_matrix = src.m_direction_matrix;
}
