#include "AlignmentConfig.h"

AlignmentConfig::AlignmentConfig(int match, int mismatch, int gap_open, int gap_extend)
        : m_match(match),
          m_mismatch(mismatch),
          m_gap_open(gap_open),
          m_gap_extend(gap_extend) {
    m_score_matrix[0][0] = 0;
    m_direction_matrix[0][0] = 0;
    for (int col{1}; col < m_score_matrix[0].size(); col++) {
        m_score_matrix[0][col] = 0; // init the first row by zero
        m_direction_matrix[0][col] = 1; //  the first row of direction matrix inited by 1, one means from left
    }
    for (int row{0}; row < m_score_matrix.size(); row++) {
        m_score_matrix[row][0] = 0; // init the first column by zero
        m_direction_matrix[row][0] = 2; // the first col of direction matrix inited by 2, two means from up
    }
}
