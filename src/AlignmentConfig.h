#ifndef NANOFQ_ALIGNMENTCONFIG_H
#define NANOFQ_ALIGNMENTCONFIG_H
// #define MAX_TARGET_LEN 2000
// #define MAX_QUERY_LEN 200

#include <array>
#include <vector>

class AlignmentConfig {
public:
    int m_max_target_len;
    int m_max_query_len;
    int m_match;
    int m_mismatch;
    int m_gap_open;
    int m_gap_extend;
    enum class Direction {Diag = 0, Left =1, Up=2 };
private:
    std::vector<std::vector<int>> m_score_matrix;
    std::vector<std::vector<Direction>> m_direction_matrix;
    // std::array<std::array<int, MAX_TARGET_LEN + 1>, MAX_QUERY_LEN + 1> m_score_matrix{};
    // std::array<std::array<Direction, MAX_TARGET_LEN + 1>, MAX_QUERY_LEN + 1> m_direction_matrix{};
public:
    AlignmentConfig(int max_target_len, int max_query_len, int match, int mismatch, int gap_open, int gap_extend);

    AlignmentConfig() = delete;

    AlignmentConfig(AlignmentConfig &&) = delete;

    AlignmentConfig(const AlignmentConfig &);

    AlignmentConfig &operator=(AlignmentConfig &&) = delete;

    AlignmentConfig &operator=(const AlignmentConfig &) = delete;

    ~AlignmentConfig() = default;

public:
    inline int get_score(int row, int col) const {
        return m_score_matrix[row][col];
    }

    inline Direction get_direction(int row, int col) const {
        return m_direction_matrix[row][col];
    }

    inline void set_score(int row, int col, int score) {
        m_score_matrix[row][col] = score;
    }

    inline void set_direction(int row, int col, Direction direction) {
        m_direction_matrix[row][col] = direction;
    }
};


#endif //NANOFQ_ALIGNMENTCONFIG_H
