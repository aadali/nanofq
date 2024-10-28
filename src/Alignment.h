#ifndef NANOFQ_ALIGNMENT_H
#define NANOFQ_ALIGNMENT_H


#include <string_view>
#include <string>
#include <array>
#include <vector>
#include <stdexcept>

#define MAX_TARGET_LEN 2000
#define MAX_QUERY_LEN 200

/**
 * @struct AlignmentResult
 * @brief Stores the result of the alignment process.
 */
struct AlignmentResult {
    std::string target_align_seq;
    std::string line;   // line between the two aligned sequence. '|' mean match, ':' mean mismatch, ' ' mean gap
    std::string query_align_seq;
    std::pair<size_t, size_t> align_stop_idx{0, 0}; // <query_align_stop_idx, target_align_stop_idx>
    std::pair<size_t, size_t> align_start_idx{0, 0};// <query_align_start_idx, target_align_start_idx>
    int max_score{0};
};

/**
 * @class Alignment
 * @brief Implements the Smith-Waterman algorithm for local sequence alignment.
 */
class Alignment {
private:
    // target_seq: the nanopore read
    std::string_view m_target_seq;
    bool m_align_finished{false};
    AlignmentResult m_align_result;
private:
    static inline bool ms_init_query_seq{false};
    static inline bool ms_init_penalty_finished{false};
    static inline bool ms_init_matrix_finished{false};
    /*
     * A simple Nanopore read structure, especially for LSK, in most case, seq3 is truncated.
     * For Rapid kit, always just consider the seq5 and ignore seq3
     *
     * <------seq5------->                                     <-------seq3------->
     * 5'-ATCGGCGGAGATTTCA.....................................TCGAGCATGACGACAGTGAG-3'
     * 3'-TAGCCGCCTCTAAAGT.....................................AGCTCGTACTGCTGTCACTC-5'
     * <-------seq3------>                                     <-------seq5------->

     query_seq: adapter, barcode etc.
     */
    static inline std::string ms_query_seq5;
    static inline std::string ms_query_seq3;
    static inline int ms_match; // default 3
    static inline int ms_mismatch; // default -3
    static inline int ms_gap_open; // default -10
    static inline int ms_gap_extend; // default -2
    static inline int ms_max_target_length;

    // the adapter or barcode or primer, even sum of these values would never be larger than 200
    static inline int ms_max_query_length{200};

    // fastq read(target) as column; barcode or adapter or primer(query) as row;
    // the length of search range in read should not larger than 1000;
    static inline std::array<std::array<int, MAX_TARGET_LEN + 1>, MAX_QUERY_LEN + 1> ms_score_matrix{0};

    // 0: from diagonal; 1: from left; 2 from up
    static inline std::array<std::array<int, MAX_TARGET_LEN + 1>, MAX_QUERY_LEN + 1> ms_direction_matrix{0};
public:
    // Construct
    explicit Alignment(const std::string &query_seq);

    Alignment() = delete;

    Alignment(Alignment &&src) = delete;

    Alignment(const Alignment &src) = delete;

    Alignment &operator=(const Alignment &rhs) = delete;

    Alignment &operator=(Alignment &&rhs) = delete;

public:
    void align();

    [[nodiscard]]inline int get_alignment_score() const {
        if (!m_align_finished) {
            throw std::runtime_error("sequences had not been aligned");
        }
        return m_align_result.max_score;
    };

    /**
    * @struct AlignmentResult
    * @brief Stores the result of the alignment process.
    */
    [[nodiscard]] std::string alignment2string() const;

    static void init_query(const std::string &query_seq5, const std::string &qyery_seq3);

    /**
     * @brief Initializes the penalty values.
     * @param match The match score.
     * @param mismatch The mismatch penalty.
     * @param gap_open The gap opening penalty.
     * @param gap_extend The gap extension penalty.
     */
    static void init_penalty(int match, int mismatch, int gap_open, int gap_extend);


    static void init_matrix();


};


#endif //NANOFQ_ALIGNMENT_H
