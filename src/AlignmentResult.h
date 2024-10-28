#ifndef NANOFQ_ALIGNMENTRESULT_H
#define NANOFQ_ALIGNMENTRESULT_H

#include <string>

class AlignmentResult {
private:
    std::string m_target_align_seq;
    std::string m_line; // line between the two aligned sequence. '|' mean match, ':' mean mismatch, ' ' mean gap
    std::string m_query_align_seq;
    std::pair<size_t, size_t> m_align_stop_idx;  // <query_align_stop_idx, target_align_stop_idx>
    std::pair<size_t, size_t> m_align_start_idx;// <query_align_start_idx, target_align_start_idx>;
    int m_max_score{0};
public:
    AlignmentResult() {};

    AlignmentResult(AlignmentResult &&) = delete;

    AlignmentResult(const AlignmentResult &) = delete;

    AlignmentResult &operator=(AlignmentResult &&) = delete;

    AlignmentResult &operator=(const AlignmentResult &) = delete;

    ~AlignmentResult() = default;

public:
    inline void push_back(char target_c, char query_c, char line_c) {
        m_target_align_seq.push_back(target_c);
        m_query_align_seq.push_back(query_c);
        m_line.push_back(line_c);
    }


    inline void set_start_idx(int x, int y) { m_align_start_idx = {x, y}; }

    inline void set_stop_idx(int x, int y) { m_align_stop_idx = {x, y}; }

    inline std::pair<size_t, size_t> &get_start_idx() { return m_align_start_idx; }

    inline std::pair<size_t, size_t> &get_stop_idx() { return m_align_stop_idx; }

    inline int get_max_score() const { return m_max_score; }

    inline void set_max_score(int mMaxScore) { m_max_score = mMaxScore; }

    bool is_empty() const;

    std::string to_string();

    void reverse_align();

};


#endif //NANOFQ_ALIGNMENTRESULT_H
