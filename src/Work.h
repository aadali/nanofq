#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <barrier>
#include <syncstream>
#include <FastqReader.h>

using read_stats_result = std::tuple<std::string, unsigned , double, double>;

class Work
{
private:
    FastqReader& m_fq;
    unsigned m_thread{1};
    bool m_gc;
    std::string_view m_outfile_path;
    std::ofstream m_outfile_stream;
    std::barrier<> m_bar;
//    std::vector<read_stats_result> m_stats_result{};
    std::vector<std::vector<read_stats_result>> m_sub_stats_result{};

public:
    Work() = delete;
    explicit Work(FastqReader& fq, unsigned thread, bool gc, std::string_view outfile_path);
    Work(const Work& w) = delete;
    Work(Work&& w) = delete;
    Work& operator=(const Work& w) = delete;
    Work& operator=(Work&& w) = delete;
    [[nodiscard]] std::vector<std::pair<unsigned , unsigned >> get_bins(unsigned length) const;
    void run_stats();
    void run_filter(const unsigned min_len, const unsigned max_len, const float min_quality, const float min_gc, const float max_gc);
    void run_find(const std::string& input_reads, bool use_index, unsigned key_length=5);
    void run_index(unsigned key_length);
    ~Work();
private:
    void stats(const unsigned ,
               const unsigned ,
               const shared_vec_reads&,
               std::vector<read_stats_result>& sub_stats_result);
    void filter(const unsigned start,
                const unsigned end,
                const unsigned min_len,
                const unsigned max_len,
                const float min_quality,
                const float min_gc,
                const float max_gc,
                const shared_vec_reads& reads);
};


#endif //WORK_H
