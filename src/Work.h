#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <barrier>
#include <syncstream>
#include <FastqReader.h>

using read_stats_result = std::tuple<std::string, size_t, double, double>;

class Work
{
private:
    FastqReader& m_fq;
    size_t m_thread{1};
    bool m_gc;
    std::string_view m_outfile_path;
    std::ofstream m_outfile_stream;
    std::barrier<> m_bar;
//    std::vector<read_stats_result> m_stats_result{};
    std::vector<std::vector<read_stats_result>> m_sub_stats_result{};

public:
    Work() = delete;
    explicit Work(FastqReader& fq, size_t thread, bool gc, std::string_view outfile_path);
    Work(const Work& w) = delete;
    Work(Work&& w) = delete;
    Work& operator=(const Work& w) = delete;
    Work& operator=(Work&& w) = delete;
    [[nodiscard]] std::vector<std::pair<size_t, size_t>> get_bins(size_t length) const;
    void run_stats();
    void run_filter(const size_t min_len, const size_t max_len, const double min_quality, const double min_gc, const double max_gc);
    void run_find(const std::string& input_reads);
    ~Work();
private:
    void stats(const size_t,
               const size_t,
               const shared_vec_reads&,
               std::vector<read_stats_result>& sub_stats_result);
    void filter(const size_t start,
                const size_t end,
                const size_t min_len,
                const size_t max_len,
                const double min_quality,
                const double min_gc,
                const double max_gc,
                const shared_vec_reads& reads);
};


#endif //WORK_H
