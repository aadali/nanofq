#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <numeric>
#include <barrier>
#include <syncstream>
#include "myUtility.h"
#include "FastqReader.h"
#include "SequenceInfo.h"

#define DEFAULT_INT std::numeric_limits<int>::max()
#define DEFAULT_FLOAT 3.14f
using read_stats_result = std::tuple<std::string, unsigned, double, double>;

class Work
{
private:
    FastqReader& m_fq;
    unsigned m_thread{1};
    bool m_gc{false};
    std::string_view m_outfile_path;
    std::ofstream m_outfile_stream;
    std::barrier<> m_bar;
    //    std::vector<read_stats_result> m_stats_result{};
    std::vector<std::vector<read_stats_result>> m_sub_stats_result{};

public:
    Work() = delete;
    explicit Work(FastqReader& fq, unsigned thread, bool gc, std::string_view outfile_path);
    explicit Work(FastqReader& fq);
    Work(const Work& w) = delete;
    Work(Work&& w) = delete;
    Work& operator=(const Work& w) = delete;
    Work& operator=(Work&& w) = delete;
    [[nodiscard]] std::vector<std::pair<unsigned, unsigned>> get_bins(unsigned length) const;
    void run_stats();
    void run_filter(unsigned min_len, unsigned max_len, float min_quality, float min_gc,
                    float max_gc);
    void run_find(const std::string& input_reads, bool use_index, unsigned key_length = 5);
    void run_index(unsigned key_length) const;
    void run_trim(const SequenceInfo& seq_info, const trim_direction&td, std::vector<AlignmentConfig>& align_configs, std::fstream& log_fstream) ;
    ~Work();

private:
    void stats(unsigned,
               unsigned,
               const shared_vec_reads&,
               std::vector<read_stats_result>& sub_stats_result);
    void filter(unsigned start,
                unsigned end,
                unsigned min_len,
                unsigned max_len,
                float min_quality,
                float min_gc,
                float max_gc,
                const shared_vec_reads& reads);

    void trim(unsigned start,
              unsigned end,
              const shared_vec_reads& reads,
              const SequenceInfo& seq_info,
              const trim_direction& td,
              AlignmentConfig& alignment_config,
              std::fstream& log_fstream
              );
};


#endif //WORK_H
