#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <numeric>
#include <barrier>
#include <syncstream>
#include <span>
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
    std::ostream& m_out{std::cout};
    std::barrier<> m_bar;
    std::vector<std::vector<read_stats_result>> m_stats_result{};

public:
    Work() = delete;
    Work(FastqReader& fq, unsigned thread, bool gc, std::ostream& out);
    // Work(FastqReader& fq, unsigned thread, bool gc, std::string_view outfile_path);

    explicit Work(FastqReader& fq);
    Work(const Work& w) = delete;
    Work(Work&& w) = delete;
    Work& operator=(const Work& w) = delete;
    Work& operator=(Work&& w) = delete;
    [[nodiscard]] std::vector<std::pair<unsigned, unsigned>> get_bins(unsigned length) const;
    void run_stats();
    void run_filter(unsigned min_len, unsigned max_len, float min_quality, float min_gc,
                    float max_gc);
    void run_find(const std::string& input_reads, bool use_index, unsigned key_length = 5) const;
    void run_index(unsigned key_length) const;
    void run_trim(const SequenceInfo& seq_info, const trim_direction&td, std::vector<AlignmentConfig>& align_configs, std::ostream& log_fstream) ;
    void save_summary(int, std::vector<int>&, const std::string& summary_file_path);
    ~Work()= default;

private:
    void stats(unsigned,
               unsigned,
               const shared_vec_reads&,
               std::vector<read_stats_result>& sub_stats_result);
    std::string summary_stats_result(int n, std::vector<int>& read_quals);
    void stats_one_thread(const Read &read);
    void filter(unsigned start,
                unsigned end,
                unsigned min_len,
                unsigned max_len,
                float min_quality,
                float min_gc,
                float max_gc,
                const shared_vec_reads& reads);
    void filter_one_thread(const Read& read, unsigned min_len, unsigned max_len, float min_quality, float min_gc, float max_gc) const;

    void trim(unsigned start,
              unsigned end,
              const shared_vec_reads& reads,
              const SequenceInfo& seq_info,
              const trim_direction& td,
              AlignmentConfig& align_config,
              std::ostream& log_fstream
              );
    void trim_one_thread(Read& read, const SequenceInfo& seq_info, const trim_direction& td, AlignmentConfig& alignment_config, std::ostream& log_fstream) const;
};


#endif //WORK_H
