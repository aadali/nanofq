#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <syncstream>
#include "ThreadPool.h"
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
    ThreadPool& m_threads_pool;
    std::mutex m_mtx{};

public:
    Work() = delete;
    Work(FastqReader& fq, ThreadPool& threads_pool);
    Work(const Work& w) = delete;
    Work(Work&& w) = delete;
    Work& operator=(const Work& w) = delete;
    Work& operator=(Work&& w) = delete;
    [[nodiscard]] std::vector<std::pair<unsigned, unsigned>> get_edges(int size) const;

    void run_stats(
        std::vector<read_stats_result>& stats_result,
        std::ostream& out,
        bool gc);

    void run_filter(
        std::atomic<size_t>& counter,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        std::ostream& out) const;

    void run_find(
        const std::string& input_reads,
        std::ostream& out,
        bool use_index,
        unsigned key_length) const;

    void run_index(unsigned key_length) const;

    void run_trim(
        std::atomic<size_t>& counter,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        std::vector<AlignmentConfig>& align_configs,
        std::ostream& log_fstream,
        std::ostream& out) const;

    void save_summary(
        int n,
        const std::vector<int>& read_quals,
        const std::vector<int>& read_length,
        std::vector<read_stats_result>& stats_result,
        const std::string& summary_file_path);

    ~Work() = default;

private:
    void stats(
        int start,
        int end,
        std::shared_ptr<std::vector<Read>>,
        std::vector<read_stats_result>&,
        std::ostream& out,
        bool gc);

    std::string summary_stats_result(
        int n,
        const std::vector<int>& read_quals,
        const std::vector<int>& read_lengths,
        std::vector<read_stats_result>& stats_result);

    static void stats_one_thread(
        const Read& read,
        std::vector<read_stats_result>& stats_result,
        std::ostream& out,
        bool gc);

    static void filter(
        int start,
        int end,
        std::shared_ptr<std::vector<Read>> reads_ptr,
        std::atomic<size_t>& counter,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        std::ostream& out);

    static void filter_one_thread(
        const Read& read,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        std::ostream& out);

    static void trim(
        int start,
        int end,
        std::shared_ptr<std::vector<Read>> reads_ptr,
        std::atomic<size_t>& counter,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& align_config,
        std::ostream& log_fstream,
        std::ostream& out);

    static void trim_one_thread(
        Read& read,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ostream& log_fstream,
        std::ostream& out);
};

#endif // WORK_H
