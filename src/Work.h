#ifndef WORK_H
#define WORK_H
#include <tuple>
#include <string>
#include <syncstream>
#include <variant>
#include <barrier>

#include "ThreadPool.h"
#include "myUtility.h"
#include "FastqReader.h"
#include "SequenceInfo.h"

#define DEFAULT_INT std::numeric_limits<int>::max()
#define DEFAULT_FLOAT 3.14f
using read_stats_result = std::tuple<std::shared_ptr<std::string>, unsigned, double, double>;
// using read_stats_result = std::tuple<std::shared_ptr<std::string>, unsigned, double, double, unsigned, double, double, bool>;

class Work {
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

    void run_stats_multi_fqs_in_multi_threads(
        const std::vector<std::filesystem::path>& paths,
        std::vector<read_stats_result>& stats_result,
        std::ostream& out,
        bool gc);

    void run_main(
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        std::vector<AlignmentConfig>& align_configs,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        std::mutex& all_mtx,
        std::mutex& passed_mtx,
        std::barrier<>& bar
    );

    void run_main_multi_fqs_in_multi_threads(
        const std::vector<std::filesystem::path>& paths,
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        std::vector<AlignmentConfig>& align_configs,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        std::mutex& all_mtx,
        std::mutex& passed_mtx
    );

    void run_filter(
        std::atomic<size_t>& counter,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        std::ostream& out) const;

    void run_filter_multi_fqs_in_multi_threads(
        const std::vector<std::filesystem::path>& paths,
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
        bool use_index) const;

    void run_index(bool force_index) const;

    void run_trim(
        std::atomic<size_t>& counter,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        std::vector<AlignmentConfig>& align_configs,
        std::ostream& log_fstream,
        std::ostream& out, std::barrier<>&bar) const;

    void run_trim_multi_fqs_in_multi_threads(
        const std::vector<std::filesystem::path>& paths,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        std::vector<AlignmentConfig>& align_configs,
        std::ostream& log_fstream,
        std::ostream& out) const;

    std::tuple<float, int, float, float> save_summary(
        int n,
        const std::vector<int>& read_quals,
        const std::vector<int>& read_length,
        std::vector<read_stats_result>& reads_stats_result,
        const std::string& summary_file_path,
        bool is_passed
        );

    void plot(
        const std::string& argv0,
        const std::string& input,
        const std::string& prefix,
        bool plot_mean_length,
        float mean_length,
        bool plot_n50,
        int n50,
        float std,
        const std::vector<std::string>& fmt,
        float mean_quality);

    ~Work() = default;

private:
    void stats(
        int start,
        int end,
        std::shared_ptr<std::vector<Read>>,
        std::vector<read_stats_result>&,
        std::ostream& out,
        bool gc );

    std::tuple<std::string, float, int, float, float> summary_stats_result(
        int n,
        const std::vector<int>& read_quals,
        const std::vector<int>& read_lengths,
        std::vector<read_stats_result>& stats_results_vec
        );

    // std::tuple<std::string, float, int, float, float> main_summary_stats_result(
    //     int n,
    //     const std::vector<int>& read_quals,
    //     const std::vector<int>& read_lengths,
    //     std::vector<read_stats_result>& stats_result_vec
    // );

    static void stats_one_thread(
        const Read& read,
        std::vector<read_stats_result>& stats_result,
        std::ostream& out,
        bool gc);

    void stats_multi_fqs_in_one_thread(
        const std::vector<std::filesystem::path>& paths,
        size_t start,
        size_t end,
        std::vector<read_stats_result>* stats_result,
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
        std::ostream& out );

    static void filter_one_thread(
        const Read& read,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        std::ostream& out);

    static void filter_multi_fqs_in_one_thread(
        const std::vector<std::filesystem::path>& paths,
        size_t start,
        size_t end,
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
        std::ostream& out,
        std::barrier<>& bar
        );

    static void trim_one_thread(
        Read& read,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ostream& log_fstream,
        std::ostream& out);

    static void trim_multi_fqs_in_one_thread(
        const std::vector<std::filesystem::path>& paths,
        size_t start,
        size_t end,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ostream& log_fstream,
        std::ostream& out);

    void main(
        int start,
        int end,
        std::shared_ptr<std::vector<Read>> reads_ptr,
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        std::mutex& all_mtx,
        std::mutex& passed_mtx,
        std::barrier<> &bar
    );

    void main_one_thread(
        Read& read,
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        std::mutex& all_mtx,
        std::mutex& passed_mtx
    );

    void main_multi_fqs_in_one_thread(
        const std::vector<std::filesystem::path>& paths,
        size_t start,
        size_t end,
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        std::mutex& all_mtx,
        std::mutex& passed_mtx
    );

    void main_core(
        Read& read,
        std::vector<read_stats_result>& all_reads_stats_result,
        std::vector<read_stats_result>& passed_reads_stats_result,
        bool gc,
        unsigned min_len,
        unsigned max_len,
        float min_quality,
        float min_gc,
        float max_gc,
        bool do_trim,
        const SequenceInfo& seq_info,
        const trim_direction& td,
        AlignmentConfig& alignment_config,
        std::ofstream& out_ofstream,
        std::ofstream& all_stats_ofstream,
        std::ofstream& passed_stats_ofstream,
        std::ofstream& failed_stats_ofstream,
        std::ofstream& trim_log_ofstream,
        bool retain_failed,
        std::ofstream& failed_ofstream,
        bool sync,
        std::mutex& all_mtx,
        std::mutex& passed_mtx
        );
};

#endif // WORK_H
