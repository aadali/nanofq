#include "Timer.h"
#include "submain.h"
#include "NanoBgzip.h"
using namespace std;


int main(int argc, char* argv[])
{
    NanoBgzip nb("-");
    nb.compress();
    // NanoBgzip nb("/home/a/big/ycq/projects/CppProjects/nanofq/test_data/output.fastq.gz");
    // nb.check_compress_type();
    // std::string fastq{
    //     "/home/a/nieyueqing.methy.pass.fastq"
    // };
    // auto start{std::chrono::steady_clock::now()};
    // std::vector<read_stats_result> stats_result{};
    // FastqReader fr{fastq, 20000};
    // ThreadPool threads_pool{1, stats_result};
    // Work worker{fr, threads_pool};
    // worker.run_stats(stats_result, std::cout, false);
    // std::vector<int> quals{25, 20, 18, 15, 12, 10};
    // std::vector<int> lengths{9000, 2000};
    // worker.save_summary(10, quals, lengths, stats_result, "./summary2.txt");
    // // std::this_thread::sleep_for(std::chrono::milliseconds(100ms));
    // auto end{std::chrono::steady_clock::now()};
    // auto during{std::chrono::duration_cast<std::chrono::milliseconds>(end - start)};
    // cout << "stats time: " << during.count() << "ms" << endl;
    // cout << stats_result.size() << endl;
}
