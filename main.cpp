#include <cstring>
#include <fstream>
#include <fmt/core.h>
#include <execution>
#include "Read.hpp"
#include "Timer.hpp"
#include "FastqReader.h"
#include <thread>
#include <syncstream>
#include "Work.h"
//#include "utility.h"
#include <atomic>


#include <cereal/types/unordered_map.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include "Alignment.h"

using namespace std;
//using namespace utility::utility;


int main() {
    Timer timer{"test"};
    // string_view fastq{"/home/a/big/ycq/projects/CppProjects/NanoFq/test_data/test.fastq"};
    string_view big_fastq{
            "/home/a/pub/ycq/data/20240610-zdyfy-FSHD-fast5s/fast5s/20230619-zdyfy-nieyueqing-fast5/sub_nieyueqing.methy.pass.fastq"
    };
    string_view big_big_fastq{
            "/home/a/pub/ycq/data/20240610-zdyfy-FSHD-fast5s/fast5s/20230619-zdyfy-nieyueqing-fast5/nieyueqing.methy.pass.fastq"
    };
    big_big_fastq = "/home/a/pub/ycq/data/ont-data-release/londoncalling2024/assembly/basecalling/ulk/PAW42495.fastq";
    FastqReader fq{big_big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/output.txt"};
//    const string hello{"hello world python java rust  "};
//    auto a = split(hello, ' ');
//    for (auto i: a) {
//        cout << i << ";" << endl;
//    }
//    cout << "hello world" << endl;
//    work.run_index(8);
//    work.run_find("a3f6101d-ba31-4654-accb-40ec9ff451b6,72f8efeb-a61c-4714-bd3c-24ff18301ffa,ff4c91d6-0490-497c-9150-3fbbe3e6b958",
//                  true);
//    work.run_find("d04fb7dd-eec6-49a4-8a76-51b72bcb8e8e,7ce48120-c654-48b3-91f6-f6ba5d400685,af974c5a-dec6-45cc-bf93-6036525e4141",
//                  true);
//    thread t1{&FastqReader::read_chunk_fastq, &fq};
//    thread t2{&Work::run_stats, &work};
//    thread t2{&Work::run_filter, &work, 9000, 10000, 12, 0.0, 1.0};
//    t1.join();
//    t2.join();
//    string hello{"hello world"};
    string target {"ATGTGTATATTTATAGCTTCCATTTATTCAAAAACCGGGTATTTTTCCAACCAAGAAAGTTGTCGGTGTCTTTGTGGTTTTCACATTATCGTGAAACGCTTTCAGCATTTTCAGCTACACACTTTCACATTTCCCATCTTCTGGCTTGTTTAAAAGCTCTAGACACAGCCAAGCACAGTGTGTATAAGTGCCTCCTCAGTGCTGGTACTCAGCCTATCCCAATATTGG"};
    string query {"AAGAAAGTTGTCGGTGTCTTTGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"};
    Alignment sw {target};
    Alignment::init_query(query, query);
    Alignment::init_matrix();
    Alignment::init_penalty(3, -3, -10, -1);
    sw.align();
    auto res = sw.alignment2string();
    cout << res << endl;
}
