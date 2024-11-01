#include "Timer.h"
#include "FastqReader.h"
#include <syncstream>

#include "Adapter.h"
#include "Work.h"
#include "myUtility.h"
#include "AlignmentConfig.h"
#include "AlignmentResult.h"


using namespace std;
string_view normal_fastq{"../test_data/normal.fastq"};
string_view normal_fastq_gz{"../test_data/normal.fastq.gz"};
string_view big_fastq{"../test_data/big.fastq"};
string_view big_fastq_gz{"../test_data/big.fastq.gz"};
string_view big_big_fastq{"../test_data/big_big.fastq"};
string_view big_big_fastq_gz{"../test_data/big_big.fastq.gz"};

//using namespace utility::utility;
void test_index() {
    Timer timer{"test index"};
    FastqReader fq{big_big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/test_output/output.txt"};
    work.run_index(8);
}

void test_find(bool use_index) {
    Timer timer{"test find"};
    FastqReader fq{big_big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/test_output/test_find.fastq"};
    work.run_find(
        "a3f6101d-ba31-4654-accb-40ec9ff451b6,72f8efeb-a61c-4714-bd3c-24ff18301ffa,ff4c91d6-0490-497c-9150-3fbbe3e6b958",
        use_index);
}

void test_stats() {
    // TODO stats bug, output number is different from input
    Timer timer{"test stats"};
    //    FastqReader fq{big_big_fastq, 100000};
    FastqReader fq{big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/test_output/test_stats.txt"};
    thread t1{&FastqReader::read_chunk_fastq, &fq};
    thread t2{&Work::run_stats, &work};
    t1.join();
    t2.join();
}

void test_filter() {
    Timer timer{"test filter"};
    FastqReader fq{big_big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/test_output/test_filter.fastq"};
    thread t1{&FastqReader::read_chunk_fastq, &fq};
    thread t2{&Work::run_filter, &work, 9000, 10000, 12, 0.0, 1.0};
    t1.join();
    t2.join();
}

void test_smith_waterman() {
    Timer timer{"test Alignment"};
    string target{
        "ATGTGTATATTTATAGCTTCCATTTATTCAAAAACCGGGTATTTTTCCAACCAAGAAAGTTGTCGGTGTCTTTGTGGTTTTCACATTATCGTGAAACGCTTTCAGCATTTTCAGCTACACACTTTCACATTTCCCATCTTCTGGCTTGTTTAAAAGCTCTAGACACAGCCAAGCACAGTGTGTATAAGTGCCTCCTCAGTGCTGGTACTCAGCCTATCCCAATATTGG"
    };
    string query{"AAGAAAGTTGTCGGTGTCTTTGTGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"};
    AlignmentConfig config{3, -3, -10, -1};
    AlignmentResult result{};
    myUtility::smith_waterman(target, query, config, result);
    cout << result.to_string() << endl;
}

void test_trim(){
    Timer timer{"test Trim"};
    auto trim_info = barcode_info::get_trim_info();
    FastqReader fq{"/home/a/big/ycq/projects/CppProjects/nanofq/test_data/nbd114.24/barcode01.fastq", 1000};
    SequenceInfo& sequence_info {trim_info.find("SQK-NBD114.24-1")->second};
    trim_direction td {myUtility::how_trim(sequence_info)};
    AlignmentConfig align_config{3, -3, -12, -1};
    Work work{fq, 1, false, "../test_data/test_output/test_trim.fastq"};
    std::fstream outfile {"../test_data/test_output/trim.log", std::ios::out};
    thread t1{&FastqReader::read_chunk_fastq, &fq};
    thread t2{&Work::run_trim, &work, std::ref(sequence_info), std::ref(td), std::ref(align_config), std::ref(outfile)};
    t1.join();
    t2.join();
}
