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
#include <atomic>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

using namespace std;


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
//    fstream infile{big_big_fastq.data(), ios::in};
//    fstream outfile{"../index.txt", ios::out};
//    size_t buffer_size{1 << 23};
//    char *line{new char[buffer_size]};
//    int line_number = 1;
//    size_t start{0};
//    size_t stop{0};
//    string id;
//    infile.seekg(ios::beg);
//    for (int i{0}; i<878; i++){
//        cout << static_cast<char>(infile.get());
//    }
//    infile.seekg(47182, ios::beg);
//    for(int i{47182};i<53513; i++){
//        cout << static_cast<char>(infile.get());
//    }
//    infile.seekg(878, ios::beg);
//    for (int i{0}; i<15682-878+1; i++){
//        cout << static_cast<char>(infile.get());
//    }
//    cout << '\n';
//    while (infile.getline(line, buffer_size, '\n')) {
//        break;
//        switch (line_number) {
//            case 1: {
//                for (int idx{0}; idx < strlen(line); idx++) {
//                    if (line[idx] == ' ') {
//                        id.assign(line, 1, idx - 1);
//                        break;
//                    }
//                }
//                if (id.empty()) id = line;
//                stop += strlen(line) + 1;
//                line_number++;
//                break;
//            }
//            case 2: {
//                stop += strlen(line) + 1;
//                line_number++;
//                break;
//            }
//            case 3: {
//                stop += strlen(line) + 1;
//                line_number++;
//                break;
//            }
//            case 4: {
//                stop += strlen(line) + 1;
//                line_number = 1;
//                cout << id << '\t' << start << '\t' << stop << endl;
//                start = stop ;
//                break;
//            }
//            default: {
//                break;
//            }
//        }
//    }
//
//    infile.close();
//    outfile.close();
//    delete[] line;
//
    FastqReader fq{big_big_fastq, 50000};
    Work work{fq, 4, false, "../test_data/output.txt"};
    work.run_index();
//    {
//        Timer time2{"cereal;"};
//        std::fstream infile{"/home/a/pub/ycq/data/ont-data-release/londoncalling2024/assembly/basecalling/ulk/PAW42495.fastq.idx", std::ios::in|ios::binary};
//        cereal::BinaryInputArchive archive{infile};
//        unordered_map<string, tuple<size_t, size_t>> read_index;
//        cereal::load(archive, read_index);
//    }

//    work.run_find("024ec006-d76e-4fb0-9f33-366997ae24b0,9ae267ae-9375-4a86-b87e-f49d9c2367e7", false);
//    work.run_index();
//    work.run_stats();
//    work.run_find("834aca2a-ca33-4bd7-b98c-0c7653dd7877,a442602c-5355-45a7-be29-11ef83f201cf,49f8fb43-50e2-4408-945f-755237825c7f,4a8d732c-6f91-4b75-8d99-95c692851718",
//                  false);
//    thread t1{&FastqReader::read_chunk_fastq, &fq};
//    thread t2{&Work::run_stats, &work};
//    thread t2{&Work::run_filter, &work, 9000, 10000, 12, 0.0, 1.0};
//    t1.join();
//    t2.join();
    // work.run_stats();
}