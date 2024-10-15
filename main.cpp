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
    // big_big_fastq = "/home/a/pub/ycq/data/ont-data-release/londoncalling2024/assembly/basecalling/ulk/PAW42495.fastq";
    fstream infile{big_fastq.data(), ios::in};
    fstream outfile{"../index.txt", ios::out};
    size_t buffer_size{1 << 23};
    char *line{new char[buffer_size]};
    int line_number = 1;
    size_t start{0};
    size_t stop{0};
    string id;
    infile.seekg(878, ios::beg);
    for (int i{0}; i<15682-878+1; i++){
        cout << static_cast<char>(infile.get());
    }
    cout << '\n';
    while (infile.getline(line, buffer_size, '\n')) {
        break;
        switch (line_number) {
            case 1: {
                for (int idx{0}; idx < strlen(line); idx++) {
                    if (line[idx] == ' ') {
                        id.assign(line, 1, idx - 1);
                        break;
                    }
                }
                if (id.empty()) id = line;
                stop += strlen(line) + 1;
                line_number++;
                break;
            }
            case 2: {
                stop += strlen(line) + 1;
                line_number++;
                break;
            }
            case 3: {
                stop += strlen(line) + 1;
                line_number++;
                break;
            }
            case 4: {
                stop += strlen(line);
                line_number = 1;
                cout << id << '\t' << start << '\t' << stop << endl;
                start = stop + 1;
                break;
            }
            default: {
                break;
            }
        }
    }

    infile.close();
    outfile.close();
    delete[] line;
//
//    FastqReader fq{big_fastq, 50000};
//    Work work{fq, 4, false, "../output.txt"};
////    work.run_find("834aca2a-ca33-4bd7-b98c-0c7653dd7877,a442602c-5355-45a7-be29-11ef83f201cf");
//    thread t1{&FastqReader::read_chunk_fastq, &fq};
//    thread t2{&Work::run_filter, &work, 9000, 10000, 12, 0.0, 1.0};
//    t1.join();
//    t2.join();
    // work.run_stats();
}
