#include "submain.h"
using namespace std;

int main(int argc, char* argv[]) {
    // cout << "hello world" << endl;
    // auto a = barcode_info::get_trim_info();
    // SequenceInfo& seq_info  {a.find("SQK-LSK114")->second};
    // cout << "hello world" << endl;
    // string today{"today is a good day"};
    // cout << "hello java" << endl;
    // test_trim();
    // test_all_seq_info();
    // std::vector<string> hello {"pdf", "jpg", "mp3", "bam"};
    // myUtility::check_one_candidate("format", "pdf", hello);
    // get_arguments(argc, argv);
    std::ofstream a {"hello.txt", std::ios::out};
    std::ostream& b {a};
    sub_main(argc, argv);
}
