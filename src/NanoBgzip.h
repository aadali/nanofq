//
// Created by a on 24-11-27.
//

#ifndef NANOBGZIP_H
#define NANOBGZIP_H

#include <string_view>
#include <iostream>
#include <fstream>
#include <memory>
#include <zlib.h>
#include <fmt/format.h>
#include <bitset>
#include <vector>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
enum class GzipType {
    B_GZIP, //
    GZIP,
    NANO_B_GZIP
};
enum class OperationType{
    COMPRESS,
    DECOMPRESS
};

class NanoBgzip {
private:
    std::string m_infile{"-"};

public:
    explicit NanoBgzip(const std::string& infile);
    void compress();
    //    void decompress();
    //    void index();
    //private:
    GzipType check_compress_type() const;
private:
    std::vector<uint8_t> compress(std::vector<uint8_t>& input_data);
};


#endif //NANOBGZIP_H
