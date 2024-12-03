#ifndef NANOBGZIP_H
#define NANOBGZIP_H

/*
 * gzip format: https://www.ietf.org/rfc/rfc1952.txt
 * bgzf format: The BGZF compression format in SAMv1.pdf https://samtools.github.io/hts-specs/SAMv1.pdf
*/


#include <iostream>
#include <fstream>
#include <filesystem>
#include <memory>
#include <sstream>
#include <zlib.h>
#include <fmt/format.h>
#include <bitset>
#include <vector>

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
constexpr size_t LINE_BUFFER{1 << 23};

namespace nanobgzip
{
    enum class GzipType { GZIP, B_GZIP, NANO_B_GZIP };


    void nano_compress(const std::string& infile,
                       const std::string& outfile,
                       const std::string& index_file,
                       int reads_number = 10,
                       unsigned key_len = 12);

    GzipType check_compress_type(const std::string& infile);

    void build_index(const std::string& file,
                    const std::string& index_file,
                     unsigned key_len);

    std::vector<uint8_t> get_uncompressed_from_block(std::ifstream& infile,
                                                     std::pair<size_t, size_t>& block_edge,
                                                     unsigned need_uncompressed_size);
}

enum class GzipType
{
    B_GZIP,
    GZIP,
    NANO_B_GZIP
};


#endif //NANOBGZIP_H
