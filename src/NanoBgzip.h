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
constexpr size_t LINE_BUFFER{1<<23};

namespace nanobgzip
{
    enum class GzipType{GZIP, B_GZIP, NANO_B_GZIP};

    struct BGZFHeader;

    struct NanoBgzipHeader;

    void nano_compress(const std::string& infile, const std::string& outfile, int reads_number = 10);
    GzipType check_compress_type(const std::string& infile);
    std::vector<std::string> get_index_in_block(const std::vector<uint8_t>& input_data);
    void build_index(const std::string& file);
    std::vector<uint8_t> nano_block_compress(std::shared_ptr<std::vector<uint8_t>> input_data,
                                             std::ostream& output_index_stream,
                                             size_t written_bytes);
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
