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

enum class GzipType
{
    B_GZIP,
    GZIP,
    NANO_B_GZIP
};


struct BGZFHeader
{
    uint8_t id1; // 31
    uint8_t id2; // 139
    uint8_t cm;
    uint8_t flg; // 4 => 00000100
    uint32_t mtime;
    uint8_t xfl;
    uint8_t os;
    uint16_t xlen;  // 6
    uint8_t si1; // B
    uint8_t si2; // C
    uint16_t slen; // 2
    uint16_t bsize;
};

struct NanoBgzipHeader
{
    uint8_t id1; // 31
    uint8_t id2; // 139
    uint8_t cm;
    uint8_t flg; // 4 => 00000100
    uint32_t mtime;
    uint8_t xfl;
    uint8_t os;
    uint16_t xlen; // 8
    uint8_t si1; // N
    uint8_t si2; // A
    uint16_t slen; // 2
    uint32_t bsize;
};

class NanoBgzip
{
public:
    NanoBgzip(){};
    NanoBgzip(const NanoBgzip&) = delete;
    NanoBgzip(NanoBgzip&&) = delete;
    NanoBgzip& operator=(const NanoBgzip&) = delete;
    NanoBgzip& operator=(NanoBgzip&&) = delete;
    static void nano_compress(const std::string& infile, const std::string& outfile, int reads_number = 10);
    static GzipType check_compress_type(const std::string& infile);

    static std::vector<std::string> get_index_in_block(const std::vector<uint8_t>& input_data);
private:
    static std::vector<uint8_t> nano_block_compress(std::vector<uint8_t>& input_data, std::ostream& output_index_stream, size_t written_bytes);
};


#endif //NANOBGZIP_H
