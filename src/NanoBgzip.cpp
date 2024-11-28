#include "NanoBgzip.h"

#include <vector>


NanoBgzip::NanoBgzip(const std::string& infile) : m_infile(infile) {}

void NanoBgzip::compress()
{
    // std::unique_ptr<std::istream> in{&std::cin};
    // if (m_infile != "-") in = std::make_unique<std::istream>(std::ifstream{m_infile});
    std::vector<uint8_t> in_buffer;
    char ch;
    int newline_count{0};
    std::ifstream infile{"/home/a/big/ycq/projects/CppProjects/nanofq/test_data/first4line.fastq"};

    while (infile.get(ch)) {
        in_buffer.push_back(ch);
        if (ch == '\n') {
            ++newline_count;
            if (newline_count == 4) {
                std::vector<uint8_t> output_buffer{compress(in_buffer)};
                cout.write(reinterpret_cast<char*>(output_buffer.data()), output_buffer.size());
            }
        }
    }
}

GzipType NanoBgzip::check_compress_type() const
{
    // https://www.ietf.org/rfc/rfc1952.txt
    std::ifstream istream{m_infile, std::ios::binary};
    uint8_t id1{0};
    uint8_t id2{0};
    istream.read(reinterpret_cast<char*>(&id1), 1);
    if (istream.gcount() == 0) {
        cerr << "File is empty" << endl;
        exit(1);
    }
    istream.read(reinterpret_cast<char*>(&id2), 1);
    if (id1 != 0x1f || id2 != 0x8b) {}
    uint8_t cm{0};
    istream.read(reinterpret_cast<char*>(&cm), 1); // get cm, compression method should be 8(deflate)
    if (cm != 8) {
        //
        cerr << "Not a valid gzip file" << endl;
        exit(1);
    }
    uint8_t flag{0};
    istream.read(reinterpret_cast<char*>(&flag), 1); // get flag
    istream.ignore(4); // ignore mtime
    istream.ignore(2); // ignore xfl and os
    auto flag_bits{std::bitset<8>(flag)};
    if (flag_bits[2]) {
        uint16_t xlen{0};
        istream.read(reinterpret_cast<char*>(&xlen), 2);
        uint8_t sid1{0};
        uint8_t sid2{0};
        uint16_t sub_fields_len{0};
        istream.read(reinterpret_cast<char*>(&sid1), 1);
        istream.read(reinterpret_cast<char*>(&sid2), 1);
        istream.read(reinterpret_cast<char*>(&sub_fields_len), 2);
        if (sid1 == 'A' && sid2 == 'P' && sub_fields_len == 2) {
            return GzipType::GZIP;
        }
        if (sid1 == 'B' && sid2 == 'C' && sub_fields_len == 2) {
            return GzipType::B_GZIP;
        }
        if (sid1 == 'N' && sid2 == 'A' && sub_fields_len == 2) {
            uint32_t bsize{0};
            istream.read(reinterpret_cast<char*>(&bsize), 4);
            istream.ignore(4);
            uint32_t isize;
            istream.read(reinterpret_cast<char*>(&isize), 4);
            return GzipType::NANO_B_GZIP;
        }
    }
    istream.close();
    return GzipType::GZIP;
}

std::vector<uint8_t> NanoBgzip::compress(std::vector<uint8_t>& input_data)
{
    /* create header, not text, not hcrc, not name, not comment, but need extra field */
    gz_header header;
    header.extra_len = 8; // the xlen is 8 bytes
    header.text = 0;
    header.hcrc = 0;
    header.name = nullptr;
    header.comment = nullptr;
    Bytef extra[8];
    /* NA - the first bytes in xlen means this is NANO_B_GZIP */
    extra[0] = 'N'; // sid1 = 'N'
    extra[1] = 'A'; // sid2 = 'A'
    uint16_t slen{2}; // uint16_t, subfield length = 2;
    extra[2] = static_cast<uint8_t>(slen & 0xff); // subfield length low byte
    extra[3] = static_cast<uint8_t>((slen & 0xff00) >> 8); // subfield length high byte
    /* uint32_t used to save the total block size, because the length of nanopore reads may be greater than multi Mb */
    uint32_t block_size{0};
    extra[4] = 0;
    extra[5] = 0;
    extra[6] = 0;
    extra[7] = 0;
    header.extra = extra;

    input_data.shrink_to_fit();
    //bigger output buffer used to store compressed data, 28 = header + sizeof(xlen) + xlen + crc32 + isize
    std::vector<uint8_t> output_data(input_data.size() * 1.2 + 28);

    z_stream strm;
    strm.zalloc = nullptr;
    strm.zfree = nullptr;
    strm.opaque = nullptr;
    strm.next_in = input_data.data();
    strm.avail_in = input_data.size();
    strm.next_out = output_data.data();
    strm.avail_out = output_data.size();
    int init_ret = deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY);
    if (init_ret != Z_OK) {
        cerr << "deflate init ERROR" << endl;
        deflateEnd(&strm);
        exit(1);
    }
    int set_header_ret = deflateSetHeader(&strm, &header);
    if (set_header_ret != 0) {
        cerr << "deflate set header failed" << endl;
        deflateEnd(&strm);
        exit(1);
    }
    // one time get N records of fastq, and call deflate once to compress all of them
    int deflate_ret = deflate(&strm,Z_FINISH);
    if (deflate_ret != Z_STREAM_END) {
        cerr << "Couldn't compress all input in one call" << endl;
        deflateEnd(&strm);
        exit(1);
    }
    deflateEnd(&strm);
    output_data.resize(strm.total_out);
    output_data.shrink_to_fit();
    block_size = output_data.size();
    uint32_t input_size = strm.total_in;
    // change the extra fields contents in header to save the total block size
    output_data[16] = static_cast<uint8_t>(block_size & 0XFF);
    output_data[17] = static_cast<uint8_t>((block_size & 0XFF00) >> 8);
    output_data[18] = static_cast<uint8_t>((block_size & 0XFF0000) >> 16);
    output_data[19] = static_cast<uint8_t>((block_size & 0XFF000000) >> 24);
    // ignore CRC32
    // save the input size in isize (the last bytes of this block)
    output_data[block_size - 4] = static_cast<uint8_t>(input_size & 0XFF);
    output_data[block_size - 3] = static_cast<uint8_t>((input_size & 0XFF00) >> 8);
    output_data[block_size - 2] = static_cast<uint8_t>((input_size & 0XFF0000) >> 16);
    output_data[block_size - 1] = static_cast<uint8_t>((input_size & 0XFF000000) >> 24);
    return output_data;
}
