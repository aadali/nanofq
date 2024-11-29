#include "NanoBgzip.h"

#include <vector>


void NanoBgzip::nano_compress(const std::string& infile, const std::string& outfile, int reads_number)
{
    std::istream* is = nullptr;
    bool using_infile = false;
    if (infile == "-") {
        is = &std::cin;
    } else {
        if (infile.ends_with(".gz")) {
            cerr << "Couldn't compress *.gz file again" << endl;
            exit(1);
        }
        is = new std::ifstream{infile};
        using_infile = true;
    }


    if (!outfile.ends_with(".gz")) {
        cerr << "outfile for compress results should be ends with .gz" << endl;
        exit(1);
    }
    std::ofstream ostream{outfile, std::ios::binary};
    std::filesystem::path index_file{std::filesystem::path{outfile}.concat(".idx")};
    if (!ostream){
        cerr << "Opened " << outfile << " failed" << endl;
        exit(1);
    }
    std::ofstream index_stream {index_file};
    if (!index_stream) {
        cerr << "Opened " << index_file << " failed for make index" << endl;
        exit(1);
    }

    std::vector<uint8_t> in_buffer;
    char ch;
    int newline_count{0};
    size_t written_bytes{0};
    while (is->get(ch)) {
        in_buffer.push_back(ch);
        if (ch == '\n') {
            ++newline_count;
            if (newline_count == 4 * reads_number) {
            // TODO consider ThreadPool
                std::vector<uint8_t> output_buffer{nano_block_compress(in_buffer, index_stream, written_bytes)};
                written_bytes += output_buffer.size();
                ostream.write(reinterpret_cast<char*>(output_buffer.data()), output_buffer.size());
                newline_count = 0;
                std::vector<uint8_t>{}.swap(in_buffer);
            }
        }
    }
    if (is->eof()) {
        if (!in_buffer.empty()) {
            std::vector<uint8_t> output_buffer{nano_block_compress(in_buffer, index_stream, written_bytes)};
            ostream.write(reinterpret_cast<char*>(output_buffer.data()), output_buffer.size());
        }
    } else {
        cerr << "Read eof of file " << infile << "Failed" << endl;
        exit(1);
    }
    if (using_infile) {
        delete is;
        is = nullptr;
    }
    ostream.close();
    index_stream.close();
}

GzipType NanoBgzip::check_compress_type(const std::string& infile)
{
    // https://www.ietf.org/rfc/rfc1952.txt
    std::ifstream istream{infile, std::ios::binary};
    uint8_t id1{0};
    uint8_t id2{0};
    istream.seekg(std::ios::beg);
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
    return GzipType::GZIP;
}

std::vector<uint8_t> NanoBgzip::nano_block_compress(std::vector<uint8_t>& input_data, std::ostream& output_index_stream,
                                                    size_t written_bytes)
{
    std::vector<std::string> index_in_block{get_index_in_block(input_data)};
    /* create header, not text, no hcrc, no name, no comment, but need extra field */
    gz_header header;
    header.extra_len = 8; // set the xlen is 8 bytes => 1 byte sid1, 1 byte sid2 , 2 bytes slen, 4 bytes block size}
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
    // one time get N records of fastq, and call deflate once to nano_block_compress all of them
    int deflate_ret = deflate(&strm,Z_FINISH);
    if (deflate_ret != Z_STREAM_END) {
        cerr << "Couldn't nano_block_compress all input in one call" << endl;
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
    output_index_stream << "#" << written_bytes << '\t' << written_bytes + block_size -1 << '\n';
    for (auto& read_position : index_in_block) {
        output_index_stream << read_position;
    }
    return output_data;
}

std::vector<std::string> NanoBgzip::get_index_in_block(const std::vector<uint8_t>& input_data)
{
    std::vector<std::string> reads_position;
    int line_number{0};

    // std::string read_name{};
    std::stringstream this_read_position{};
    bool find_read_name{false};
    int record_start_idx{0};
    int record_end_idx{0};
    for (int i{0}; i < input_data.size(); ++i) {
        if (line_number == 0 && !find_read_name) {
            if (input_data[i] != ' ' && input_data[i] != '\n') {
                this_read_position << input_data[i];
            } else {
                find_read_name = true;
            }
        }
        if (input_data[i] == '\n' || i == input_data.size() - 1) {
            ++line_number;
            if (line_number == 4) {
                record_end_idx = i;
                this_read_position << '\t' << record_start_idx << '\t' << record_end_idx << '\n';
                reads_position.emplace_back(this_read_position.str().substr(1, this_read_position.str().size() - 1));
                record_start_idx = record_end_idx + 1;
                line_number = 0;
                find_read_name = false;
                this_read_position.str("");
                this_read_position.clear();
            }
        }
    }
    return reads_position;
}
