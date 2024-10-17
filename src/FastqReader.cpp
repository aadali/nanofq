#include <filesystem>
#include <stdio.h>

#include <fmt/core.h>
//#include <cereal/types/unordered_map.hpp>
//#include <cereal/types/tuple.hpp>
//#include <cereal/types/string.hpp>
//#include <cereal/archives/binary.hpp>

#include "FastqReader.h"
#include "utility.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)
#define FASTQ_BUFFER_SIZE (1<<23)

using std::cout;
using std::endl;
using std::filesystem::last_write_time;
std::mutex FastqReader::ms_mtx{};
std::condition_variable FastqReader::ms_cond{};

FastqReader::FastqReader(std::string_view input_file, unsigned chunk)
        : m_input_file(input_file), m_chunk(chunk) {
    std::vector<std::shared_ptr<Read>> tmp;
    tmp.reserve(m_chunk);
    m_reads = std::make_shared<std::vector<std::shared_ptr<Read>>>(tmp);
    m_buffer = new char[FASTQ_BUFFER_SIZE];
    m_infile_gz = gzopen(input_file.data(), "rb");
    if (!m_infile_gz){
        throw std::runtime_error(fmt::format("Failed when opened file: {}", input_file));
    }
}

FastqReader::~FastqReader() {
    if (m_infile_gz) gzclose(m_infile_gz);
    if (m_buffer) {
        delete[] m_buffer;
        m_buffer = nullptr;
    }
}

int FastqReader::read_chunk_fastq() {
    std::string id, desc, sequence, quality;
    kseq_t *seq = kseq_init(m_infile_gz);
    int l;
    while (true) {
        l = kseq_read(seq);
        if (l == -1) break; // end of file
        if (l == -2) {
            throw std::runtime_error(fmt::format("Error: bad FASTQ format for read {}", seq->name.s));
        }
        if (l == -3) {
            throw std::runtime_error(fmt::format("Error reading {}", m_input_file));
        }
        bool fastq_format{seq->qual.l > 0 && seq->seq.l > 0 && seq->seq.l == seq->qual.l};
        if (!fastq_format) {
            throw std::runtime_error(fmt::format("\n\nError: could not parse input read \nproblem occurred at read {}", seq->name.s));
        }
        std::unique_lock<std::mutex> lock{ms_mtx};
        m_reads->emplace_back(std::make_shared<Read>(seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s));
        if (m_reads->size() == m_chunk) {
            std::cout << "first finished" << std::endl;
            ms_cond.wait(lock, [this]() { return m_reads->empty(); });
        }
    }
    kseq_destroy(seq);
    m_finish = true;
    return 0;
}

std::optional<shared_vec_reads> FastqReader::get_reads() {
    if (m_reads->size() == m_chunk || m_finish) {
        std::unique_lock<std::mutex> lock{ms_mtx};
        shared_vec_reads a{std::move(m_reads)};
        std::vector<std::shared_ptr<Read>> tmp;
        tmp.reserve(m_chunk);
        m_reads = std::make_shared<std::vector<std::shared_ptr<Read>>>(tmp);
        ms_cond.notify_all();
        return a;
    }
    return {};
}

std::unordered_set<std::string> FastqReader::get_searching_read_names(const std::string &input_reads) {
    std::unordered_set<std::string> read_names;
    const unsigned read_name_buf{256};
    if (exists(std::filesystem::path{input_reads.data()})) {
        // input_reads is a file
        // one read name per line in this file
        char read_name[read_name_buf];
        std::fstream infile{input_reads.data(), std::ios::in};
        if (!infile) {
            throw std::runtime_error(fmt::format("Failed when opened file: {}", input_reads.data()));
        }
        while (infile.getline(read_name, read_name_buf, '\n')) {
            bool empty_line{true};
            for (char &c: read_name) {
                if (isprint(c) && c != ' ') {
                    empty_line = false;
                    break;
                }
            }
            // empty line is supported
            // line starts with @ will be treated as comment
            if (empty_line || read_name[0] == '@') { continue; }
            if (read_name[strlen(read_name) - 1] == '\r') {
                read_name[strlen(read_name) - 1] = '\0';
            }
            read_names.emplace(read_name);
        }
    } else {
        std::vector<std::string_view> read_names_view{utility::split(input_reads, ",")};
        for (auto read_name: read_names_view) {
            read_names.emplace(read_name);
        }
    }
    for (auto &r: read_names) {
        cout << r << endl;
    }
    return read_names;
}

void FastqReader::find_reads(const std::string &input_reads, std::ostream &out, bool use_index) {
    std::ifstream infile_text {m_input_file.data(), std::ios::in};
    std::unordered_set<std::string> read_names{FastqReader::get_searching_read_names(input_reads)};
    if (use_index) {
        index();
        std::unordered_map<std::string, std::pair<size_t, size_t>> reads_index;
        std::fstream input_file_index{std::filesystem::path{m_input_file.data()}.concat(".idx").c_str(), std::ios::in};
        char index_line[1024];
        std::string read_name;
        size_t start{0}, stop{0};
        while (input_file_index.getline(index_line, 1024, '\n')) {
            if (strlen(index_line) < 4) break; // this is a simple judge
            std::vector<unsigned int> tab_pos;
            for (int i{0}; i < strlen(index_line); i++) {
                if (index_line[i] == '\t') { tab_pos.push_back(i); }
            }
            std::string this_line{index_line};
            start = std::stoull(this_line.substr(tab_pos[0] + 1, tab_pos[1] - tab_pos[0] - 1));
            stop = std::stoull(this_line.substr(tab_pos[1] + 1, this_line.size() - tab_pos[1] - 1));
//            start = atoll(this_line.substr(tab_pos[0] + 1, tab_pos[1] - tab_pos[0] - 1).data());
//            stop = atoll(this_line.substr(tab_pos[1] + 1, this_line.size() - tab_pos[1] - 1).data());
            reads_index.try_emplace(std::string(this_line.substr(0, tab_pos[0])), start, stop);
        }
        input_file_index.close();
        for (const std::string &id: read_names) {
            if (!reads_index.contains(id)) {
                std::cerr << fmt::format("There is no read named {} in fastq", id) << endl;
                continue;
            }
            auto [begin, end] = reads_index.at(id);
            infile_text.seekg(begin, std::ios::beg);
            for (size_t idx{0}; idx < end - begin; idx++) {
                out << static_cast<char>(infile_text.get());
            }
        }
        return;
    }
    infile_text.seekg(std::ios::beg);
    int line_number{1};
    bool find_read{false};
    std::string id;
   infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n');
    while (infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        if (!find_read && read_names.empty())break;
        if (m_buffer[strlen(m_buffer) - 1] == '\r') m_buffer[strlen(m_buffer) - 1] = '\0';
        switch (line_number) {
            case 1: {
                if (strlen(m_buffer) == 0) { break; }
                for (int idx{0}; idx < strlen(m_buffer); idx++) {
                    if (m_buffer[idx] == ' ') {
                        id.assign(m_buffer, 1, idx - 1);
                        break;
                    }
                }
                if (id.empty()) { id = m_buffer; };
                if (read_names.contains(id)) {
                    read_names.erase(id);
                    out << m_buffer << '\n';
                    find_read = true;
                }
                line_number++;
                break;
            }
            case 2:
            case 3: {
                if (find_read) out << m_buffer << '\n';
                line_number++;
                break;
            }

            case 4: {
                if (find_read) out << m_buffer << '\n';
                line_number = 1;
                find_read = false;
                break;
            }
            default: {
                break;
            }
        }
    }
    infile_text.close();
}

void FastqReader::index() {
    std::filesystem::path input_file{m_input_file};
    std::filesystem::path input_file_idx{input_file.concat(".idx")};
    if (exists(input_file_idx)) {
        if (last_write_time(input_file) > last_write_time(input_file_idx)) { // input_file is newer than index
            index(input_file_idx.c_str());
        }
    } else {
        index(input_file_idx.c_str());
    }
}

void FastqReader::index(std::string_view output_file_path) {
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::fstream output_index_stream{output_file_path.data(), std::ios::out};
    if (!output_index_stream) {
        throw std::runtime_error(fmt::format("Failed when opened file: {}", output_file_path.data()));
    }
//    cereal::BinaryOutputArchive bin_index{output_index_stream};
//    std::unordered_map<std::string, std::tuple<size_t, size_t>> reads_index{};
    infile_text.seekg(std::ios::beg);
    std::string id;
    int line_number{1};
    size_t start{0};
    size_t stop{0};
    while (infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        switch (line_number) {
            case 1: {
                for (int idx{0}; idx < strlen(m_buffer); idx++) {
                    if (m_buffer[idx] == ' ') {
                        id.assign(m_buffer, 1, idx - 1);
                        break;
                    }
                }
                if (id.empty()) id.assign(m_buffer, 1, strlen(m_buffer) - 1);
                stop += strlen(m_buffer) + 1;
                line_number++;
                break;
            }
            case 2:
            case 3: {
                stop += strlen(m_buffer) + 1;
                line_number++;
                break;
            }
            case 4: {
                stop += strlen(m_buffer) + 1;
                line_number = 1;
//                reads_index.try_emplace(id, start, stop);
                output_index_stream << fmt::format("{}\t{}\t{}\n", id, start, stop);
                start = stop;
                id.clear();
                break;
            }
            default: {
                break;
            }
        }
    }
//    bin_index(reads_index);
    output_index_stream.close();
    infile_text.close();
}


