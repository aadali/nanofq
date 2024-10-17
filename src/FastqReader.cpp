#include <filesystem>

#include <fmt/core.h>
//#include <cereal/types/unordered_map.hpp>
//#include <cereal/types/tuple.hpp>
//#include <cereal/types/string.hpp>
//#include <cereal/archives/binary.hpp>

#include "FastqReader.h"

using std::cout;
using std::endl;
using std::filesystem::last_write_time;
#define FASTQ_BUFFER_SIZE (1<<23)
std::mutex FastqReader::ms_mtx{};
std::condition_variable FastqReader::ms_cond{};

FastqReader::FastqReader(std::string_view input_file, unsigned chunk)
        : m_input_file(input_file), m_chunk(chunk) {
    std::vector<std::shared_ptr<Read>> tmp;
    tmp.reserve(m_chunk);
    m_reads = std::make_shared<std::vector<std::shared_ptr<Read>>>(tmp);
    m_buffer = new char[FASTQ_BUFFER_SIZE];
    m_infile = std::fstream{input_file.data(), std::ios::in};
    if (!m_infile) {
        throw std::runtime_error(fmt::format("Failed when opened file: {}", input_file));
    }
}

FastqReader::~FastqReader() {
    if (m_infile.is_open()) m_infile.close();
    if (m_buffer) {
        delete m_buffer;
        m_buffer = nullptr;
    }
}

void FastqReader::read_chunk_fastq() {
    std::string id, desc, sequence, quality;
    int line_number{1};
    while (m_infile.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        unsigned buf_len = strlen(m_buffer);
        if (m_buffer[buf_len - 1] == '\r') m_buffer[buf_len - 1] = '\0';
        switch (line_number) {
            case 1: {
                if (strlen(m_buffer) == 0) { break; }
                if (m_buffer[0] != '@') {
                    throw std::runtime_error(fmt::format("id must starts with @, which is {}", m_buffer));
                }
                for (int idx{0}; idx < buf_len; idx++) {
                    if (m_buffer[idx] == ' ') {
                        id.assign(m_buffer, 1, idx - 1);
                        desc.assign(m_buffer, idx + 1, buf_len - idx - 1);
                        break;
                    }
                }
                if (id.empty()) { id.assign(m_buffer, 1, strlen(m_buffer) - 1); }
                line_number++;
                break;
            }
            case 2: {
                if (buf_len == 0) {
                    throw std::runtime_error(fmt::format("empty sequence for {}", id));
                }
                sequence = m_buffer;
                line_number++;
                break;
            }
            case 3: {
                line_number++;
                break;
            }
            case 4: {
                line_number++;
                if (buf_len == 0) {
                    throw std::runtime_error(fmt::format("empty quality for {}", id));
                }
                quality = m_buffer;
                break;
            }
            default: {
                break;
            }
        }
        if (line_number == 5) {
            line_number = 1;
            if (sequence.size() != quality.size()) {
                throw std::runtime_error(fmt::format("different length between sequence({}) and quality({}) for {}",
                                                     sequence.size(), quality.size(), id));
            }
            std::unique_lock<std::mutex> lock{ms_mtx};
            m_reads->emplace_back(std::make_shared<Read>(id, desc, sequence, quality));
            if (m_reads->size() == m_chunk) {
                std::cout << "first finished" << std::endl;
                ms_cond.wait(lock, [this]() { return m_reads->empty(); });
            }
        }
    }
    m_finish = true;
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
        std::string input_reads2;
        if (!input_reads.ends_with(',')) {
            input_reads2 = input_reads + ',';
        }
        std::string read_name;
        auto start{input_reads2.begin()};
        for (auto it{std::begin(input_reads2)}; it != std::end(input_reads2); it++) {
            if (*it == ',') {
                read_name.assign(start, it);
                start = it + 1;
                if (!read_names.contains(read_name)) {
                    read_names.emplace(read_name);
                }
            }
        }
    }
    for (auto &r: read_names) {
        cout << r << endl;
    }
    return read_names;
}

void FastqReader::find_reads(const std::string &input_reads, std::ostream &out, bool use_index) {
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
            m_infile.seekg(begin, std::ios::beg);
            for (size_t idx{0}; idx < end - begin; idx++) {
                out << static_cast<char>(m_infile.get());
            }
        }
        return;
    }
    m_infile.seekg(std::ios::beg);
    int line_number{1};
    bool find_read{false};
    std::string id;
    while (m_infile.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
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
    std::fstream output_index_stream{output_file_path.data(), std::ios::out};
    if (!output_index_stream) {
        throw std::runtime_error(fmt::format("Failed when opened file: {}", output_file_path.data()));
    }
//    cereal::BinaryOutputArchive bin_index{output_index_stream};
//    std::unordered_map<std::string, std::tuple<size_t, size_t>> reads_index{};
    m_infile.seekg(std::ios::beg);
    std::string id;
    int line_number{1};
    size_t start{0};
    size_t stop{0};
    while (m_infile.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
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
}

