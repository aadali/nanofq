#include <filesystem>
#include <cstdio>
#include <charconv>

#include <fmt/core.h>

#include "FastqReader.h"
#include "myUtility.h"

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
    if (!m_infile_gz) {
        std::cerr << REDS + fmt::format("Failed opening file: {}", input_file) + COLOR_END << std::endl;;
        exit(1);
    }
    m_seq = kseq_init(m_infile_gz);
}

FastqReader::~FastqReader() {
    if (m_infile_gz) gzclose(m_infile_gz);
    if (m_buffer) {
        delete[] m_buffer;
        m_buffer = nullptr;
    }
    // if (m_seq) { kseq_destroy(m_seq); }
}

Read FastqReader::read_one_fastq() {
    int l;
    l = kseq_read(m_seq);
    if (l == -1) {
        // end of file
        std::string id = finished_read_name, qual = id, seq = id, desc = id;
        Read read{id, desc, seq, qual};
        return read;
    }
    if (l == -2) {
        std::cerr << REDS + fmt::format("Error: bad FASTQ format for read {}", m_seq->name.s) + COLOR_END << std::endl;;
        exit(1);
    }
    if (l == -3) {
        std::cerr << REDS + fmt::format("Error reading {}", m_input_file) + COLOR_END << std::endl;;
        exit(1);
    }
    bool fastq_format{m_seq->qual.l > 0 && m_seq->seq.l > 0 && m_seq->seq.l == m_seq->qual.l};
    if (!fastq_format) {
        std::cerr << REDS + fmt::format("\n\nError: could not parse input read \nproblem occurred at read {}",
                                        m_seq->name.s) + COLOR_END << std::endl;;
        exit(1);
    }
    Read read{m_seq->name.s, m_seq->comment.s, m_seq->seq.s, m_seq->qual.s};
    // m_reads->emplace_back(std::make_shared<Read>(std::move(read)));
    return read;
}


void FastqReader::read_chunk_fastq() {
    // std::string id, desc, sequence, quality;
    // int l;
    while (true) {
        Read read{read_one_fastq()};
        if (read.get_id() == finished_read_name) break;
        std::unique_lock<std::mutex> lock{ms_mtx};
        m_reads->emplace_back(std::make_shared<Read>(std::move(read)));
        if (m_reads->size() == m_chunk) {
            // std::cout << "first finished" << std::endl;
            ms_cond.wait(lock, [this](){ return m_reads->empty(); });
        }
    }
    kseq_destroy(m_seq);
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

std::unordered_set<std::string> FastqReader::get_searching_read_names(const std::string& input_reads) {
    std::unordered_set<std::string> read_names;
    const unsigned read_name_buf{256};
    if (exists(std::filesystem::path{input_reads.data()})) {
        // input_reads is a file
        // one read name per line in this file
        char read_name[read_name_buf];
        std::fstream infile{input_reads.data(), std::ios::in};
        if (!infile) {
            std::cerr << REDS + fmt::format("Failed when opened file: {}", input_reads.data()) + COLOR_END << std::endl;
            exit(1);
        }
        while (infile.getline(read_name, read_name_buf, '\n')) {
            bool empty_line{true};
            for (char& c : read_name) {
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
        std::vector<std::string_view> read_names_view{myUtility::split(input_reads, ",")};
        for (auto read_name : read_names_view) {
            read_names.emplace(read_name);
        }
    }
    for (auto& r : read_names) {
        cout << r << endl;
    }
    return read_names;
}

void FastqReader::find_reads(const std::string& input_reads, std::ostream& out, bool use_index, unsigned key_length) {
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    std::filesystem::path index_file_path_prefix{m_input_file.data()};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(key_length);
        std::unordered_multimap<std::string, std::pair<size_t, size_t>> reads_index;
        std::fstream input_file_index{index_file_path_prefix.concat(".idx").c_str(), std::ios::in};
        char index_line[1024];
        input_file_index.getline(index_line, 1024, '\n');
        unsigned used_key_length;
        std::from_chars(index_line + 1, index_line + strlen(index_line), used_key_length);
        size_t start{0}, stop{0};
        while (input_file_index.getline(index_line, 1024, '\n')) {
            if (strlen(index_line) < 4) break; // this is a simple judge
            std::vector<unsigned int> tab_pos;
            auto this_read_index{myUtility::split(index_line, "\t")};
            std::from_chars(this_read_index[1].data(), this_read_index[1].data() + this_read_index[1].size(), start);
            std::from_chars(this_read_index[2].data(), this_read_index[2].data() + this_read_index[2].size(), stop);
            reads_index.emplace(std::string{this_read_index[0]}, std::make_pair(start, stop));
        }
        input_file_index.close();
        char read_name[4096];
        for (const std::string& id : read_names) {
            string key{id.size() <= used_key_length ? id : id.substr(0, used_key_length)};
            auto [begin_ele, end_ele] = reads_index.equal_range(key);
            bool find_read_name{false};
            for (auto ele{begin_ele}; ele != end_ele; ele++) {
                if (find_read_name) break;
                infile_text.seekg(ele->second.first, std::ios::beg);
                infile_text.get(read_name, 4096);
                if (std::string_view{read_name}.substr(1, strlen(read_name) - 1) == id) {
                    find_read_name = true;
                    out << read_name;
                    for (size_t idx{strlen(read_name)}; idx < ele->second.second - ele->second.first; idx++) {
                        out << static_cast<char>(infile_text.get());
                    }
                }
            }
            if (!find_read_name) {
                std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) + COLOR_END <<
                    endl;
            }
        }
        return;
    }
    if (exists(index_file_path_prefix.concat(".idx")) &&
        last_write_time(index_file_path_prefix.concat(".idx")) > last_write_time(std::filesystem::path{m_input_file})) {
        /*
        * if user didn't set --use_index, check whether the index file exists firstly, if true and index file is newer than input
        * file, so just use it, else do the following iteration searching
        * */
        find_reads(input_reads, out, true, key_length);
        return;
    }
    while (true) {
        // iteration searching
        int l;
        l = kseq_read(m_seq);
        if (l == -1 || read_names.empty()) {
            // end of file or find all reads
            if (!read_names.empty()) {
                for (std::string name : read_names) {
                    std::cerr << WARNS + fmt::format("There is no read named {} in this fastq", name) << COLOR_END <<
                        std::endl;
                }
            }
            kseq_destroy(m_seq);
            return;
        }
        if (l == -2) {
            std::cerr << REDS + fmt::format("Error: bad FASTQ format for read {}", m_seq->name.s) + COLOR_END <<
                std::endl;;
            exit(1);
        }
        if (l == -3) {
            std::cerr << REDS + fmt::format("Error reading {}", m_input_file) + COLOR_END << std::endl;;
            exit(1);
        }
        bool fastq_format{m_seq->qual.l > 0 && m_seq->seq.l > 0 && m_seq->seq.l == m_seq->qual.l};
        if (!fastq_format) {
            std::cerr << REDS + fmt::format("\n\nError: could not parse input read \nproblem occurred at read {}",
                                            m_seq->name.s) + COLOR_END << std::endl;;
            exit(1);
        }
        if (read_names.contains(m_seq->name.s)) {
            out << fmt::format("@{} {}\n{}\n+\n{}\n",
                               m_seq->name.s,
                               m_seq->comment.s ? m_seq->comment.s : "",
                               m_seq->seq.s,
                               m_seq->qual.s);
            read_names.erase(m_seq->name.s);
        }
    }
}

void FastqReader::index(unsigned key_len) {
    std::filesystem::path input_file{m_input_file};
    std::filesystem::path input_file_idx{input_file.concat(".idx")};
    if (exists(input_file_idx)) {
        if (last_write_time(input_file) > last_write_time(input_file_idx)) {
            // input_file is newer than index
            index_fastq(input_file_idx.c_str(), key_len);
        }
    } else {
        index_fastq(input_file_idx.c_str(), key_len);
    }
}

void FastqReader::index_fastq(std::string_view output_file_path, unsigned key_len) {
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::fstream output_index_stream{output_file_path.data(), std::ios::out};
    if (!output_index_stream) {
        std::cerr << REDS + fmt::format("Failed when opened file: {}", output_file_path.data()) + COLOR_END <<
            std::endl;
        exit(1);
    }
    output_index_stream << '#' << key_len << '\n';
    //    cereal::BinaryOutputArchive bin_index{output_index_stream};
    std::unordered_map<std::string, std::vector<size_t>> reads_index{};
    infile_text.seekg(std::ios::beg);
    std::string id;
    int line_number{1};
    size_t start{0};
    size_t stop{0};
    while (infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        switch (line_number) {
        case 1: {
            string_view read_name_prefix{myUtility::get_read_name_prefix(m_buffer, key_len)};
            id = read_name_prefix;
            //                if (!reads_index.contains(id)){
            //                    reads_index.insert(std::make_pair(id, std::vector<size_t>{}));
            //                }
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
            output_index_stream << id << '\t' << start << '\t' << stop << '\n';
            //                reads_index.try_emplace(id, start, stop);
            //                reads_index.at(id).push_back(start);
            //                reads_index.at(id).push_back(stop);
            start = stop;
            //                id.clear();
            break;
        }
        default: {
            break;
        }
        }
    }
    output_index_stream.close();
    infile_text.close();
}
