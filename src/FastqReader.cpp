#include "FastqReader.h"

using std::cout;
using std::endl;
using std::filesystem::last_write_time;

FastqReader::FastqReader(const std::string& input_file, int chunk)
    : m_input_file(input_file),
      m_input_file_index(std::filesystem::path{input_file}.concat(".index").c_str()),
      m_chunk_size(chunk) {
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
    if (m_seq) {
        kseq_destroy(m_seq);
    }
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
        std::cerr << REDS + fmt::format("Error: bad FASTQ format for read {}", m_seq->name.s) + COLOR_END << std::endl;
        exit(1);
    }
    if (l == -3) {
        std::cerr << REDS + fmt::format("Error reading {}", m_input_file) + COLOR_END << std::endl;
        exit(1);
    }
    bool fastq_format{ m_seq->seq.l == m_seq->qual.l && m_seq->seq.l > 0};
    if (!fastq_format) {
        std::cerr << REDS + fmt::format("\n\nError: could not parse input read \nproblem occurred at read {}",
                                        m_seq->name.s) + COLOR_END << std::endl;;
        exit(1);
    }
    Read read{m_seq->name.s, m_seq->comment.s, m_seq->seq.s, m_seq->qual.s};
    return read;
}



std::optional<std::vector<std::filesystem::path>> FastqReader::get_fastqs() const {
    std::vector<std::filesystem::path> fastqs_paths;
    auto input_path{std::filesystem::path{m_input_file}};
    if (!std::filesystem::is_directory(input_path)) return {};
    for (const auto&p: std::filesystem::directory_iterator(input_path)){
        if (std::filesystem::is_directory(input_path)){
            std::cerr << REDS << "Input is directory, file and directory found simultaneously, but only files allowed in directory" << COLOR_END << std::endl;
            exit(1);
        }
        fastqs_paths.push_back(p);
    }
    return fastqs_paths;
}


std::shared_ptr<std::vector<Read>> FastqReader::read_chunk_fastq() {
    std::shared_ptr<std::vector<Read>> reads_ptr = std::make_shared<std::vector<Read>>();
    reads_ptr->reserve(m_chunk_size);
    while (true) {
        reads_ptr->emplace_back(read_one_fastq());
        if (reads_ptr->size() == m_chunk_size || reads_ptr->back().get_id() == finished_read_name) {
            if (reads_ptr->back().get_id() == finished_read_name) {
                reads_ptr->pop_back();
                m_finish = true;
            }
            break;
        };
    }
    reads_ptr->shrink_to_fit();
    return reads_ptr;
}


std::unordered_set<std::string> FastqReader::get_searching_read_names(const std::string& input_reads) {
    std::unordered_set<std::string> read_names;
    constexpr unsigned read_name_buf{256};
    if (exists(std::filesystem::path{input_reads.data()})) {
        // input_reads is a file, one read name per line in this file
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
            // line starts with # will be treated as comment
            if (empty_line || read_name[0] == '#') { continue; }
            if (read_name[strlen(read_name) - 1] == '\r') {
                read_name[strlen(read_name) - 1] = '\0';
            }
            read_names.emplace(read_name);
        }
    } else {
        std::vector<std::string_view> read_names_view{myutility::split(input_reads, ",")};
        for (auto read_name : read_names_view) {
            read_names.emplace(read_name);
        }
    }
    return read_names;
}

void FastqReader::find_reads(const std::string& input_reads, std::ostream& out, bool use_index) {
    std::ifstream infile_text{m_input_file, std::ios::in};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(false);
        auto reads_index{read_index()};
        for (const std::string& id : read_names) {
            std::pair<size_t, size_t> idxes;
            if (!reads_index.contains(id)) {
                std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) + COLOR_END <<
                    endl;
                continue;
            }
            idxes = reads_index[id];
            size_t start_idx{idxes.first};
            size_t stop_idx{idxes.second};
            infile_text.seekg(start_idx, std::ios::beg);
            for (size_t idx{start_idx}; idx < stop_idx + 1; ++idx) {
                out << static_cast<char>(infile_text.get());
            }
        }
        return;
    }
    std::filesystem::path index_file{m_input_file_index};
    if (exists(index_file) && last_write_time(index_file) > last_write_time(std::filesystem::path{m_input_file})) {
        /*
        * if user didn't set --use_index, check whether the index file exists firstly, if true and index file is newer than input
        * file, so just use it, else do the following iteration searching
        * */
        find_reads(input_reads, out, true);
        return;
    }
    search_read_one_by_one(read_names, out);
}

void FastqReader::find_reads_in_gz(const std::string& input_reads, std::ostream& out, bool use_index) {
    std::ifstream infile_text{m_input_file, std::ios::binary};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(false);
        auto reads_index{read_gz_index()};
        auto block_edges{reads_index.first};
        auto reads_index_in_block{reads_index.second};
        for (const std::string& id : read_names) {
            if (!reads_index_in_block.contains(id)) {
                std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) + COLOR_END <<
                    std::endl;
                continue;
            }
            auto indexes = reads_index_in_block[id];
            bool find_read_name{false};
            auto& [block_index, read_start_index, read_end_index] = indexes;
            auto uncompressed_data{
                nanobgzip::get_uncompressed_from_block(infile_text, block_edges[block_index],
                                                       read_end_index + 1)
            };
            for (size_t idx{read_start_index}; idx < read_end_index + 1; ++idx) {
                out << static_cast<char>(uncompressed_data[idx]);
            }
        }
        return;
    }
    std::filesystem::path index_file{m_input_file_index};
    if (exists(index_file) && last_write_time(index_file) > last_write_time(std::filesystem::path{m_input_file})) {
        find_reads(input_reads, out, true);
        return;
    }
    search_read_one_by_one(read_names, out);
}


void FastqReader::search_read_one_by_one(std::unordered_set<std::string>& read_names, std::ostream& out) {
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
            // kseq_destroy(m_seq);
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

void FastqReader::index(bool force_index) {
    std::filesystem::path input_file{m_input_file};
    std::filesystem::path input_file_idx{m_input_file_index};
    if (force_index) {
        // when use subcommand index, forcibly make index file at any case
        m_input_file.ends_with(".gz") ? index_fastq_gz() : index_fastq();
        return;
    }
    if (exists(input_file_idx)) {
        // when use subcommand find, make index depended on user's set
        if (last_write_time(input_file) > last_write_time(input_file_idx)) {
            // input_file is newer than index
            m_input_file.ends_with(".gz") ? index_fastq_gz() : index_fastq();
        }
    } else {
        m_input_file.ends_with(".gz") ? index_fastq_gz() : index_fastq();
    }
}

void FastqReader::find(const std::string& input_reads, std::ostream& out, bool use_index) {
    if (m_input_file.ends_with(".gz")) {
        this->find_reads_in_gz(input_reads, out, use_index);
    } else {
        this->find_reads(input_reads, out, use_index);
    }
}

Read FastqReader::fastq_record_ok(int l, kseq_t* seq, const char* file) {
    if (l == -1){
        std::string id = finished_read_name, qual = id, seq1 = id, desc =id;
        Read read{id, desc, seq1, qual};
        return read;
    }
    if (l == -2){
        std::cerr << REDS + fmt::format("Error: bad FASTQ format for read {}", seq->name.s) + COLOR_END << std::endl;
        exit(1);
    }
    if (l == -3) {
        std::cerr <<  REDS + fmt::format("Error reading {}", file) + COLOR_END << std::endl;
        exit(1);
    }
    bool fastq_format {seq->seq.l == seq->qual.l && seq->seq.l > 0};
    if (!fastq_format){
        std::cerr << REDS + fmt::format("\n\nError: could not parse input read \nproblem occurred at read {}", seq->name.s) + COLOR_END << std::endl;
        exit(1);
    }
    Read read{seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s};
    return read;
}

void FastqReader::index_fastq() {
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::ofstream output_index_stream{m_input_file_index, std::ios::out};
    if (!output_index_stream) {
        std::cerr << REDS + fmt::format("Failed when opened file: {}", m_input_file_index) + COLOR_END <<
            std::endl;
        exit(1);
    }
    // std::unordered_map<std::string, std::pair<size_t, size_t>> reads_index{};
    infile_text.seekg(std::ios::beg);
    std::stringstream read_name;
    std::string id;
    int line_number{1};
    // 0-based coordinate, close interval
    size_t start{0}; // the start index in file
    size_t length{0}; // the length of record, include newline of last line
    while (infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        switch (line_number) {
        case 1: {
            // std::string_view read_name_prefix{myutility::get_read_name_prefix(m_buffer, key_len)};
            for (size_t i{0}; i < strlen(m_buffer); i++) {
                char this_char{m_buffer[i]};
                if (this_char == ' ' || this_char == '\n') {
                    break;
                }
                read_name << this_char;
            }
            id = read_name.str().substr(1);
            length += strlen(m_buffer) + 1;
            line_number++;
            break;
        }
        case 2:
        case 3: {
            length += strlen(m_buffer) + 1;
            line_number++;
            break;
        }
        case 4: {
            length += strlen(m_buffer) + 1;
            line_number = 1;
            output_index_stream << id << '\t' << start << '\t' << start + length - 1 << '\n';
            /* readName\tstartIndex\tendIndex
             * 0-based, close interval
            */
            start += length;
            length = 0;
            read_name.str("");
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

void FastqReader::index_fastq_gz() {
    if (!m_input_file.ends_with(".gz")) {
        std::cerr << REDS + fmt::format("Error: input file path {} must ends with .gz when use index_fastq_gz",
                                        m_input_file) + COLOR_END << std::endl;
        exit(1);
    }
    // check_compress_type(std::string{output_file_path});
    if (nanobgzip::check_compress_type(std::string{m_input_file}) != nanobgzip::GzipType::NANO_B_GZIP) {
        std::cerr << REDS <<
            "Error: couldn't index common gzip file. \nYou can refer to the following command to convert it into a NanoBgzip file and build an index at the same time"
            << COLOR_END << std::endl;
        auto idx = m_input_file.rfind(".gz");
        std::cerr << fmt::format("\nzcat {} | nanofq compress - {}nanobgzip.gz\n\n", m_input_file,
                                 m_input_file.substr(0, idx + 1));
        // std::cerr << "\nzcat input.fastq.gz | nanofq compress - output.fastq.gz\n\n";
        std::cerr << REDS <<
            fmt::format(
                "The above command will turn common {} into NanoBgzip file: {}.nanobgzip.gz and create index file simulaneusly\n",
                m_input_file, m_input_file.substr(0, idx + 1));
        exit(1);
    }
    nanobgzip::build_index(std::string{m_input_file}, m_input_file_index);
}

std::unordered_map<std::string, std::pair<size_t, size_t>> FastqReader::read_index() const {
    auto index_path{std::filesystem::path{m_input_file_index}};
    if (!exists(index_path)) {
        std::cerr << "No such index file: " << index_path << endl;
        exit(1);
    }
    if (m_input_file_index.ends_with(".gz.index")) {
        std::cerr << "You are using gz.index file and text fastq that is not suitable" << endl;
        exit(1);
    }
    std::unordered_map<std::string, std::pair<size_t, size_t>> reads_index;
    std::ifstream infile{m_input_file_index, std::ios::in};
    char index_line[512];
    size_t start{0}, end{0};
    while (infile.getline(index_line, 512, '\n')) {
        if (strlen(index_line) < 4) {
            std::cerr << "Too short for index line" << std::endl;
            exit(1);
        }
        auto this_read_index{myutility::split(index_line, "\t")};
        std::string_view read_name{this_read_index[0]};
        std::string_view start_sv{this_read_index[1]};
        std::string_view end_sv{this_read_index[2]};
        std::from_chars(start_sv.data(), start_sv.data() + start_sv.size(), start);
        std::from_chars(end_sv.data(), end_sv.data() + end_sv.size(), end);
        reads_index[std::string{read_name}].first = start;
        reads_index[std::string{read_name}].second = end;
    }
    infile.close();
    return reads_index;
}

nanobgzip_reads_index FastqReader::read_gz_index() const {
    if (auto index_path{std::filesystem::path{m_input_file_index}}; !exists(index_path)) {
        std::cerr << "No such file: " << index_path << endl;
        exit(1);
    }
    if (!m_input_file_index.ends_with(".gz.index")) {
        std::cerr << "read_gz_index need the index file ends_with .gz.index" << endl;
        exit(1);
    }
    std::ifstream infile{m_input_file_index, std::ios::in};
    /* type of nanobgzip_reads_index => std::pair
     * .first => std::vector<std::pair<size_t, size_t>>: the vector of std::paired<block_start_index, block_end_index>
     * .second => std::unordered_map<std::string,std::tuple<unsigned, size_t, size_t>>
                keys: read name
                values: std::tuple<the_block_index_that_read_belongs_to_in_first, read_start_index_in_this_block, read_end_index_in_the_block>
    */
    nanobgzip_reads_index reads_index;
    char index_line[512];
    while (infile.getline(index_line, 512, '\n')) {
        if (strlen(index_line) < 4) {
            std::cerr << "Too short for nanobgzip index line" << std::endl;
            exit(1);
        }
        auto this_read_index{myutility::split(index_line, "\t")};
        if (index_line[0] == '#') {
            size_t block_start{0}, block_end{0};
            std::string_view block_start_sv{this_read_index[0].substr(1, this_read_index[0].size() - 1)};
            std::string_view block_end_sv{this_read_index[1]};
            std::from_chars(block_start_sv.data(), block_start_sv.data() + block_start_sv.size(), block_start);
            std::from_chars(block_end_sv.data(), block_end_sv.data() + block_end_sv.size(), block_end);
            reads_index.first.emplace_back(block_start, block_end);
            continue;
        }
        std::string_view read_name{this_read_index[0]};
        std::string_view read_start_sv{this_read_index[1]};
        std::string_view read_end_sv{this_read_index[2]};
        size_t read_start{0}, read_end{0};
        std::from_chars(read_start_sv.data(), read_start_sv.data() + read_start_sv.size(), read_start);
        std::from_chars(read_end_sv.data(), read_end_sv.data() + read_end_sv.size(), read_end);
        reads_index.second.try_emplace(std::string{read_name}, reads_index.first.size() - 1, read_start, read_end);
    }
    infile.close();
    return reads_index;
}
