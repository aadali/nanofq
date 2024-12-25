#include "FastqReader.h"

using std::cout;
using std::endl;
using std::filesystem::last_write_time;

FastqReader::FastqReader(const std::string& input_file, int chunk)
    : m_input_file(input_file),
      m_input_file_index(std::filesystem::path{input_file}.concat(".index").c_str()),
      m_chunk_size(chunk)
{
    m_buffer = new char[FASTQ_BUFFER_SIZE];
    m_infile_gz = gzopen(input_file.data(), "rb");
    if (!m_infile_gz) {
        std::cerr << REDS + fmt::format("Failed opening file: {}", input_file) + COLOR_END << std::endl;;
        exit(1);
    }
    m_seq = kseq_init(m_infile_gz);
}

FastqReader::~FastqReader()
{
    if (m_infile_gz) gzclose(m_infile_gz);
    if (m_buffer) {
        delete[] m_buffer;
        m_buffer = nullptr;
    }
    if (m_seq) {
        kseq_destroy(m_seq);
    }
}

Read FastqReader::read_one_fastq()
{
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
    return read;
}


std::shared_ptr<std::vector<Read>> FastqReader::read_chunk_fastq()
{
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


std::unordered_set<std::string> FastqReader::get_searching_read_names(const std::string& input_reads)
{
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
    for (auto& r : read_names) {
        cout << r << endl;
    }
    return read_names;
}

void FastqReader::find_reads(const std::string& input_reads, std::ostream& out, bool use_index, unsigned key_len)
{
    std::ifstream infile_text{m_input_file, std::ios::in};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(key_len, false);
        auto reads_index{read_index()};
        size_t used_key_len{reads_index["used_key_len"][0]};

        char read_name[512];
        for (const std::string& id : read_names) {
            std::string key{id.size() <= used_key_len ? id : id.substr(0, used_key_len)};
            std::vector<size_t> idxes = reads_index[key];
            if (idxes.empty()) {
                std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) + COLOR_END <<
                    endl;
                continue;
            }
            size_t start_idx{0};
            size_t stop_idx{1};
            size_t start_offset{idxes[0]};
            size_t stop_offset{idxes[1]};
            bool find_read_name{false};
            while (stop_idx + 1 <= reads_index.size()) {
                if (find_read_name) break;
                infile_text.seekg(start_offset, std::ios::beg);
                infile_text.get(read_name, 512);
                if (std::string_view{read_name}.substr(1, strlen(read_name) - 1) == id) {
                    find_read_name = true;
                    out << read_name;
                    for (size_t idx{strlen(read_name)}; idx < stop_offset - start_offset+1; idx++) {
                        out << static_cast<char>(infile_text.get());
                    }
                }
                ++start_idx;
                ++stop_idx;
                if (!find_read_name) {
                    std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) + COLOR_END <<
                        endl;
                }
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
        find_reads(input_reads, out, true, key_len);
        return;
    }
    search_read_one_by_one(read_names, out);
}

void FastqReader::find_reads_in_gz(const std::string& input_reads, std::ostream& out, bool use_index,
                                   unsigned key_len)
{
    std::ifstream infile_text{m_input_file, std::ios::binary};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(key_len, false);
        auto reads_index{read_gz_index()};
        auto block_edges{reads_index.first};
        auto reads_index_in_block{reads_index.second};
        auto used_key_len = std::get<0>(reads_index.second["used_key_len"][0]);
        for (const std::string& id : read_names) {
            std::string key{id.size() <= used_key_len ? id : id.substr(0, used_key_len)};
            try{
                auto indexes= reads_index_in_block.at(key);
                bool find_read_name{false};
                for (auto&[block_index, read_start_index, read_end_index] : indexes){
                    if (find_read_name) break;
                    auto uncompressed_data {nanobgzip::get_uncompressed_from_block(infile_text, block_edges[block_index], read_end_index+1)};
                    std::string this_read_name{reinterpret_cast<char*>(uncompressed_data.data() + 1), id.size()};
                    if (this_read_name == id){
                        find_read_name = true;
                        for (size_t idx{read_start_index}; idx<read_end_index+1; ++idx){
                            cout << static_cast<char>(uncompressed_data[idx]);
                        }
                    }
                }
            } catch (const std::out_of_range& e){
                std::cerr << WARNS + fmt::format("There is no read named {} in this fastq file", id) +COLOR_END <<endl;
            }
        }
        return;
    }
    std::filesystem::path index_file{m_input_file_index};
    if (exists(index_file) && last_write_time(index_file) > last_write_time(std::filesystem::path{m_input_file})){
        find_reads(input_reads, out, true, key_len);
        return;
    }
    search_read_one_by_one(read_names, out);
}


void FastqReader::search_read_one_by_one(std::unordered_set<std::string>& read_names, std::ostream& out)
{
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

void FastqReader::index(unsigned key_len, bool force_index)
{
    // TODO when key_len, use the full length read name to index
    std::filesystem::path input_file{m_input_file};
    std::filesystem::path input_file_idx{m_input_file_index};
    if (force_index){
        m_input_file.ends_with(".gz") ? index_fastq_gz(key_len) : index_fastq(key_len);
        return;
    }
    if (exists(input_file_idx)) {
        if (last_write_time(input_file) > last_write_time(input_file_idx)) {
            // input_file is newer than index
            m_input_file.ends_with(".gz") ? index_fastq_gz(key_len) : index_fastq(key_len);
        }
    } else {
        m_input_file.ends_with(".gz") ? index_fastq_gz(key_len) : index_fastq(key_len);
    }
}

void FastqReader::find(const std::string& input_reads, std::ostream& out, bool use_index, unsigned key_len) {
    if (m_input_file.ends_with(".gz")){
        this->find_reads_in_gz(input_reads, out, use_index, key_len);
    } else {
        this->find_reads(input_reads, out, use_index, key_len);
    }
}

void FastqReader::index_fastq(unsigned key_len)
{
    key_len = key_len == 0 ? 999 : key_len;
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::ofstream output_index_stream{m_input_file_index, std::ios::out};
    if (!output_index_stream) {
        std::cerr << REDS + fmt::format("Failed when opened file: {}", m_input_file_index) + COLOR_END <<
            std::endl;
        exit(1);
    }
    std::unordered_map<std::string, std::vector<size_t>> reads_index{};
    infile_text.seekg(std::ios::beg);
    std::string id;
    int line_number{1};
    // 0-based coordinate, close interval
    size_t start{0}; // the start index in file
    size_t length{0}; // the length of record, include newline of last line
    while (infile_text.getline(m_buffer, FASTQ_BUFFER_SIZE, '\n')) {
        switch (line_number) {
        case 1: {
            std::string_view read_name_prefix{myutility::get_read_name_prefix(m_buffer, key_len)};
            id = read_name_prefix;
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

void FastqReader::index_fastq_gz(unsigned key_len)
{
    key_len = key_len == 0 ? 999 : key_len;
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
        std::cerr << "\nzcat input.fastq.gz | nanofq compress - output.fastq.gz\n\n";
        std::cerr << REDS <<
            "The above command will turn comman input.fastq.gz into NanoBgzip output.fastq.gz and create index file simultaneously\n";
        std::cerr <<
            "Or use bgzip to compress the text file and samtools fqidx to build index. Refer the following command:\n"
            << COLOR_END << std::endl;
        std::cerr << "zcat input.fastq.gz | bgzip -c > output.fastq.gz && samtools fqidx output.fastq.gz" << std::endl;
        exit(1);
    }
    nanobgzip::build_index(std::string{m_input_file}, m_input_file_index, key_len);
}

std::unordered_map<std::string, std::vector<size_t>> FastqReader::read_index() const
{
    auto index_path{std::filesystem::path{m_input_file_index}};
    if (!exists(index_path)) {
        std::cerr << "No such file: " << index_path << endl;
        exit(1);
    }
    if (m_input_file_index.ends_with(".gz.index")) {
        std::cerr << "read_index is not suitable for gz file" << endl;
        exit(1);
    }
    std::unordered_map<std::string, std::vector<size_t>> reads_index;
    std::ifstream infile{m_input_file_index, std::ios::in};
    char index_line[512];
    infile.getline(index_line, 512, '\n');
    unsigned used_key_len;
    std::from_chars(index_line + 1, index_line + strlen(index_line), used_key_len);
    reads_index["used_key_len"] = std::vector<size_t>{used_key_len}; // store the key length
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
        reads_index[std::string{read_name}].push_back(start);
        reads_index[std::string{read_name}].push_back(end);
    }
    infile.close();
    return reads_index;
}

nanobgzip_reads_index FastqReader::read_gz_index() const
{
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
     * .second => std::unordered_map<std::string, std::vector<std::tuple<unsigned, size_t, size_t>>>
                keys: read name or  the first used_key_len chars of reads depended on the key_len parameter when built index
                values: vector of std::tuple<the_block_index_that_read_belongs_to_in_first, read_start_index_in_this_block, read_end_index_in_the_block>
                In most cases, the size of this vector should be 1.
                But sometimes multi read names shared one prefix (key_len), so the vector size will be larger than 1
                Note this situation
    */
    nanobgzip_reads_index reads_index;
    char index_line[512];
    infile.getline(index_line, 512, '\n');
    unsigned used_key_len;
    std::from_chars(index_line + 1, index_line + strlen(index_line), used_key_len);
    reads_index.first.emplace_back(0, 0); // store the used_key_len
    reads_index.second["used_key_len"].emplace_back(used_key_len, used_key_len, used_key_len);
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
        reads_index.second[std::string{read_name}].emplace_back(reads_index.first.size() - 1, read_start, read_end);
    }
    infile.close();
    return reads_index;
}
