#include "FastqReader.h"

using std::cout;
using std::endl;
using std::filesystem::last_write_time;

FastqReader::FastqReader(std::string_view input_file, int chunk)
    : m_input_file(input_file), m_chunk_size(chunk)
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

void FastqReader::find_reads(const std::string& input_reads, std::ostream& out, bool use_index, unsigned key_length)
{
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    std::unordered_set<std::string> read_names{get_searching_read_names(input_reads)};
    std::filesystem::path index_file_path_prefix{m_input_file.data()};
    if (use_index) {
        /* when user need use index, firstly check whether the index file exists.
        If true and index is newer than input file, just use it, else make index and use it */
        index(key_length);
        std::unordered_map<std::string, std::vector<size_t>> reads_index;
        std::ifstream input_file_index{index_file_path_prefix.concat(".index").c_str(), std::ios::in};
        char index_line[512];
        input_file_index.getline(index_line, 512, '\n');
        unsigned used_key_length;
        std::from_chars(index_line + 1, index_line + strlen(index_line), used_key_length);
        size_t start{0}, stop{0};
        while (input_file_index.getline(index_line, 512, '\n')) {
            if (strlen(index_line) < 4) {
                std::cerr << "Too short for index line" << std::endl;
                exit(1);
            } // this is a simple judge
            std::vector<unsigned int> tab_pos;
            auto this_read_index{myutility::split(index_line, "\t")};
            std::string read_name{this_read_index[0]};
            std::string_view start_sv{this_read_index[1]};
            std::string_view stop_sv{this_read_index[2]};
            std::from_chars(start_sv.data(), start_sv.data() + start_sv.size(), start);
            std::from_chars(stop_sv.data(), stop_sv.data() + stop_sv.size(), stop);
            reads_index[read_name].push_back(start);
            reads_index[read_name].push_back(stop);
        }
        input_file_index.close();

        char read_name[512];
        for (const std::string& id : read_names) {
            std::string key{id.size() <= used_key_length ? id : id.substr(0, used_key_length)};
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
                    for (size_t idx{strlen(read_name)}; idx < stop_offset - start_offset; idx++) {
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
    if (exists(index_file_path_prefix.concat(".index")) &&
        last_write_time(index_file_path_prefix.concat(".index")) >
        last_write_time(std::filesystem::path{m_input_file})) {
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

void FastqReader::index(unsigned key_len)
{
    if (m_input_file.ends_with(".gz")) {
        index_fastq_gz(key_len);
        return;
    }
    std::filesystem::path input_file{m_input_file};
    std::filesystem::path input_file_idx{input_file.concat(".index")};
    if (exists(input_file_idx)) {
        if (last_write_time(input_file) > last_write_time(input_file_idx)) {
            // input_file is newer than index
            index_fastq(key_len);
        }
    } else {
        index_fastq(key_len);
    }
}

void FastqReader::index_fastq(unsigned key_len)
{
    std::ifstream infile_text{m_input_file.data(), std::ios::in};
    auto index_file{std::filesystem::path{m_input_file}.concat(".index")};
    std::ofstream output_index_stream{index_file, std::ios::out};
    if (!output_index_stream) {
        std::cerr << REDS + fmt::format("Failed when opened file: {}", index_file.c_str()) + COLOR_END <<
            std::endl;
        exit(1);
    }
    std::unordered_map<std::string, std::vector<size_t>> reads_index{};
    infile_text.seekg(std::ios::beg);
    std::string id;
    int line_number{1};
    size_t start{0}; // the start index in file, [include]
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
            // readName\tstartIndex\tendIndex
            length = 0;
            start += length;
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
    if (!m_input_file.ends_with(".gz")) {
        std::cerr << REDS + fmt::format("Error: input file path {} must ends with .gz when use index_fastq_gz",
                                        m_input_file) + COLOR_END << std::endl;
        exit(1);
    }
    // check_compress_type(std::string{output_file_path});
    if (nanobgzip::check_compress_type(std::string{m_input_file}) != nanobgzip::GzipType::NANO_B_GZIP) {
        std::cerr << REDS <<
            "Error: couldn't index common gzip file. \nYou can refer to the following command to convert it into a NanoBgzip file and build an index at the same time"
            << COLOR_END << endl;
        std::cerr << "\nzcat input.fastq.gz | nanofq compress - output.fastq.gz\n";
        std::cerr <<
            "#The above command will turn comman input.fastq.gz into NanoBgzip output.fastq.gz and create index file simultaneously";
        std::cerr <<
            "Or use bgzip to compress the text file and samtools fqidx to build index. Refer the following command:\n";
        std::cerr << "zcat input.fastq.gz | bgzip -c > output.fastq.gz && samtools fqidx output.fastq.gz" << std::endl;
        exit(1);
    }
    nanobgzip::build_index(std::string{m_input_file}, key_len);
}
