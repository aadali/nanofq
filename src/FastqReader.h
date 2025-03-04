#ifndef FASTQREADER_H
#define FASTQREADER_H

#include <string_view>
#include <vector>
#include <fstream>
#include <memory>
#include <zlib.h>
#include <filesystem>
#include <cstdio>
#include <charconv>
#include <optional>
#include "fmt/core.h"
#include "Read.h"
#include "NanoBgzip.h"
#include "kseq.h"
#include "fmt/chrono.h"
const std::string finished_read_name{"FINISHED FINISHED FINISHED"};
KSEQ_INIT(gzFile, gzread)
using shared_vec_reads = std::shared_ptr<std::vector<std::shared_ptr<Read>>>;
using shared_read = std::shared_ptr<Read>;
constexpr size_t FASTQ_BUFFER_SIZE{1 << 23}; // the longest read length exceeds 4Mb
using nanobgzip_reads_index = std::pair<std::vector<std::pair<size_t, size_t>>,
                                        std::unordered_map<std::string, std::tuple<unsigned, size_t, size_t>>>;

class FastqReader
{
private:
    std::string m_input_file;
    std::string m_input_file_index;
    gzFile m_infile_gz{nullptr};
    kseq_t* m_seq{nullptr};
    int m_chunk_size;
    char* m_buffer;
    bool m_finish{false};
    bool m_is_directory{false};
    std::shared_ptr<std::unordered_map<std::string, std::string>> m_reads_seq;

public:
    FastqReader() = default;

    FastqReader(const std::string& input_file, int chunk, bool is_directory);

    FastqReader(const FastqReader&) = delete;

    FastqReader(FastqReader&&) = delete;

    FastqReader& operator=(const FastqReader&) = delete;

    FastqReader& operator=(FastqReader&&) = delete;

    ~FastqReader();

    std::shared_ptr<std::vector<Read>> read_chunk_fastq();

    Read read_one_fastq();

    std::optional<std::vector<std::filesystem::path>> get_fastqs() const;

    inline bool read_finish() const { return m_finish; };

    void index(bool force_index);

    void find(
        const std::string& input_reads,
        std::ostream& out,
        bool use_index);

    static Read fastq_record_ok(int l, kseq_t* seq, const char* file);

private:
    static std::unordered_set<std::string> get_searching_read_names(const std::string& input_reads);

    void search_read_one_by_one(
        std::unordered_set<std::string>& read_names,
        std::ostream& out);

    void index_fastq();

    void index_fastq_gz();

    std::unordered_map<std::string, std::pair<size_t,size_t>> read_index() const;

    nanobgzip_reads_index read_gz_index() const;

    void find_reads(
        const std::string& input_reads,
        std::ostream& out,
        bool use_index);

    void find_reads_in_gz(
        const std::string& input_reads,
        std::ostream& out,
        bool use_index);
};

#endif // FASTQREADER_H
