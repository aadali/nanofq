#ifndef FASTQREADER_H
#define FASTQREADER_H

#include <string_view>
#include <vector>
#include <fstream>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <optional>
#include <zlib.h>
#include "Read.h"
#include "kseq.h"
const std::string finished_read_name{"FINISHED FINISHED FINISHED"};
KSEQ_INIT(gzFile, gzread)
using shared_vec_reads = std::shared_ptr<std::vector<std::shared_ptr<Read>>>;
using shared_read = std::shared_ptr<Read>;

class FastqReader
{
private:
    std::string_view m_input_file;
    gzFile m_infile_gz{nullptr};
    kseq_t* m_seq;
    int m_chunk_size;
    char *m_buffer;
    bool m_finish{false};

public:
    FastqReader() = delete;

    FastqReader(std::string_view input_file, int chunk);

    FastqReader(const FastqReader &) = delete;

    FastqReader(FastqReader &&) = delete;

    FastqReader &operator=(const FastqReader &) = delete;

    FastqReader &operator=(FastqReader &&) = delete;

    ~FastqReader();

    std::shared_ptr<std::vector<Read>> read_chunk_fastq();

    Read read_one_fastq();

    inline bool read_finish() const { return m_finish; };

    void find_reads(const std::string &input_reads, std::ostream &out, bool use_index, unsigned key_length = 8);
    void index(unsigned key_len);

private:
    static std::unordered_set<std::string> get_searching_read_names(const std::string &input_reads);
    void index_fastq(std::string_view output_file_path, unsigned key_len);
    void index_fastq_gz(std::string_view output_file_path, unsigned key_len);
};

#endif // FASTQREADER_H
