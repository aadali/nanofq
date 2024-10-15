#ifndef FASTQREADER_H
#define FASTQREADER_H

#include <string_view>
#include <vector>
#include <fstream>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <optional>
#include "Read.hpp"

using shared_vec_reads = std::shared_ptr<std::vector<std::shared_ptr<Read>>>;
using shared_read = std::shared_ptr<Read>;

class FastqReader {
private:
    bool m_finish{false};
    std::fstream m_infile;
    shared_vec_reads m_reads{};
    std::string_view m_input_file;
    size_t m_chunk;
    char *m_buffer;
    static std::mutex ms_mtx;
    static std::condition_variable ms_cond;

public:
    FastqReader() = delete;

    FastqReader(std::string_view input_file, size_t chunk);

    FastqReader(const FastqReader &) = delete;

    FastqReader(FastqReader &&) = delete;

    FastqReader &operator=(const FastqReader &) = delete;

    FastqReader &operator=(FastqReader &&) = delete;

    ~FastqReader();

    void read_chunk_fastq();

    bool inline read_finish() const {
        return m_finish;
    };

    std::optional<shared_vec_reads> get_reads();
    void find_reads(const std::string& input_reads, std::ostream& out);
    void index();

private:
    static std::unordered_set<std::string> get_searching_read_names(const std::string& input_reads) ;
};


#endif //FASTQREADER_H
