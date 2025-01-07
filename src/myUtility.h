#ifndef NANOFQ_UTILITY_H
#define NANOFQ_UTILITY_H

#include <string>
#include <optional>
#include <filesystem>
#include <string_view>
#include <vector>
#include <iostream>
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
class SequenceInfo;

#define DEFAULT_INT std::numeric_limits<int>::max()
#define DEFAULT_FLOAT 3.14f

const std::string REDS = "\033[1;31m";
const std::string COLOR_END = "\033[0m";
const std::string WARNS = "\033[1;33m";

struct trim_direction {
    bool trim_top5end{false};
    bool trim_top3end{false};
    bool trim_bot5end{false};
    bool trim_bot3end{false};
};

namespace myutility {
    std::string rev_com(const std::string& seq);

    [[nodiscard]] std::vector<std::string_view> split(
        std::string_view str,
        std::string_view delim);

    std::string_view get_read_name_prefix(
        std::string_view header,
        unsigned key_length);

    std::optional<std::vector<std::filesystem::path>> get_fastqs(
        const std::string& input_path);

    void smith_waterman(
        std::string_view target_seq,
        std::string_view query_seq,
        AlignmentConfig& config,
        AlignmentResult& result);

    void update_sequence_info(
        SequenceInfo& seq_info,
        int top5end_len,
        float top5end_percent,
        float top5end_identity,
        int top3end_len,
        float top3end_percent,
        float top3end_identity,
        int bot5end_len,
        float bot5end_percent,
        float bot5end_identity,
        int bot3end_len,
        float bot3end_percent,
        float bot3end_identity);

    trim_direction how_trim(const SequenceInfo& seq_info);

    std::string get_all_seq_info();

    std::string check_file(
        const std::string& parameter,
        const std::string& file,
        bool need_directory);

    template <typename T>
    T check_number(
        const std::string& parameter,
        const std::string& number,
        T min,
        T max);

    std::string check_one_candidate(
        const std::string& parameter,
        const std::string& candidate,
        const std::vector<std::string>& right_candidate);
};

#endif
