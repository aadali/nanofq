#ifndef NANOFQ_UTILITY_H
#define NANOFQ_UTILITY_H

#include <string>
#include <string_view>
#include <vector>
#include <fmt/core.h>
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
#include "SequenceInfo.h"

#define DEFAULT_INT std::numeric_limits<int>::max()
#define DEFAULT_FLOAT 3.14f
using namespace std;

struct trim_direction {
    bool trim_top5end{false};
    bool trim_top3end{false};
    bool trim_bot5end{false};
    bool trim_bot3end{false};
};

struct utility {
    static std::string rev_com(const std::string& seq);

    [[nodiscard]] static vector<string_view> split(string_view str, string_view delim);

    static string_view get_read_name_prefix(string_view header, unsigned key_length);

    static void smith_waterman(string_view target_seq, string_view query_seq, AlignmentConfig& config,
                               AlignmentResult& result);

    static void update_sequence_info(SequenceInfo& seq_info,
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
                                     float bot3end_identity
    );

    static trim_direction how_trim(const SequenceInfo& seq_info);
};

#endif
