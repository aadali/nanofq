#ifndef NANOFQ_UTILITY_H
#define NANOFQ_UTILITY_H

#include <string>
#include <sstream>
#include <string_view>
#include <vector>
#include <iostream>
#include <fmt/core.h>
#include "AlignmentConfig.h"
#include "AlignmentResult.h"
#include "SequenceInfo.h"

#define DEFAULT_INT std::numeric_limits<int>::max()
#define DEFAULT_FLOAT 3.14f
using namespace std;

const std::string REDS = "\033[1;31m";
const std::string COLOR_END =  "\033[0m";
const std::string WARNS = "\033[1;33m";

struct trim_direction {
    bool trim_top5end{false};
    bool trim_top3end{false};
    bool trim_bot5end{false};
    bool trim_bot3end{false};
};

struct myUtility {
    static std::string rev_com(const std::string& seq);

    [[nodiscard]] static vector<string_view> split(string_view str, string_view delim);

    template <typename T>
    static std::string join(const std::string& separator, const std::vector<T>& list){
        std::stringstream oos;
        for (int i {0}; i<list.size(); i++){
            if (i == list.size()-1){
                oos << list[i];
            } else {
                oos << list[i] << separator;
            }
        }
        return oos.str();
    };

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
                                     float bot3end_identity);

    static trim_direction how_trim(const SequenceInfo& seq_info);

    static std::string get_all_seq_info();

    static std::string check_file(const std::string& parameter, const std::string& file, bool need_directory);
    template <typename T>
    static T check_number(const std::string& parameter, const std::string& number, T min, T max);

    static std::string check_one_candidate(const std::string& parameter,
                                 const std::string& candidate,
                                 const std::vector<std::string>& right_candidate);
};

#endif
