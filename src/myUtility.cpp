#include <algorithm>
#include <ranges>
#include <regex>
#include <sstream>
#include <filesystem>
#include "myUtility.h"
#include "Adapter.h"

std::string myUtility::rev_com(const std::string& seq) {
    std::string sequence{seq};
    std::ranges::transform(sequence,
                           std::begin(sequence),
                           [](const char& c){
                               switch (c) {
                               case 'A':
                                   return 'T';
                               case 'T':
                                   return 'A';
                               case 'G':
                                   return 'C';
                               case 'C':
                                   return 'G';
                               case 'V':
                                   return 'B';
                               default:
                                   return 'N';
                               }
                           });

    std::ranges::reverse(sequence);
    return sequence;
}

[[nodiscard]] vector<string_view> myUtility::split(string_view str, string_view delim) {
    vector<string_view> result;
    size_t pos = 0;
    while ((pos = str.find(delim)) != string::npos) {
        string_view token{str.substr(0, pos)};
        result.push_back(token);
        str.remove_prefix(pos + delim.size());
    }
    result.push_back(str);
    return result;
}

template <typename... Args>
std::string myUtility::join(const std::string& separator, Args... args) {
    std::ostringstream oss;
    ((oss << args << separator), ...);
    std::string result = oss.str();
    if (!result.empty()) {
        result.erase(result.size() - separator.size());
    }
    return result;
}

template <typename T>
std::string myUtility::join(const std::string& separator, const std::vector<T>& v) {
    std::stringstream oss;
    for (int i{0}; i < v.size(); i++) {
        if (i == v.size() - 1) {
            oss << v[i];
        } else {
            oss << v[i] << separator;
        }
    }
    return oss.str();
}

string_view myUtility::get_read_name_prefix(string_view header, unsigned key_length) {
    size_t space_idx{header.find(' ')};
    if (space_idx == std::string::npos) {
        return header.size() < key_length + 1 ? header.substr(1) : header.substr(1, key_length);
    } else {
        return space_idx <= key_length ? header.substr(1, space_idx - 1) : header.substr(1, key_length);
    }
}

void myUtility::smith_waterman(string_view target_seq, string_view query_seq, AlignmentConfig& config,
                               AlignmentResult& result) {
    if (!result.is_empty()) {
        throw std::runtime_error("AlignmentResult should be empty");
    }
    if (target_seq.size() > MAX_TARGET_LEN) {
        throw std::runtime_error(fmt::format("max length of target seq should less than {}", MAX_TARGET_LEN));
    }
    if (query_seq.size() > MAX_QUERY_LEN) {
        throw std::runtime_error(fmt::format("max length of query seq should less than {}", MAX_QUERY_LEN));
    }
    int diagonal_score{0}, up_score{0}, left_score{0};
    for (int row{1}; row < query_seq.size(); row++) {
        for (int col{1}; col < target_seq.size(); col++) {
            /* calculate score from diagonal */
            // int score{target_seq[col - 1] == query_seq[row - 1] ? config.m_match : config.m_mismatch};
            int score{
                target_seq[col - 1] == query_seq[row - 1] ||
                (query_seq[row - 1] == 'V' && (target_seq[col - 1] == 'G' || target_seq[col - 1] == 'A' || target_seq[
                    col
                    - 1] == 'C')) ||
                (query_seq[row - 1] == 'B' && (target_seq[col - 1] == 'G' || target_seq[col - 1] == 'T' || target_seq[
                    col
                    - 1] == 'C'))
                    ? config.m_match
                    : config.m_mismatch
            };
            diagonal_score = config.get_score(row - 1, col - 1) + score;

            /* calculate score from left
             * 1 in direction matrix mean from left */
            score = config.get_direction(row, col - 1) == AlignmentConfig::Direction::Left
                        ? config.m_gap_extend
                        : config.m_gap_open;
            left_score = config.get_score(row, col - 1) + score;

            /* calculate score from up
             * 2 in direction matrix mean from up*/
            score = config.get_direction(row - 1, col) == AlignmentConfig::Direction::Up
                        ? config.m_gap_extend
                        : config.m_gap_open;
            up_score = config.get_score(row - 1, col) + score;

            // update the value of this cell
            config.set_score(row, col, std::max({diagonal_score, left_score, up_score, 0}));

            // record which direction the value of this cell is from
            if (config.get_score(row, col) == diagonal_score) {
                config.set_direction(row, col, AlignmentConfig::Direction::Diag);
            } else if (config.get_score(row, col) == left_score) {
                config.set_direction(row, col, AlignmentConfig::Direction::Left);
            } else {
                config.set_direction(row, col, AlignmentConfig::Direction::Up);
            }

            // update the max score and the row, column index of max score
            if (config.get_score(row, col) > result.get_max_score()) {
                result.set_max_score(config.get_score(row, col));
                result.set_stop_idx(row, col);
            }
        }
    }

    // trace back depend on the direction matrix
    auto [row, col] = result.get_stop_idx();
    while (row > 0 && col > 0 && config.get_score(row, col) > 0) {
        if (config.get_direction(row, col) == AlignmentConfig::Direction::Diag) {
            // '|' mean match, ':' mean mismatch
            result.push_back(target_seq[col - 1],
                             query_seq[row - 1],
                             target_seq[col - 1] == query_seq[row - 1] ? '|' : ':');
            --row;
            --col;
        } else if (config.get_direction(row, col) == AlignmentConfig::Direction::Left) {
            // ' ' mean gap
            result.push_back(target_seq[col - 1], '-', ' ');
            --col;
        } else {
            result.push_back('-', query_seq[row - 1], ' ');
            --row;
        }
    }
    result.reverse_align();
    result.set_start_idx(row, col);
}


void myUtility::update_sequence_info(SequenceInfo& seq_info, int top5end_len, float top5end_percent,
                                     float top5end_identity, int top3end_len, float top3end_percent,
                                     float top3end_identity, int bot5end_len, float bot5end_percent,
                                     float bot5end_identity, int bot3end_len, float bot3end_percent,
                                     float bot3end_identity) {
    /*
     * parameter: *len, *percent, *identity should be checked at the stage of parsing argument from CLI
     * all of them couldn't be negative
     * *percent, *identity must >0 and <1.0
     * if parameter not be set,
        default *len = std::numeric_limits<int>::max()
        default *percent = DEFAULT_FLOAT // 3.14
        default *identity = DEFAULT_FLOAT // 3.14
    */
    if (top5end_len != DEFAULT_INT) {
        get<0>(seq_info.m_top5end) = top5end_len;
    }
    if (static_cast<int>(top5end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(seq_info.m_top5end) = top5end_percent;
    }
    if (static_cast<int>(top5end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(seq_info.m_top5end) = top5end_identity;
    }

    if (top3end_len != DEFAULT_INT) {
        get<0>(seq_info.m_top3end) = top3end_len;
    }
    if (static_cast<int>(top3end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(seq_info.m_top3end) = top3end_percent;
    }
    if (static_cast<int>(top3end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(seq_info.m_top3end) = top3end_identity;
    }


    if (bot5end_len != DEFAULT_INT) {
        get<0>(seq_info.m_bot5end) = bot5end_len;
    }
    if (static_cast<int>(bot5end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(seq_info.m_bot5end) = bot5end_percent;
    }
    if (static_cast<int>(bot5end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(seq_info.m_bot5end) = bot5end_identity;
    }

    if (bot3end_len != DEFAULT_INT) {
        get<0>(seq_info.m_bot3end) = bot3end_len;
    }
    if (static_cast<int>(bot3end_percent) != static_cast<int>(DEFAULT_FLOAT)) {
        get<1>(seq_info.m_bot3end) = bot3end_percent;
    }
    if (static_cast<int>(bot3end_identity) != static_cast<int>(DEFAULT_FLOAT)) {
        get<2>(seq_info.m_bot3end) = bot3end_identity;
    }
}

trim_direction myUtility::how_trim(const SequenceInfo& seq_info) {
    trim_direction td;
    if (!seq_info.m_top5end_query.empty() && get<0>(seq_info.m_top5end) > 0) {
        td.trim_top5end = true;
    }
    if (!seq_info.m_top3end_query.empty() && get<0>(seq_info.m_top3end) > 0) {
        td.trim_top3end = true;
    }
    if (!seq_info.m_bot5end_query.empty() && get<0>(seq_info.m_bot5end) > 0) {
        td.trim_bot5end = true;
    }
    if (!seq_info.m_bot3end_query.empty() && get<0>(seq_info.m_bot3end) > 0) {
        td.trim_bot3end = true;
    }
    return td;
}

std::string myUtility::get_all_seq_info() {
    std::stringstream info;
    std::unordered_map<std::string, SequenceInfo> all_trim_info = barcode_info::get_trim_info();
    for (std::string kit : std::vector{"SQK-LSK114", "SQK-ULK114", "SQK-RAD114", "SQK-PCS114"}) {
        SequenceInfo& seq_info{all_trim_info.at(kit)};
        info << seq_info.seq_info() << '\n';
    }
    for (std::string kit : std::vector{"SQK-NBD114.24", "SQK-RBK114.24", "SQK-PCB114.24"}) {
        for (int i{0}; i < 24; i++) {
            info << all_trim_info.at(fmt::format("{}-{}", kit, i + 1)).seq_info() << '\n';
        }
    }

    for (std::string kit : std::vector{"SQK-NBD114.96", "SQK-RBK114.96"}) {
        for (int i{0}; i < 24; i++) {
            info << all_trim_info.at(fmt::format("{}-{}", kit, i + 1)).seq_info() << '\n';
        }
    }
    return info.str();
}

std::string myUtility::check_file(const std::string& parameter, const std::string& file, bool need_directory) {
    std::filesystem::path fp{file};
    if (!exists(fp)) {
        std::cerr << fmt::format("Error parameter: {}. No such file or directory for {}", parameter, file) << std::endl;
        exit(1);
    }
    if (need_directory && !is_directory(fp)) {
        std::cerr << fmt::format("{} is file. Directory needed", file) << std::endl;
        exit(1);
    }
    return file;
}

template <typename T>
T myUtility::check_number(const std::string& parameter, const std::string& number, T min, T max) {
    std::regex pat{R"([-+]?\d*\.?\d+)"};
    if (bool mat{std::regex_match(number, pat)}; mat) {
        T n{static_cast<T>(stof(number))};
        if (n < min || n > max) {
            std::cerr << fmt::format("Error parameter: {}. The range should be ({}, {})", parameter, min, max);
            exit(1);
        }
        return n;
    }
    std::cerr << fmt::format("Error parameter: {}. this parameter should be number") << endl;
    exit(1);
}

std::string myUtility::check_one_candidate(const std::string& parameter,
                                           const std::string& candidate,
                                           const std::vector<std::string>& right_candidate) {
    for (const std::string& item : right_candidate) {
        if (item == candidate) {
            return candidate;
        }
    }
    std::cerr << fmt::format("Error parameter: {}. {} is illegal, allowed options are [{}]", parameter, candidate,
                             join<std::string>(", ", right_candidate)) << std::endl;
    exit(1);
}
