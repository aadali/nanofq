#include "myUtility.h"
#include "SequenceInfo.h"
#include "Adapter.h"

namespace myutility {
std::string rev_com(const std::string& seq)
{
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

[[nodiscard]] std::vector<std::string_view> split(std::string_view str, std::string_view delim)
{
    std::vector<std::string_view> result;
    size_t pos = 0;
    while ((pos = str.find(delim)) != std::string::npos) {
        std::string_view token{str.substr(0, pos)};
        result.push_back(token);
        str.remove_prefix(pos + delim.size());
    }
    result.push_back(str);
    return result;
}


std::string_view get_read_name_prefix(std::string_view header, unsigned key_length)
{
    size_t space_idx{header.find(' ')};
    if (space_idx == std::string::npos) {
        return header.size() < key_length + 1 ? header.substr(1) : header.substr(1, key_length);
    }
    return space_idx <= key_length ? header.substr(1, space_idx - 1) : header.substr(1, key_length);
}

void smith_waterman(std::string_view target_seq, std::string_view query_seq, AlignmentConfig& config,
                               AlignmentResult& result)
{
    if (!result.is_empty()) {
        std::cerr << REDS + "AlignmentResult should be empty" + COLOR_END << std::endl;;
        exit(1);
    }
    if (target_seq.size() > MAX_TARGET_LEN) {
        std::cerr << REDS + fmt::format("max length of target seq should less than {}", MAX_TARGET_LEN) + COLOR_END <<
            std::endl;;
        exit(1);
    }
    if (query_seq.size() > MAX_QUERY_LEN) {
        std::cerr << REDS + fmt::format("max length of query seq should less than {}", MAX_QUERY_LEN) + COLOR_END <<
            std::endl;;
        exit(1);
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


void update_sequence_info(SequenceInfo& seq_info, int top5end_len, float top5end_percent,
                                     float top5end_identity, int top3end_len, float top3end_percent,
                                     float top3end_identity, int bot5end_len, float bot5end_percent,
                                     float bot5end_identity, int bot3end_len, float bot3end_percent,
                                     float bot3end_identity)
{
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

trim_direction how_trim(const SequenceInfo& seq_info)
{
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

std::string get_all_seq_info()
{
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
}
