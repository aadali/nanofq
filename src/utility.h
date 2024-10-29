#ifndef NANOFQ_UTILITY_H
#define NANOFQ_UTILITY_H

#include <string>
#include <string_view>
#include <vector>
#include <fmt/core.h>
#include "kseq.h"
#include "zlib.h"
#include "AlignmentConfig.h"
#include "AlignmentResult.h"


using namespace std;

struct utility {
    [[nodiscard]] static vector<string_view> split(string_view str, string_view delim) {
        vector<string_view> result;
        size_t pos = 0;
        string_view token;
        while ((pos = str.find(delim)) != string::npos) {
            token = str.substr(0, pos);
            result.push_back(token);
            str.remove_prefix(pos + delim.size());
        }
        result.push_back(str);
        return result;
    }

    static string_view get_read_name_prefix(string_view header, unsigned key_length) {
        size_t space_idx{header.find(' ')};
        if (space_idx == std::string::npos) {
            return header.size() < key_length + 1 ? header.substr(1) : header.substr(1, key_length);
        } else {
            return space_idx <= key_length ? header.substr(1, space_idx - 1) : header.substr(1, key_length);
        }
    }

    static void smith_waterman(const string &target_seq, const string &query_seq, AlignmentConfig &config,
                               AlignmentResult &result) {
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
                int score{target_seq[col - 1] == query_seq[row - 1] ? config.m_match : config.m_mismatch};
                diagonal_score = config.get_score(row - 1, col - 1) + score;

                /* calculate score from left
                 * 1 in direction matrix mean from left */
                score = config.get_direction(row, col - 1) == AlignmentConfig::Direction::Left ? config.m_gap_extend : config.m_gap_open;
                left_score = config.get_score(row, col - 1) + score;


                /* calculate score from up
                 * 2 in direction matrix mean from up*/
                score = config.get_direction(row - 1, col) == AlignmentConfig::Direction::Up ? config.m_gap_extend : config.m_gap_open;
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
};

#endif