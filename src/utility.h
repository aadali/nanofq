#pragma once

#include <string>
#include <string_view>
#include <vector>
#include "kseq.h"
#include "zlib.h"


using namespace std;
namespace utility {
    [[nodiscard]] vector<string_view> split(string_view str, string_view delim) {
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

    string_view get_read_name_prefix(string_view header, unsigned key_length) {
        size_t space_idx{header.find(' ')};
        if (space_idx == std::string::npos) {
            return header.size() < key_length + 1 ? header.substr(1) : header.substr(1, key_length);
        } else {
            return space_idx <= key_length ? header.substr(1, space_idx - 1) : header.substr(1, key_length);
        }
    }
}