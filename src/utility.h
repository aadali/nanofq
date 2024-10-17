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
}