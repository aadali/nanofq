#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#include <string>
#include <fmt/core.h>
#include "argparse.h"
#include "myUtility.h"

const int MINT = 1;
const int MAXT = 16;
const int MINC = 10000;
const int MAXC = 100000;
const int MINL = 1;
const int MAXL = std::numeric_limits<int>::max();
const double MIN_PERCENT = 0.0;
const double MAX_PERCENT = 1.0;
const int MINB = 1;
const int MAX24B = 24;
const int MAX96B = 96;
const int MIN_TARGET = 10;
const int MAX_TARGET = 2000;

void check_file(const std::string& file, bool need_directory);

template <typename T>
void check_number_in_range(const std::string& parameter, const T& number, T min, T max,
                           argparse::ArgumentParser& parser, bool integer) {
    std::string type{integer ? "integer" : "float"};
    if (number < min || number > max) {
        std::cerr << fmt::format("{} should be a {}, and in range ({}, {})", parameter, type, min, max) << std::endl;
        std::cerr << parser << std::endl;
        exit(1);
    }
}

template <typename T>
void check_choices(const std::string& parameter, std::vector<T>& choices, std::vector<T>& allowed_choices,
                   argparse::ArgumentParser& parser) {
    for (T& candidate : choices) {
        if (std::find(allowed_choices.begin(), allowed_choices.end(), candidate) == allowed_choices.end()) {
            std::cerr << fmt::format("{} allowed choice should be in [{}]", parameter,
                                     myUtility::join<std::string>(", ", allowed_choices)) << '\n';
            std::cerr << parser << std::endl;
            exit(1);
        }
    }
}

argparse::ArgumentParser& get_arguments(int argc, char* argv[]);
#endif //PARSE_ARGUMENTS_H