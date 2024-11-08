#ifndef PARSE_ARGUMENTS_H
#define PARSE_ARGUMENTS_H

#include <string>
#include "argparse.h"

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
void check_number_in_range(const std::string& parameter, const T& number, T min, T max, argparse::ArgumentParser&parser, bool integer);

template <typename T>
void check_choices(const std::string& parameter, std::vector<T>& choices, std::vector<T>& allowed_choices,
                   argparse::ArgumentParser& parser);

argparse::ArgumentParser& get_arguments(int argc, char* argv[]);
#endif //PARSE_ARGUMENTS_H
