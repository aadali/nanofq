#pragma once

#include <string>
#include <iostream>
#include <chrono>

class Timer {
private:
    std::string mMessage;
    std::chrono::time_point<std::chrono::steady_clock> mStart{std::chrono::steady_clock::now()};

public:
    Timer() = default;

    explicit Timer(std::string name) : mMessage(std::move(name)) {}

    ~Timer() {
        std::cout
                << mMessage
                << ": "
                << std::chrono::duration_cast<std::chrono::duration<double>>(
                        std::chrono::steady_clock::now() - mStart).count()
                << 's'
                << std::endl;
    }
};

