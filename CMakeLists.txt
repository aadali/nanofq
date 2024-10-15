cmake_minimum_required(VERSION 3.28)
project(nanofq2)

set(CMAKE_CXX_STANDARD 20)
add_executable(nanofq2
        main.cpp
        src/Read.hpp
        src/Timer.hpp
        src/FastqReader.h
        src/FastqReader.cpp
        src/Work.h
        src/Work.cpp
)

target_link_libraries(nanofq2 fmt)
include_directories(./src)
