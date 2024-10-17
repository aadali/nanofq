cmake_minimum_required(VERSION 3.28)
project(nanofq)

set(CMAKE_CXX_STANDARD 20)
add_executable(nanofq
        main.cpp
        src/Read.hpp
        src/Timer.hpp
        src/FastqReader.h
        src/FastqReader.cpp
        src/Work.h
        src/Work.cpp
        src/utility.h
        src/kseq.h
)

target_link_libraries(nanofq fmt)
target_link_libraries(nanofq z)
include_directories(./src)
