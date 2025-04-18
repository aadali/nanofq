cmake_minimum_required(VERSION 3.28)
project(nanofq)

set(CMAKE_CXX_STANDARD 20)
add_executable(nanofq
        main.cpp
        src/Read.h
        src/Read.cpp
        src/Timer.h
        src/FastqReader.h
        src/FastqReader.cpp
        src/Work.h
        src/Work.cpp
        src/myUtility.h
        src/myUtility.cpp
        src/kseq.h
        src/AlignmentConfig.h
        src/AlignmentConfig.cpp
        src/AlignmentResult.h
        src/AlignmentResult.cpp
        src/Adapter.h
        src/SequenceInfo.h
        src/SequenceInfo.cpp
        src/ArgumentParse.cpp
        src/ArgumentParse.h
        src/submain.h
        src/ThreadPool.h
        src/NanoBgzip.cpp
        src/NanoBgzip.h
)

#target_link_libraries(nanofq fmt)
target_link_libraries(nanofq z)
include_directories(./src)
add_subdirectory(fmt)
target_link_libraries(nanofq fmt::fmt)

add_custom_command(TARGET nanofq POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_SOURCE_DIR}/scripts/plot.py"
        "$<TARGET_FILE_DIR:nanofq>/plot.py"
        COMMENT "Copying plot.py to the binary directory"
)
