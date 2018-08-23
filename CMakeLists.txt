cmake_minimum_required(VERSION 2.8)

project(MatrixChain
        LANGUAGES CXX)

# Add any platform-specific and global build settings here:
#
include_directories(.)

add_executable(mattest
    test.cpp
    matrix.cpp
    )
