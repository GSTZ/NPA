cmake_minimum_required(VERSION 4.1)
project(npaOdd)

set(CMAKE_CXX_STANDARD 26)

find_package(OpenMP REQUIRED)

add_executable(npaOdd main.cpp
        CG.cpp
        input.cpp
        npa.cpp
        CGC.h)

target_compile_options(npaOdd PRIVATE -O3 -fopenmp)
target_link_options(npaOdd PRIVATE -fopenmp)