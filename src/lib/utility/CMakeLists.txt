cmake_minimum_required(VERSION 2.8)

set(SOURCES
    Binomial.cpp
    Binomial.hpp
    CountingSort.hpp
    Lut.cpp
    Lut.hpp
    TempFile.hpp
)

add_library(utility ${SOURCES})
