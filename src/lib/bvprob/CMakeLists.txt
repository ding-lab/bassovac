cmake_minimum_required(VERSION 2.8)

project(bvprob_lib)

set(SOURCES
    Bassovac.cpp
    Bassovac.hpp
    Fasta.cpp
    Fasta.hpp
    FastaReader.cpp
    FastaReader.hpp
    IOError.hpp
    PBin.cpp
    PBin.hpp
    ResultFormatter.cpp
    ResultFormatter.hpp
    Sample.cpp
    Sample.hpp
    Tokenizer.hpp
)

add_library(bvprob ${SOURCES})
target_link_libraries(bvprob m ${Samtools_LIBRARIES} ${Boost_LIBRARIES})
