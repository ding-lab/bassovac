project(breakdancer)

set(SOURCES
    BamEntry.cpp
    BamEntry.hpp
    BamFilter.hpp
    BamIntersector.cpp
    BamIntersector.hpp
    BamReaderBase.cpp
    BamReaderBase.hpp
    BamReader.hpp
    CigarParser.cpp
    CigarParser.hpp
    PileupBuffer.cpp
    PileupBuffer.hpp
    Pileup.cpp
    Pileup.hpp
    RegionLimitedBamReader.hpp
    SamConvert.cpp
    SamConvert.hpp
)

add_library(io ${SOURCES})
target_link_libraries(io ${Samtools_LIBRARIES} z m)
