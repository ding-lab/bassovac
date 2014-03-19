cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS io utility ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

def_test(BamEntry)
def_test(BamIntersector)
def_test(BamReader)
def_test(CigarParser)
def_test(Pileup)
def_test(PileupBuffer)
