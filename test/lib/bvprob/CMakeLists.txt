cmake_minimum_required(VERSION 2.8)

project(test_bvprob)

set(TEST_LIBS bvprob io utility ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

#def_test(Bassovac)
def_test(ExpectedResult)
def_test(FastaReader)
def_test(PBin)
def_test(Sample)
