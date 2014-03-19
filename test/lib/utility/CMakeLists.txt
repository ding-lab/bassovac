cmake_minimum_required(VERSION 2.8)

set(TEST_LIBS utility ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
include_directories(${GTEST_INCLUDE_DIRS})

def_test(Binomial)
def_test(CountingSort)
def_test(Lut)
