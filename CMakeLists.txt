cmake_minimum_required(VERSION 2.8)

project(bassovac)

set(CMAKE_INSTALL_PREFIX "/usr")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build-common/cmake)
include(TestHelper)
include(VersionHelper)
include(ProjectHelper)
include(PlatformChecks)
include(BuildDeps)

configure_file(version/version.h.in version/version.h @ONLY)
include_directories(${PROJECT_BINARY_DIR}/version)

add_custom_target(deps ALL)

set(REQUIRED_BOOST_LIBS program_options filesystem system iostreams)
set(BOOST_URL ${CMAKE_SOURCE_DIR}/vendor/boost)
set(SAMTOOLS_URL ${CMAKE_SOURCE_DIR}/vendor/samtools-0.1.19)
add_definitions("-DBOOST_CHRONO_HEADER_ONLY")
build_boost(${BOOST_URL} ${CMAKE_BINARY_DIR}/vendor/boost ${REQUIRED_BOOST_LIBS})
build_samtools(${SAMTOOLS_URL} ${CMAKE_BINARY_DIR}/vendor/samtools)

message("-- Boost include directory: ${Boost_INCLUDE_DIRS}")
message("-- Boost libraries: ${Boost_LIBRARIES}")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${Samtools_INCLUDE_DIRS})

find_package(Threads REQUIRED)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall -pedantic -fPIC")
endif()

if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -fPIC ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -no-prec-div -fp-model fast=1 -ansi-alias")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -no-prec-div -fp-model fast=1 -ansi-alias -g")
else()
    message(
        WARNING
        "Unknown compiler: ${CMAKE_C_COMPILER_ID}, we'll see how this goes..."
    )
endif()

# Find C++11 flags
find_cxx11_flags(CXX11_FLAGS CXX11_FOUND)
if(NOT CXX11_FOUND)
    message(FATAL_ERROR
        "Failed to find C++11 support with compiler:
         ${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}).")
endif(NOT CXX11_FOUND)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX11_FLAGS}")

# make sure to pick up headers from library dirs
include_directories("src/lib")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# unit tests
find_package(GTest)
if(GTEST_FOUND)
    message("Google Test framework found, building unit tests")
    enable_testing(true)
    add_projects(test/lib)
else()
    message("Google Test framework not found, no tests will be built")
    message("GTest is available at http://code.google.com/p/googletest/")
    message("Ubuntu users can likely sudo apt-get install gtest")
endif()

# main project
add_projects(src/lib)
add_projects(src/exe)

include(PackageInfo.cmake)
add_subdirectory(integration-test)
