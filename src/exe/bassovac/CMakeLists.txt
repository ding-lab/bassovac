cmake_minimum_required(VERSION 2.8)

project(bassovac_app)

set(SOURCES 
    BassovacApp.cpp
    BassovacApp.hpp
    main.cpp
)

set(EXECUTABLE_NAME bassovac)
add_executable(${EXECUTABLE_NAME} ${SOURCES})
target_link_libraries(${EXECUTABLE_NAME} bvprob io utility ${Boost_PROGRAM_OPTIONS_LIBRARY})
set_target_properties(${EXECUTABLE_NAME} PROPERTIES PACKAGE_OUTPUT_NAME ${EXECUTABLE_NAME}${EXE_VERSION_SUFFIX})
install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin/)
