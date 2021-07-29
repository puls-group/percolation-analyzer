cmake_minimum_required(VERSION 3.15)


project(percolation-analyzer VERSION 0.1.0 DESCRIPTION "" LANGUAGES C CXX)

#include(CTest)
#enable_testing()

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# add_subdirectory(cmake)

if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(BIN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})

#include_directories("${INCLUDE_DIR}")

# include external dependencies
add_subdirectory(${LIB_DIR})

# Add compile scripts for own submodules
add_subdirectory(${SRC_DIR})

add_subdirectory(test/Catch2)

add_executable(test_runner test/testing.cpp)
target_include_directories(test_runner PUBLIC ${INCLUDE_DIR})
target_link_libraries(test_runner PRIVATE Catch2::Catch2 percolation-analyzer-cpp)


add_executable(sample_graph src/sample_graph_builder.cpp)
target_include_directories(sample_graph PUBLIC ${INCLUDE_DIR})
target_link_libraries(sample_graph PRIVATE molecular-graph-cpp)
