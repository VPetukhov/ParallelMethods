cmake_minimum_required(VERSION 3.0)
project(ParallelMethids)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

find_package(Boost 1.54.0 COMPONENTS unit_test_framework REQUIRED)

set(INCLUDE_DIRS ${Boost_INCLUDE_DIR} ${PROJECT_SOURCE_DIR} Utils)
include_directories(${INCLUDE_DIRS})

set(LIBRARIES ${Boost_LIBRARIES})

FILE(GLOB UTIL_SOURCES Utils/*.cpp)
add_library(Utils ${UTIL_SOURCES})

set(SOURCE_FILES main.cpp)
add_executable(Smirnov3 ${SOURCE_FILES})

target_link_libraries(Smirnov3 Utils ${LIBRARIES})

set (TESTS Tests)
add_executable (${TESTS} ${PROJECT_SOURCE_DIR}/Tests.cpp)
target_link_libraries (${TESTS} Utils ${LIBRARIES})
add_test (${TESTS} ${TESTS})