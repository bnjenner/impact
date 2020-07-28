cmake_minimum_required(VERSION "3.18.0")

# Name of project and that it is C++ only.
project (impact CXX)

# Supplementary Files
file(GLOB sources "utils/*.cpp")

# Dependencies
find_package (ZLIB)
find_package (SeqAn REQUIRED)
find_package(Armadillo REQUIRED)

# Build Setup
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -pthread -std=c++14")

add_executable (impact main.cpp)
target_link_libraries (impact ${SEQAN_LIBRARIES})
target_link_libraries(impact ${ARMADILLO_LIBRARIES})
