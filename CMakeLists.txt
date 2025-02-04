cmake_minimum_required(VERSION "3.18.0")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Name of project and that it is C++ only.
project (impact CXX)


# Dependencies
find_package(ZLIB REQUIRED)
#find_package(Armadillo REQUIRED 9.800.4)
find_package(Threads REQUIRED)

# Statically link program and depedencies
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/lib")
include_directories("${PROJECT_SOURCE_DIR}/lib/bamtools/src")
include_directories("${PROJECT_SOURCE_DIR}/lib/seqan")

# Build Setup
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_definitions (${SEQAN_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -pthread -std=c++14")

file(GLOB IMPACT_FILES
    ${PROJECT_SOURCE_DIR}/src/impact.cpp
)


add_subdirectory(lib)

# TIDDIT executable
add_executable(impact ${IMPACT_FILES})
target_link_libraries (impact ${SEQAN_LIBRARIES})
target_link_libraries(impact ${ARMADILLO_LIBRARIES})
target_link_libraries(impact ${ZLIB_LIBRARIES})
target_link_libraries(impact  BamTools)
