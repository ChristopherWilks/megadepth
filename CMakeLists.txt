cmake_minimum_required(VERSION 2.6)
project(bamcount)

set(CMAKE_CXX_STANDARD 14)

# without following line, I get CMake error:
# "SeqAn requires C++14 since v2.2.0, but your compiler does not support it"
set(CXX14_FLAG -std=c++14)

find_package (ZLIB)

# installed with conda install -c conda-forge seqan-library
find_package (SeqAn REQUIRED)
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (${CMAKE_SOURCE_DIR}/docopt)
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_executable(bamcount
        bamcount.cpp
        docopt/docopt.cpp)
target_link_libraries (bamcount ${SEQAN_LIBRARIES})
