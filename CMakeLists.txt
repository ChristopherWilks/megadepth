cmake_minimum_required(VERSION 2.6)
set( CMAKE_VERBOSE_MAKEFILE on )
project(megadepth)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMEGADEPTH_VERSION=\"\\\"`cat ../VERSION`\\\"\"")
add_executable(megadepth_dynamic megadepth.cpp)
add_executable(megadepth_static megadepth.cpp)
add_executable(megadepth_statlib megadepth.cpp)
include_directories(htslib libBigWig)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
    set_property(TARGET megadepth_dynamic PROPERTY COMPILE_OPTIONS "-pthread")
    set_property(TARGET megadepth_dynamic PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
    target_link_libraries(megadepth_dynamic "${CMAKE_THREAD_LIBS_INIT}")
endif()

target_link_libraries(megadepth_dynamic z hts BigWig -L${CMAKE_SOURCE_DIR}/htslib -L${CMAKE_SOURCE_DIR}/libBigWig)
#requires static libraries for both zlib and pthread
target_link_libraries(megadepth_static -static ${CMAKE_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SOURCE_DIR}/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/zlib/libz.a -lpthread)
#this build a dynamic binary, but with htslib, libBigWig, and libz statically linked in
target_link_libraries(megadepth_statlib ${CMAKE_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SOURCE_DIR}/libBigWig/libBigWig.a -lz -lcurl -lpthread)
