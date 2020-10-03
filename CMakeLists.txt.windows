cmake_minimum_required(VERSION 2.6)
set( CMAKE_VERBOSE_MAKEFILE on )
project(megadepth)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMEGADEPTH_VERSION=\"\\\"`cat ../VERSION`\\\"\"")
add_executable(megadepth_dynamic megadepth.cpp)
add_executable(megadepth_static megadepth.cpp)
add_executable(megadepth_statlib megadepth.cpp)
include_directories(libdeflate htslib libBigWig)

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
target_link_libraries(megadepth_static -static ${CMAKE_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SOURCE_DIR}/libBigWig/libBigWig.a ${CMAKE_SOURCE_DIR}/zlib/libz.a -lpthread ${CMAKE_SOURCE_DIR}/libdeflate/libdeflate.a)
#this build a dynamic binary, but with htslib, libBigWig, and libz statically linked in
#remember order is backwards, earliest needed libraries go *last*
target_link_libraries(megadepth_statlib ${CMAKE_SOURCE_DIR}/htslib/libhts.a ${CMAKE_SOURCE_DIR}/libBigWig/libBigWig.a -lz -lcurl -lpthread ${CMAKE_SOURCE_DIR}/libdeflate/libdeflate.a)


##windows build related from here on
#SET(CMAKE_SYSTEM_NAME Windows)
#set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCURL_STATICLIB")
add_executable(megadepth_static.exe megadepth.cpp getline.c)
include_directories(megadepth_static.exe libcurl_windows/include zlib_windows htslib_windows libBigWig_windows)
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

target_link_libraries(megadepth_static.exe -static ${CMAKE_SOURCE_DIR}/htslib_windows/libhts.a ${CMAKE_SOURCE_DIR}/libBigWig_windows/libBigWig.a ${CMAKE_SOURCE_DIR}/libdeflate_windows/libdeflate.a ${CMAKE_SOURCE_DIR}/libcurl_windows/lib/libcurl.a ${CMAKE_SOURCE_DIR}/nghttp2/lib/libnghttp2.a ${CMAKE_SOURCE_DIR}/openssl/lib/libssl.a ${CMAKE_SOURCE_DIR}/openssl/lib/libcrypto.a ${CMAKE_SOURCE_DIR}/libssh2/lib/libssh2.a ${CMAKE_SOURCE_DIR}/brotli/lib/libbrotlidec-static.a ${CMAKE_SOURCE_DIR}/brotli/lib/libbrotlicommon-static.a ${CMAKE_SOURCE_DIR}/zlib_windows/libz.a -lcrypt32 -lwldap32 -pthread -lwsock32 -lws2_32)
