cmake_minimum_required(VERSION 3.16)
project(vdpsu-bench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# VDPSU_NATIVE128: select the 128-bit OpenFHE build for full statistical security.
# When ON, finds OpenFHE at /Users/msunkim/opt/openfhe-128 (built with NATIVE_SIZE=128)
# and compiles with -DVDPSU_NATIVE128 so params.h picks the 2^62 plaintext modulus.
option(VDPSU_NATIVE128 "Use the 128-bit NATIVEINT OpenFHE install" OFF)

if(VDPSU_NATIVE128)
    message(STATUS "VDPSU: using NATIVEINT=128 OpenFHE at /Users/msunkim/opt/openfhe-128")
    set(OpenFHE_DIR /Users/msunkim/opt/openfhe-128/lib/OpenFHE)
    add_definitions(-DVDPSU_NATIVE128)
endif()

# Find OpenFHE
find_package(OpenFHE REQUIRED)

if(VDPSU_NATIVE128)
    link_directories(/Users/msunkim/opt/openfhe-128/lib)
endif()

# Find NTL
find_library(NTL_LIB ntl PATHS /usr/local/lib /opt/homebrew/lib)
find_library(GMP_LIB gmp PATHS /usr/local/lib /opt/homebrew/lib)

# Find Crypto++
find_library(CRYPTOPP_LIB cryptopp PATHS /usr/local/lib /opt/homebrew/lib)

# OpenMP (required by OpenFHE)
link_directories(/usr/local/opt/libomp/lib /opt/homebrew/opt/libomp/lib)

# Enable OpenMP for our code too
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include -Wno-unused-command-line-argument")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/third-party/include
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
    /usr/local/include
    /opt/homebrew/include
)

# Library sources
set(LIB_SOURCES
    src/prf.cc
    src/mht.cc
    src/phf.cc
    src/bloom_filter.cc
    src/fhe_utils.cc
    src/client.cc
    src/cloud.cc
)

# Build library
add_library(vdpsu STATIC ${LIB_SOURCES})
target_link_libraries(vdpsu
    ${OpenFHE_LIBRARIES}
    ${NTL_LIB}
    ${GMP_LIB}
    ${CRYPTOPP_LIB}
)

# Benchmark executable
add_executable(bench bench/main.cc)
target_link_libraries(bench vdpsu)

# Test executable
add_executable(vdpsu_tests test/test_main.cc)
target_link_libraries(vdpsu_tests vdpsu)

# Parameter comparison test
add_executable(test_params test/test_params.cc)
target_link_libraries(test_params vdpsu)

# Large plaintext modulus test
add_executable(test_large_p test/test_large_p.cc)
target_link_libraries(test_large_p vdpsu)
