cmake_minimum_required(VERSION 3.10)

project(unipept-libsais LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(include libsais/include)
set(SRC_FILES src/main.c src/bitpacking.c libsais/src/libsais64.c libsais/src/libsais16x64.c libsais/src/libsais32x64.c)

add_executable(build_ssa ${SRC_FILES})
