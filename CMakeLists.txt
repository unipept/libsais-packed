cmake_minimum_required(VERSION 3.10)

project(libsais-packed LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(include libsais/include)
set(SRC_FILES src/main.c src/bitpacking.c libsais/src/libsais64.c libsais/src/libsais16x64.c libsais/src/libsais32x64.c)

add_executable(libsais-packed ${SRC_FILES})

add_library(libsais STATIC)

target_sources(libsais PRIVATE
    libsais/include/libsais64.h
    libsais/include/libsais16x64.h
    libsais/include/libsais32x64.h
    libsais/src/libsais64.c
    libsais/src/libsais16x64.c
    libsais/src/libsais32x64.c
)
