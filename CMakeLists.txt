cmake_minimum_required(VERSION 3.10)

project(libsais VERSION 2.8.5 LANGUAGES C DESCRIPTION "libsais is a library for linear time suffix array, longest common prefix array and burrows wheeler transform construction based on induced sorting algorithm.")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(LIBSAIS_BUILD_SHARED_LIB "Build libsais as a shared library" OFF)

if(LIBSAIS_BUILD_SHARED_LIB)
    set(LIBSAIS_LIBRARY_TYPE SHARED)
else()
    set(LIBSAIS_LIBRARY_TYPE STATIC)
endif()

add_library(libsais ${LIBSAIS_LIBRARY_TYPE})

target_sources(libsais PRIVATE
    include/libsais16x64.h
    src/libsais16x64.c
)

if(LIBSAIS_BUILD_SHARED_LIB)
    target_compile_definitions(libsais PUBLIC LIBSAIS_SHARED)
    target_compile_definitions(libsais PRIVATE LIBSAIS_EXPORTS)
endif()

target_include_directories(libsais PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# add_executable(libsais64 src/libsais16x64.c src/main.c)
# target_include_directories(libsais64 PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#     $<INSTALL_INTERFACE:include>
# )