cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 23)

if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    message(FATAL_ERROR "Requires at least gcc 12")
endif()

project(constixel_tests)

string(FIND ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)
if(NOT ${IS_CLANG} EQUAL -1)
    set(CMAKE_CXX_FLAGS "-O2 -fconstexpr-steps=33554432 -Weffc++ -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-O2 -fconstexpr-ops-limit=268435456 -Weffc++ -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion -Wno-volatile")
else()
    message("Unknown compiler!")
endif()

add_executable(constixel_tests constixel_tests.cpp ${PROJECT_SOURCE_DIR}/../fontbm/src/external/lodepng/lodepng.cpp)
