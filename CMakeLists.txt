cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 23)

project(sixel_test)

string(FIND ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)
if(NOT ${IS_CLANG} EQUAL -1)
    set(CMAKE_CXX_FLAGS "-O3 -fconstexpr-steps=33554432")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-O3 -fconstexpr-ops-limit=268435456")
else()
    message("Unknown compiler!")
endif()

add_executable(sixel_test sixel_test.cpp)
