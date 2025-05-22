cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 20)

project(genfont)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.gitmodules")
    execute_process(COMMAND git submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_executable(genfont genfont.cpp ${PROJECT_SOURCE_DIR}/fontbm/src/external/lodepng/lodepng.cpp)

