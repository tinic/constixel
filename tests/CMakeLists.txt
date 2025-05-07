cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 23)

#set(CMAKE_C_COMPILER "gcc-14")
#set(CMAKE_CXX_COMPILER "g++-14")

if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    message(FATAL_ERROR "Requires at least gcc 12")
endif()

project(constixel_tests)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.gitmodules")
    execute_process(COMMAND git submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

string(FIND ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)

if(NOT ${IS_CLANG} EQUAL -1)
    set(CMAKE_CXX_FLAGS "-Os")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-Os")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "/EHsc")
else()
    message("Unknown compiler!")
endif()

if(NOT ${IS_CLANG} EQUAL -1)
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "-fconstexpr-steps=33554432 -Weverything -Weffc++ -Wall -Wpedantic -Wextra -Wfloat-conversion -Wdouble-promotion -Wuninitialized -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion -Wno-c++98-compat -Wno-pre-c++14-compat -Wno-c++98-compat-local-type-template-args -Wno-c++20-compat -Wno-c++98-c++11-compat-binary-literal")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "-fconstexpr-ops-limit=268435456 -Weffc++ -Wall -Wpedantic -Wextra -Wfloat-conversion -Wdouble-promotion -Wuninitialized -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion -Wno-volatile -Wuseless-cast")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "/W4 /w14242 /w14254 /w1426 /w14265 /w14287 /we4289 /w14296 /w14311 /w14545 /w14546 /w14547 /w14549 /w14555 /w14619 /w14640 /w14826 /w14905 /w14906 /w14928")
else()
    message("Unknown compiler!")
endif()

add_executable(constixel_tests constixel_tests.cpp ${PROJECT_SOURCE_DIR}/../genfonts/fontbm/src/external/lodepng/lodepng.cpp)

target_include_directories(constixel_tests PRIVATE ${PROJECT_SOURCE_DIR}/.. ${PROJECT_SOURCE_DIR}/../genfonts/)
