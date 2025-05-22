cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_STANDARD 23)

set(EXAMPLES "text" "graph" "embedded" "specimen" "minimal" "drawingapi" "industrial" "watchui" "mandelbrot")

project(examples)

string(FIND ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)
foreach(EXAMPLE IN LISTS EXAMPLES)
    add_executable(example_${EXAMPLE} ${EXAMPLE}/example_${EXAMPLE}.cpp)
    target_include_directories(example_${EXAMPLE} PRIVATE ${PROJECT_SOURCE_DIR}/..)
    if(NOT ${IS_CLANG} EQUAL -1)
        target_compile_options(example_${EXAMPLE} PRIVATE -Os -fconstexpr-steps=33554432 -Wall -Wextra -Wpedantic -Weffc++)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(example_${EXAMPLE} PRIVATE -Os -march=native -fconstexpr-ops-limit=268435456 -Wall -Wextra -Wpedantic -Weffc++)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        target_compile_options(example_${EXAMPLE} PRIVATE /EHsc /Ot /arch:AVX2 /constexpr:steps 268435456)
    else()
        message("Unknown compiler!")
    endif()
endforeach()
