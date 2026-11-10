cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_STANDARD 23)

set(EXAMPLES
    "text" "graph" "embedded" "specimen" "minimal" "drawingapi" "industrial"
    "watchui" "mandelbrot" "visualization3d" "shader_fill"
    "energy" "thermostat" "weather" "server_status"
)

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

# Render the visual examples to files without requiring a graphics terminal.
# Kept separate from the default build so source-controlled previews are never
# overwritten as a side effect of compiling an example.
set(PREVIEW_EXAMPLES
    "text" "graph" "drawingapi" "industrial" "watchui" "mandelbrot" "visualization3d" "shader_fill"
    "energy" "thermostat" "weather" "server_status"
)
add_custom_target(example_previews
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/previews"
)
foreach(EXAMPLE IN LISTS PREVIEW_EXAMPLES)
    add_dependencies(example_previews example_${EXAMPLE})
    add_custom_command(TARGET example_previews POST_BUILD
        COMMAND $<TARGET_FILE:example_${EXAMPLE}>
                "${CMAKE_CURRENT_BINARY_DIR}/previews/example_${EXAMPLE}.png"
        VERBATIM
    )
endforeach()
