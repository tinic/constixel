cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 23)

#set(CMAKE_C_COMPILER "gcc-14")
#set(CMAKE_CXX_COMPILER "g++-14")

if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12)
    message(FATAL_ERROR "Requires at least gcc 12")
endif()

project(constixel_tests)

# ----- UNIT TESTS -------------------------------------------------------------------------

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../.gitmodules")
    execute_process(COMMAND git submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

string(FIND ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)

if(NOT ${IS_CLANG} EQUAL -1)
    set(CMAKE_CXX_FLAGS "-Os -fbracket-depth=1024")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-Os")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "/EHsc")
else()
    message("Unknown compiler!")
endif()

if(NOT ${IS_CLANG} EQUAL -1)
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "-Os -fno-rtti -fno-exceptions -fno-threadsafe-statics -fconstexpr-steps=33554432 -Weverything -Weffc++ -Wall -Wpedantic -Wextra -Wfloat-conversion -Wdouble-promotion -Wuninitialized -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion -Wno-c++98-compat -Wno-pre-c++14-compat -Wno-c++98-compat-local-type-template-args -Wno-c++20-compat -Wno-c++98-c++11-compat-binary-literal -Wno-unsafe-buffer-usage")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "-Os -march=native -fno-rtti -fno-exceptions -fno-threadsafe-statics -fconstexpr-ops-limit=268435456 -Weffc++ -Wall -Wpedantic -Wextra -Wfloat-conversion -Wdouble-promotion -Wuninitialized -Wshadow -Wnon-virtual-dtor -Wunused -Woverloaded-virtual -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wstrict-aliasing -Wcast-align -Wconversion -Wno-volatile -Wuseless-cast")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set_source_files_properties(constixel_tests.cpp PROPERTIES COMPILE_FLAGS "/EHsc /Ot /arch:AVX2 /permissive- /constexpr:steps 268435456 /Wall /external:anglebrackets /external:W0 /wd4710 /wd4711 /wd4514 /wd5045 /wd4820 /wd4324")
else()
    message("Unknown compiler!")
endif()

add_executable(constixel_tests constixel_tests.cpp ${PROJECT_SOURCE_DIR}/../genfonts/fontbm/src/external/lodepng/lodepng.cpp)

target_compile_definitions(constixel_tests PRIVATE CMAKE_PROJECT_PATH="${PROJECT_SOURCE_DIR}")

target_include_directories(constixel_tests PRIVATE ${PROJECT_SOURCE_DIR}/.. ${PROJECT_SOURCE_DIR}/../genfonts/)

# ----- FONTS -------------------------------------------------------------------------

set(HEADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../fonts")   # adjust as needed
file(GLOB HDRS "${HEADER_DIR}/*.h")                     # or *.hpp

# --- split into _aa and _mono ------------------------------------------------
set(HDRS_AA   "")
set(HDRS_MONO "")
foreach(hdr IN LISTS HDRS)
    get_filename_component(name "${hdr}" NAME)            # e.g. foo_aa.h
    if(name MATCHES "_aa\\.h$")
        list(APPEND HDRS_AA   "${hdr}")
    elseif(name MATCHES "_mono\\.h$")
        list(APPEND HDRS_MONO "${hdr}")
    endif()
endforeach()

# --- generate header ---------------------------------------------------------
set(GENERATED_HDR "${CMAKE_BINARY_DIR}/fontslist.h")

set(_includes_aa   "")
set(_includes_mono "")
set(_names_aa      "")
set(_names_mono    "")

foreach(hdr IN LISTS HDRS_AA)
    get_filename_component(name    "${hdr}" NAME)
    get_filename_component(name_we "${hdr}" NAME_WE)
    string(APPEND _includes_aa "#include \"${name}\"\n")
    list(APPEND _names_aa "${name_we}")
endforeach()

foreach(hdr IN LISTS HDRS_MONO)
    get_filename_component(name    "${hdr}" NAME)
    get_filename_component(name_we "${hdr}" NAME_WE)
    string(APPEND _includes_mono "#include \"${name}\"\n")
    list(APPEND _names_mono "${name_we}")
endforeach()

string(JOIN ", " _names_aa_joined   ${_names_aa})
string(JOIN ", " _names_mono_joined ${_names_mono})

string(APPEND _includes_aa   "\n#define AA_HEADERS   ${_names_aa_joined}\n")
string(APPEND _includes_mono "\n#define MONO_HEADERS ${_names_mono_joined}\n")

file(WRITE "${GENERATED_HDR}" "${_includes_aa}${_includes_mono}")

# --- make the header usable --------------------------------------------------


add_executable(constixel_fonts constixel_fonts.cpp)
target_include_directories(constixel_fonts PRIVATE
    ${CMAKE_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/..
    ${PROJECT_SOURCE_DIR}/../fonts
)
