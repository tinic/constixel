cmake_minimum_required(VERSION 3.21)

if (PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_STANDARD 23)
endif()

project(constixel LANGUAGES)

add_library(constixel INTERFACE)

target_include_directories(constixel INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(constixel INTERFACE cxx_std_20)

install(TARGETS constixel
        EXPORT ConstixelConfig)

install(FILES constixel.hpp
        DESTINATION include)

install(EXPORT ConstixelConfig
        NAMESPACE constixel::
        DESTINATION lib/cmake/constixel)

add_library(constixel::constixel ALIAS constixel)

if (PROJECT_IS_TOP_LEVEL)
    enable_testing()
    add_subdirectory(tests)
    add_subdirectory(examples)
endif()
