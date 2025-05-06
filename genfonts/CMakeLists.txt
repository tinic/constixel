cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 23)

project(genfont)

add_executable(genfont genfont.cpp ${PROJECT_SOURCE_DIR}/fontbm/src/external/lodepng/lodepng.cpp)

