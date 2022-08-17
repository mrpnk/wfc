cmake_minimum_required(VERSION 3.23)
project(WFC_2D)

set(CMAKE_CXX_STANDARD 23)
set(CXX_STANDARD 23)

file(GLOB cppSources src/*.cpp)
add_executable(WFC_2D ${cppSources})
