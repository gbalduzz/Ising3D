cmake_minimum_required(VERSION 3.0)
project(Ising3D)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)
add_library(src src/ising_lattice.cpp)

add_executable(ising main_ising.cpp)
target_link_libraries(ising src)
