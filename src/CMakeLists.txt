cmake_minimum_required(VERSION 3.16)
project(analyze_lammps_adios LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ADIOS2 REQUIRED)

add_executable(analyze_lammps_adios analyze_lammps_adios.cpp)
target_link_libraries(analyze_lammps_adios PRIVATE adios2::cxx)
