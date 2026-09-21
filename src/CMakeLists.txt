cmake_minimum_required(VERSION 3.16)
project(ocp-framework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ADIOS2 REQUIRED)

add_executable(ocp-framework ocp-framework.cpp)
target_link_libraries(ocp-framework PRIVATE adios2::cxx)
