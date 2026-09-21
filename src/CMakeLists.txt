cmake_minimum_required(VERSION 3.16)
project(ocp-framework LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ADIOS2 REQUIRED)

add_executable(ocp-framework ocp-framework.cpp)
target_link_libraries(ocp-framework PRIVATE adios2::cxx)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_executable(write-test-dump tests/write-dump.cpp)
  target_link_libraries(write-test-dump PRIVATE adios2::cxx)
  add_test(NAME observables
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check-observables.py
      $<TARGET_FILE:ocp-framework> $<TARGET_FILE:write-test-dump>)
endif()
