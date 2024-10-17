cmake_minimum_required(VERSION 3.22)
project(FEM_Plasticity C)

# Set the C Standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_DIR}/lib)

option(BUILD_WITH_SOLVER "Build with solver module" OFF)
option(BUILD_WITH_MATH "Build with math module" OFF)
option(BUILD_WITH_UTILS "Build with utils module" OFF)
option(BUILD_ALL "Build all the modules" ON)

add_subdirectory(src)
add_subdirectory(test)
