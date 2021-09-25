cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# set the project name and version
project(AlgRemez VERSION 1.0.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# add a directory for cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(ExternalProject)

# Link OpempMP CMake automatically adds 
# _OPENMP to the CXX flags
option(ENABLE_OPENMP "enable OpenMP" ON)
if(ENABLE_OPENMP)
  find_package(OpenMP)
endif()

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

add_subdirectory(src)
add_subdirectory(test)
