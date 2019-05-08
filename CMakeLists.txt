############################ BASE ######################################

cmake_minimum_required (VERSION 3.13 FATAL_ERROR)

project(PeleC CXX C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for the project
option(ENABLE_DOCUMENTATION "Build documentation" OFF)
option(ENABLE_SPHINX_API_DOCS "Link Doxygen API docs to Sphinx" OFF)
option(ENABLE_ALL_WARNINGS "Show most warnings for most compilers" OFF)

#Enabling tests overrides the executable options
option(ENABLE_TESTS "Enable testing suite" OFF)
option(ENABLE_VERIFICATION "Enable verification suite" OFF)
option(ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)

#Options for the executable in a single build dir
option(PELEC_ENABLE_MPI "Enable MPI" OFF)
#option(PELEC_ENABLE_OPENMP "Enable OpenMP" OFF)
option(PELEC_ENABLE_MASA "Enable MASA for MMS" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

########################### MASA #####################################

if(ENABLE_TESTS)
  set(PELEC_ENABLE_MASA ON)
endif()

if(PELEC_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
  endif()
endif()

########################### PeleC #####################################

if(ENABLE_VERIFICATION AND NOT ENABLE_TESTS)
  message(FATAL_ERROR "-- Testing must be on to enable verification suite")
endif()

if(ENABLE_VERIFICATION)
  message(STATUS "Warning: Verification tests expect a specific Python environment and take a long time to run")
endif()

message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

include(${CMAKE_SOURCE_DIR}/CMake/set_compile_flags.cmake)
set_compile_flags()

find_package(PythonInterp REQUIRED)

include(${CMAKE_SOURCE_DIR}/CMake/build_fcompare.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/build_pelec.cmake)

if(ENABLE_TESTS)
  set(PELEC_ENABLE_MPI ON)
  set(ENABLE_FCOMPARE ON)
  enable_testing()
  include(CTest)
  find_package(MPI REQUIRED)
  build_fcompare()
  add_subdirectory(Testing)
else()
  include(${CMAKE_BINARY_DIR}/exe_options.cmake)
  if(PELEC_ENABLE_MPI)
    find_package(MPI REQUIRED)
  endif()
  #if(PELEC_ENABLE_OPENMP)
  #  find_package(OPENMP REQUIRED)
  #endif()
  build_pelec(PeleC-${PELEC_DIM}D ${CMAKE_BINARY_DIR}/exe_options.cmake)
  if(ENABLE_FCOMPARE)
    build_fcompare()
  endif()
endif()

if(ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
