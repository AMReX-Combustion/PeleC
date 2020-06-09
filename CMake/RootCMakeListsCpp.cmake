############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(PeleC CXX C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for the project
option(PELEC_ENABLE_DOCUMENTATION "Build documentation" OFF)
option(PELEC_ENABLE_EB "Enable EB" OFF)
option(PELEC_ENABLE_REACTIONS "Enable reactions" ON)
option(PELEC_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)
option(PELEC_ENABLE_TESTS "Enable regression and unit tests" OFF)
option(PELEC_ENABLE_VERIFICATION_TESTS "Enable verification tests" OFF)
option(PELEC_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(PELEC_ENABLE_FCOMPARE_FOR_TESTS "Check test plots against gold files" OFF)

#Options for performance
option(PELEC_ENABLE_MPI "Enable MPI" OFF)
option(PELEC_ENABLE_OPENMP "Enable OpenMP" OFF)
option(PELEC_ENABLE_CUDA "Enable CUDA" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(PELEC_ENABLE_REACTIONS)
  set(PELEC_ENABLE_EXPLICIT_REACT ON)
endif()

if(PELEC_ENABLE_VERIFICATION_TESTS)
  set(PELEC_ENABLE_TESTS ON)
  message(STATUS "Warning: Verification tests expect a specific Python environment and take a long time to run")
endif()

if(PELEC_ENABLE_TESTS AND PELEC_ENABLE_FCOMPARE_FOR_TESTS)
  set(PELEC_ENABLE_FCOMPARE ON)
endif()

if(PELEC_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "9.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 9.0 or newer.")
  endif()
endif()

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})

########################### MASA #####################################

if(PELEC_ENABLE_TESTS)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
  endif()
endif()

########################### PeleC #####################################

if(PELEC_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "PeleC Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

#Build pelec executables and link to amrex library
add_subdirectory(ExecCpp)

if(PELEC_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory("Submodules/GoogleTest")
  add_subdirectory(TestsCpp)
endif()

if(PELEC_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
