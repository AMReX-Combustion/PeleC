############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(PeleC CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)
include(PeleCUtils)

########################## OPTIONS #####################################

#General options for all executables in the project
set(PELEC_DIM "3" CACHE STRING "Number of physical dimensions")
option(PELEC_ENABLE_AMREX_PARTICLES "Enable particles in AMReX library" ON)
option(PELEC_ENABLE_CLANG_TIDY "Enable clang-tidy analysis" OFF)
option(PELEC_ENABLE_CPPCHECK "Enable cppcheck analysis" OFF)
option(PELEC_ENABLE_MASA "Enable tests that require MASA" OFF)
option(PELEC_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(PELEC_ENABLE_FCOMPARE_FOR_TESTS "Check test plots against gold files" OFF)
option(PELEC_SAVE_GOLDS "Enable saving of gold files to a specified directory" OFF)
option(PELEC_ENABLE_SANITIZE_FOR_TESTS "Currently only disables certain long running MMS tests if set" OFF)
option(PELEC_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)
option(PELEC_ENABLE_TINY_PROFILE "Enable tiny profiler in AMReX" OFF)
option(PELEC_ENABLE_HDF5 "Enable plot file output using HDF5" OFF)
option(PELEC_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5" OFF)
option(PELEC_ENABLE_ASCENT "Enable Ascent in-situ visualization" OFF)
set(PELEC_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

#Options for performance
option(PELEC_ENABLE_MPI "Enable MPI" OFF)
option(PELEC_ENABLE_OPENMP "Enable OpenMP" OFF)
option(PELEC_ENABLE_CUDA "Enable CUDA" OFF)
option(PELEC_ENABLE_HIP "Enable HIP" OFF)
option(PELEC_ENABLE_SYCL "Enable SyCL" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(PELEC_ENABLE_FCOMPARE_FOR_TESTS OR PELEC_ENABLE_AMREX_PARTICLES)
  set(PELEC_ENABLE_FCOMPARE ON)
endif()

if(PELEC_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  if(PELEC_ENABLE_MASA)
    message(FATAL_ERROR "MASA enabled with CUDA is unsupported.")
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70)
  endif()
endif()

if(PELEC_DIM EQUAL 1)
  message(FATAL_ERROR "PeleC does not enable 1D.")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(PELEC_ENABLE_FPE_TRAP_FOR_TESTS OFF)
  message(WARNING "Disabling FPE trapping for tests when using AppleClang.")
endif()

if(PELEC_ENABLE_HDF5_ZFP AND (NOT PELEC_ENABLE_HDF5))
   message(FATAL_ERROR "HDF5 must be enabled to use ZFP compression in HDF5 plot files")
endif()

########################### AMReX #####################################

message(STATUS "AMReX Configure Section")
set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX")
include(SetAmrexOptions)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(SetAmrexCompileFlags)

########################### AMReX-Hydro ###############################

message(STATUS "AMReX-Hydro Configure Section")
set(AMREX_HYDRO_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX-Hydro")
include(SetAmrexHydroOptions)
add_subdirectory(${AMREX_HYDRO_SUBMOD_LOCATION})

########################## SUNDIALS ###################################

message(STATUS "SUNDIALS Configure Section")
set(SUNDIALS_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/sundials")
include(SetSundialsOptions)
#BUILD_TESTING is an old CMake keyword so don't clear it in sundials configure
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
add_subdirectory(${SUNDIALS_SUBMOD_LOCATION})

########################### MASA #####################################

if(PELEC_ENABLE_MASA)
  message(STATUS "MASA Configure Section")
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
  endif()
endif()

########################### PeleC #####################################

if(PELEC_ENABLE_MPI)
  message(STATUS "MPI Configure Section")
  find_package(MPI REQUIRED CXX)
endif()

# General information about machine, compiler, and build type
message(STATUS "PeleC Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

message(STATUS "PELEC_PRECISION = ${PELEC_PRECISION}")
if("${PELEC_PRECISION}" STREQUAL "SINGLE" AND PELEC_ENABLE_MASA)
  message(FATAL_ERROR "Single precision currently not allowed in AMReX with MASA enabled")
endif()

init_code_checks()

include(SetRpath)

add_subdirectory(Exec)

enable_testing()
include(CTest)
#add_subdirectory("Submodules/GoogleTest")
add_subdirectory(Tests)
