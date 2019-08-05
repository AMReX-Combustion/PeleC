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
option(ENABLE_FEXTREMA "Enable building fextrema when not testing" OFF)
option(TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(TEST_WITH_FEXTREMA "Check test plots against maxima and minima files" OFF)

#Options for the executable in a single build dir
option(PELEC_ENABLE_MPI "Enable MPI" OFF)
#option(PELEC_ENABLE_OPENMP "Enable OpenMP" OFF)
option(PELEC_ENABLE_MASA "Enable MASA for MMS" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Options for libraries we build
#SET(BUILD_SHARED_LIBS OFF)

if(ENABLE_TESTS)
  set(PELEC_ENABLE_MPI ON)
  set(PELEC_ENABLE_MASA ON)
  if(TEST_WITH_FCOMPARE)
    set(ENABLE_FCOMPARE ON)
  endif()
  if(TEST_WITH_FEXTREMA)
    set(ENABLE_FEXTREMA ON)
  endif()
endif()

########################### MASA #####################################

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

include(${CMAKE_SOURCE_DIR}/CMake/set_compile_flags.cmake)
set_compile_flags()

find_package(Python REQUIRED)
find_package(Threads REQUIRED) # Needed this for the Travis CI system
if(PELEC_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

include(${CMAKE_SOURCE_DIR}/CMake/build_amrex.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/build_pelec.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/build_plot_tool.cmake)

# General information about machine, compiler, and build type
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Regular flags we have added
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

# User can set CMAKE_BUILD_TYPE case insensitively but we want it uppercase
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE)

# Build type flags in which CMake adds for us
message(STATUS "CMAKE_CXX_FLAGS_${BUILD_TYPE} = ${CMAKE_CXX_FLAGS_${BUILD_TYPE}}")
message(STATUS "CMAKE_C_FLAGS_${BUILD_TYPE} = ${CMAKE_C_FLAGS_${BUILD_TYPE}}")
message(STATUS "CMAKE_Fortran_FLAGS_${BUILD_TYPE} = ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")

# Use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# When building, don't use the install RPATH already (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

if(ENABLE_FCOMPARE)
  build_plot_tool(fcompare)
endif()
if(ENABLE_FEXTREMA)
  build_plot_tool(fextrema)
endif()

if(ENABLE_TESTS)
  enable_testing()
  include(CTest)
  build_amrex()
  add_subdirectory(Testing)
else()
  include(${CMAKE_BINARY_DIR}/exe_options.cmake)
  build_amrex()
  build_pelec(PeleC-${PELEC_DIM}D ${CMAKE_BINARY_DIR}/exe_options.cmake)
endif()

if(ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
