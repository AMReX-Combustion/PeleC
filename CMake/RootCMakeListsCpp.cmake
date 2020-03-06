############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(PeleC CXX C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for the project
option(PELEC_ENABLE_MASA "Enable MASA for MMS" OFF)
option(PELEC_ENABLE_EB "Enable EB" OFF)
option(PELEC_ENABLE_REACTIONS "Enable reactions" OFF)
option(PELEC_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)
option(PELEC_ENABLE_TESTING "Enable regression tests" OFF)
option(PELEC_ENABLE_VERIFICATION "Enable verification tests" OFF)
option(PELEC_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(PELEC_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)

#Options for the executable
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

#Create target names
set(pelec_unit_test_exe_name "pelec_unit_tests")

if(PELEC_ENABLE_TESTS)
  set(PELEC_ENABLE_MPI ON)
  set(PELEC_ENABLE_MASA ON)
  if(PELEC_TEST_WITH_FCOMPARE)
    set(PELEC_ENABLE_FCOMPARE ON)
  endif()
endif()

if(PELEC_ENABLE_VERIFICATION AND NOT PELEC_ENABLE_TESTS)
  message(FATAL_ERROR "-- Testing must be on to enable verification suite")
endif()

if(PELEC_ENABLE_VERIFICATION)
  message(STATUS "Warning: Verification tests expect a specific Python environment and take a long time to run")
endif()

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "Submodules/AMReX")
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
add_subdirectory(${AMREX_SUBMOD_LOCATION})

if(PELEC_ENABLE_FCOMPARE OR PELEC_TEST_WITH_FCOMPARE)
  add_subdirectory(${AMREX_SUBMOD_LOCATION}/Tools/Plotfile)
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

find_package(Python REQUIRED)
find_package(Threads REQUIRED) # Needed this for the Travis CI system
if(PELEC_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

# General information about machine, compiler, and build type
message(STATUS "PeleC Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

#Create directory unique to executable to store generated files
set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files)
file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

#Generate the AMReX_buildInfo.cpp file with Python
add_custom_target(generate_build_info ALL
   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}/Tools/C_scripts/makebuildinfo_C.py
   --amrex_home "${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"                        
   --COMP ${CMAKE_CXX_COMPILER_ID} --COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
   --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"
   WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
   COMMENT "Generating AMReX_buildInfo.cpp"
)                  
  
#Build pelec executables and link to amrex library
add_subdirectory(ExecCpp)

if(PELEC_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(TestingCpp)
endif()

if(PELEC_ENABLE_DOCUMENTATION)
   add_subdirectory(Docs)
endif()
