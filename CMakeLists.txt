cmake_minimum_required (VERSION 3.14 FATAL_ERROR)

project(PeleC CXX C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

option(PELEC_USE_CPP "Use C++ kernels instead of Fortran kernels" OFF)

if(PELEC_USE_CPP)
  include(${CMAKE_SOURCE_DIR}/CMake/RootCMakeListsCpp.cmake)
else()
  include(${CMAKE_SOURCE_DIR}/CMake/RootCMakeLists.cmake)
endif()
