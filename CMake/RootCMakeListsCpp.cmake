############################ BASE ######################################

cmake_minimum_required (VERSION 3.14 FATAL_ERROR)
project(PeleC CXX C Fortran)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)

########################## OPTIONS #####################################

#General options for the project
option(PELEC_ENABLE_MASA "Enable MASA for MMS" OFF)
option(PELEC_ENABLE_EB "Enable EB" OFF)
option(PELEC_ENABLE_ALL_WARNINGS "Enable all compiler warnings" OFF)

#Options for the executable
option(PELEC_ENABLE_MPI "Enable MPI" OFF)
option(PELEC_ENABLE_OPENMP "Enable OpenMP" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Create target names
set(pelec_exe_name "pelec")
set(pelec_unit_test_exe_name "${pelec_exe_name}_unit_tests")

#Create main target executable
add_executable(${pelec_exe_name} "")

########################### AMReX #####################################

set(AMREX_SUBMOD_LOCATION "Submodules/AMReX")
include(${CMAKE_SOURCE_DIR}/CMake/SetAmrexOptions.cmake)
add_subdirectory(${AMREX_SUBMOD_LOCATION})

########################### MASA #####################################

if(PELEC_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA QUIET REQUIRED)
  if(MASA_FOUND)
    message(STATUS "Found MASA = ${MASA_DIR}")
    #Link our executable to the MASA libraries, etc
    target_link_libraries(${pelec_exe_name} PRIVATE ${MASA_LIBRARY})
    target_compile_definitions(${pelec_exe_name} PRIVATE USE_MASA DO_PROBLEM_POST_TIMESTEP DO_PROBLEM_POST_INIT)
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_INCLUDE_DIRS})
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_MOD_DIRS})
  endif()
endif()

########################### PeleC #####################################

find_package(Python REQUIRED)
find_package(Threads REQUIRED) # Needed this for the Travis CI system
if(PELEC_ENABLE_MPI)
  find_package(MPI REQUIRED)
  target_link_libraries(${pelec_exe_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
endif()

include(${CMAKE_SOURCE_DIR}/CMake/SetCompileFlags.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

# General information about machine, compiler, and build type
message(STATUS "PeleC Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

#Create directory unique to executable to store generated files
set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files/${pelec_exe_name}_generated_files)
file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

#Generate the AMReX_buildInfo.cpp file with Python
add_custom_target(generate_build_info_${pelec_exe_name} ALL
   COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}/Tools/C_scripts/makebuildinfo_C.py
   --amrex_home "${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"                        
   --COMP ${CMAKE_CXX_COMPILER_ID} --COMP_VERSION ${CMAKE_CXX_COMPILER_VERSION}
   --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/${AMREX_SUBMOD_LOCATION}"
   WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
   COMMENT "Generating AMReX_buildInfo.cpp"
)                  
  
#Set the dependencies on targets so the generated source code files are there before we try to build the executable 
add_dependencies(${pelec_exe_name} generate_build_info_${pelec_exe_name})

#Build pelec and link to amrex library
add_subdirectory(SourceCpp)

#Define what we want to be installed during a make install 
install(TARGETS ${pelec_exe_name}
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
