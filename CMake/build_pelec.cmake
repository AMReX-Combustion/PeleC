function(build_pelec pelec_exe_name pelec_exe_options_file)

  unset(PELEC_EXTRA_SOURCES)
  unset(PELEC_DIM)
  unset(PELEC_ENABLE_EB)
  unset(PELEC_ENABLE_MASA)
  unset(PELEC_ENABLE_REACTIONS)
  unset(PELEC_ENABLE_MOL)
  unset(PELEC_ENABLE_PARTICLES)
  unset(PELEC_EOS_MODEL)
  unset(PELEC_REACTIONS_MODEL)
  unset(PELEC_CHEMISTRY_MODEL)
  unset(PELEC_TRANSPORT_MODEL)

  get_filename_component(exe_directory ${pelec_exe_options_file} DIRECTORY)
  include(${pelec_exe_options_file})

  #message("-- PELEC_DIM = ${PELEC_DIM}D")

  if(PELEC_ENABLE_EB)
    set(EB "eb")
  else()
    unset(EB)
  endif()

  #Expose functions we want to be able to call
  include(${CMAKE_SOURCE_DIR}/CMake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/pelephysics_sources.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/pelec_sources.cmake)

  #Clear source file list from any previous executables
  set_property(GLOBAL PROPERTY GlobalSourceList "") 

  #Check for incompatiblities
  if((${PELEC_DIM} GREATER 3) OR (${PELEC_DIM} LESS 1))
    message(FATAL_ERROR "PELEC_DIM must be either 1, 2 or 3.")
  endif()

  if(PELEC_ENABLE_MOL AND NOT PELEC_ENABLE_EB)
    message(FATAL_ERROR "PELEC_ENABLE_MOL does not work without PELEC_ENABLE_EB")
  endif()

  if(${PELEC_DIM} EQUAL 1 AND PELEC_ENABLE_MOL)
    message(FATAL_ERROR "PELEC_ENABLE_MOL does not work with PELEC_DIM=1")
  endif()

  if(${PELEC_DIM} EQUAL 1 AND PELEC_ENABLE_EB)
    message(FATAL_ERROR "PELEC_ENABLE_EB does not work with PELEC_DIM=1")
  endif()

  if("${PELEC_TRANSPORT_MODEL}" STREQUAL "EGLib")
    set(USE_FUEGO ON)
  endif()

  if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
    set(PELEC_TRANSPORT_TYPE "IDEAL_GAS")
  elseif("${PELEC_EOS_MODEL}" STREQUAL "GammaLaw")
    set(PELEC_TRANSPORT_TYPE "IDEAL_GAS")
  else()
    set(PELEC_TRANSPORT_TYPE "REAL_GAS")
  endif()

  #Aggregate amrex and pelephysics source files
  get_pelephysics_sources()
  get_pelec_sources(${pelec_exe_name})
  
  #Put source list from global property into local list
  get_property(PELE_SOURCES GLOBAL PROPERTY GlobalSourceList)

  #Create the full path to the extra case-specific source files
  #Each PELEC_EXTRA_SOURCE must be an explicit path to each source file at the moment
  foreach(PELEC_EXTRA_SOURCE ${PELEC_EXTRA_SOURCES})
    list(APPEND MY_EXTRA_SOURCES ${PELEC_EXTRA_SOURCE})
  endforeach()

  #Create an executable based on all the source files we aggregated
  add_executable(${pelec_exe_name} ${PELE_SOURCES} ${MY_EXTRA_SOURCES})
  target_link_libraries(${pelec_exe_name} PRIVATE amrex${PELEC_DIM}d${EB})

  #AMReX definitions
  target_compile_definitions(${pelec_exe_name} PRIVATE BL_SPACEDIM=${PELEC_DIM})
  target_compile_definitions(${pelec_exe_name} PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(${pelec_exe_name} PRIVATE BL_USE_F_BASELIB)
  target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_SPACEDIM=${PELEC_DIM})
  target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_USE_F_BASELIB)
  target_compile_definitions(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)

  # CMake BUILD_TYPE should already define this
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(${pelec_exe_name} PRIVATE NDEBUG)
  #endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions(${pelec_exe_name} PRIVATE BL_Darwin)
    target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_Darwin)
  endif()

  #PeleC definitions
  if(PELEC_ENABLE_EB)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELE_USE_EB)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_EB)
    target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_USE_EB)
  endif()

  if("${PELEC_TRANSPORT_MODEL}" STREQUAL "EGLib")
    target_compile_definitions(${pelec_exe_name} PRIVATE EGLIB_TRANSPORT)
  elseif("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
    target_compile_definitions(${pelec_exe_name} PRIVATE SIMPLE_TRANSPORT)
  elseif("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
    target_compile_definitions(${pelec_exe_name} PRIVATE CONSTANT_TRANSPORT)
  endif()

  if(PELEC_ENABLE_MOL)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_MOL)
  endif()

  if(PELEC_ENABLE_REACTIONS)
    target_compile_definitions(${pelec_exe_name} PRIVATE REACTIONS)
  endif()

  if(PELEC_ENABLE_PARTICLES)
    target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_PARTICLES)
  endif()
  
  #AMReX include directories
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Amr)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/AmrCore)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/F_Interfaces/AmrCore)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Boundary)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/C_scripts)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/F_Interfaces/Base)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules/amrex${PELEC_DIM}d${EB}_fortran_modules)
  if(PELEC_ENABLE_EB)
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/EB)
  endif()

  #PelePhysics include directories
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport)
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport/${PELEC_TRANSPORT_MODEL})
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Reactions)
  if(("${PELEC_TRANSPORT_MODEL}" STREQUAL "EGLib") OR ("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple"))
    target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Evaluation)
    target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism)
  endif()
  if(DEFINED PELEC_CHEMISTRY_MODEL)
    target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism/Models/${PELEC_CHEMISTRY_MODEL})
  endif()

  #PeleC include directories
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source)
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/param_includes)
  
  #Link our executable to the MPI libraries, etc
  if(PELEC_ENABLE_MPI)
    target_link_libraries(${pelec_exe_name} PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
    target_compile_definitions(${pelec_exe_name} PRIVATE BL_USE_MPI)
    target_compile_definitions(${pelec_exe_name} PRIVATE AMREX_USE_MPI)
  endif()

  #Link our executable to the MASA libraries, etc
  if(PELEC_ENABLE_MASA)
    target_link_libraries(${pelec_exe_name} PRIVATE ${MASA_LIBRARY} ${MASA_FORTRAN_LIBRARY})
    target_compile_definitions(${pelec_exe_name} PRIVATE USE_MASA DO_PROBLEM_POST_TIMESTEP DO_PROBLEM_POST_INIT)
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_INCLUDE_DIRS})
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_MOD_DIRS})
  endif()

  #if(PELEC_ENABLE_OPENMP)
  #  string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")
  #  string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
  #  string(APPEND CMAKE_Fortran_FLAGS " ${OpenMP_Fortran_FLAGS}")
  #endif()

  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(${pelec_exe_name} PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/${pelec_exe_name}_fortran_modules")

  #Create directory unique to executable to store generated files
  set(GENERATED_FILES_DIR ${CMAKE_BINARY_DIR}/generated_files/${pelec_exe_name}_generated_files)
  file(MAKE_DIRECTORY ${GENERATED_FILES_DIR})

  set(PARAMETER_DIRS "")
  string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Eos/_parameters")
  string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Reactions/_parameters")
  string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport/_parameters")
  string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Eos/${PELEC_EOS_MODEL}/_parameters")
  string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport/${PELEC_TRANSPORT_MODEL}/_parameters")
  if("${PELEC_REACTIONS_MODEL}" STREQUAL "Fuego")
    string(APPEND PARAMETER_DIRS " ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Reactions/${PELEC_REACTIONS_MODEL}/_parameters")
  endif()

  if(PYTHON_FOUND)
     #Generate the extern.f90 file with Python
     add_custom_target(generate_extern_${pelec_exe_name} ALL
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/F_scripts/write_probin.py
        -t ${CMAKE_SOURCE_DIR}/Source/extern_probin.template
        -o extern.f90 -n extern --pa "${PARAMETER_DIRS}"
        WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/extern.f90
        COMMENT "Generating extern.f90"
     )

     #Generate the AMReX_buildInfo.cpp file with Python
     add_custom_target(generate_build_info_${pelec_exe_name} ALL
        COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/C_scripts/makebuildinfo_C.py
        --amrex_home "${CMAKE_SOURCE_DIR}/Submodules/AMReX"                        
        --COMP ${CMAKE_C_COMPILER_ID} --COMP_VERSION ${CMAKE_C_COMPILER_VERSION}
        --FCOMP ${CMAKE_Fortran_COMPILER_ID} --FCOMP_VERSION ${CMAKE_Fortran_COMPILER_VERSION}
        --GIT "${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/Submodules/AMReX ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics"
        WORKING_DIRECTORY ${GENERATED_FILES_DIR} BYPRODUCTS ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
        COMMENT "Generating AMReX_buildInfo.cpp"
     )                  
  endif()
  
  #Set the dependencies on targets so the generated source code files are there before we try to build the executable 
  add_dependencies(${pelec_exe_name} generate_extern_${pelec_exe_name} generate_build_info_${pelec_exe_name})
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${pelec_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
  
endfunction(build_pelec pelec_exe_name pelec_exe_options_file)
