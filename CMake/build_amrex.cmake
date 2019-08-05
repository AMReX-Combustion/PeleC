function(build_amrex_library AMREX_DIM AMREX_ENABLE_EB)
  # Set library suffixes for EB enabled
  if(AMREX_ENABLE_EB)
    set(EB "eb")
  else()
    unset(EB)
  endif()

  #Expose functions we want to be able to call
  include(${CMAKE_SOURCE_DIR}/CMake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/amrex_sources.cmake)

  #Clear source file list from any previous executables
  set_property(GLOBAL PROPERTY GlobalSourceList "")

  #Aggregate amrex and pelephysics source files
  get_amrex_sources()

  #Put source list from global property into local list
  get_property(AMREX_SOURCES GLOBAL PROPERTY GlobalSourceList)

  #Create an executable based on all the source files we aggregated
  add_library(amrex${AMREX_DIM}d${EB} ${AMREX_SOURCES})

  #AMReX definitions
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE BL_SPACEDIM=${AMREX_DIM})
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE BL_USE_F_BASELIB)
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_SPACEDIM=${AMREX_DIM})
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_USE_F_BASELIB)
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)

  # CMake BUILD_TYPE should already define this for us
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE NDEBUG)
  #endif()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE BL_Darwin)
    target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_Darwin)
  endif()

  #PeleC definitions
  if(AMREX_ENABLE_EB)
    target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_USE_EB)
  endif()

  #AMReX include directories
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Amr)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/AmrCore)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/F_Interfaces/AmrCore)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Boundary)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/C_scripts)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/F_Interfaces/Base)
  target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_BINARY_DIR}/fortran_modules/amrex${AMREX_DIM}d${EB}_fortran_modules)
  if(AMREX_ENABLE_EB)
    target_include_directories(amrex${AMREX_DIM}d${EB} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/EB)
  endif()

  #Link our executable to the MPI libraries, etc
  if(PELEC_ENABLE_MPI)
    target_link_libraries(amrex${AMREX_DIM}d${EB} PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
    target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE BL_USE_MPI)
    target_compile_definitions(amrex${AMREX_DIM}d${EB} PRIVATE AMREX_USE_MPI)
  endif()

  target_link_libraries(amrex${AMREX_DIM}d${EB} PRIVATE Threads::Threads)

  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(amrex${AMREX_DIM}d${EB} PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/amrex${AMREX_DIM}d${EB}_fortran_modules")

  #Define what we want to be installed during a make install 
  install(TARGETS amrex${AMREX_DIM}d${EB}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
endfunction(build_amrex_library AMREX_DIM AMREX_ENABLE_EB)

function(build_amrex)
  if(ENABLE_TESTS)
    # Build all libraries if testing is enabled
    foreach(AMREX_ENABLE_EB IN ITEMS TRUE FALSE)
      foreach(AMREX_DIM IN ITEMS 1 2 3)
        if(${AMREX_DIM} EQUAL 1 AND AMREX_ENABLE_EB)
          continue()
        endif()
        build_amrex_library(${AMREX_DIM} ${AMREX_ENABLE_EB})
      endforeach()
    endforeach()
  else()
    # Otherwise only build the necessary library for the exe
    set(AMREX_ENABLE_EB ${PELEC_ENABLE_EB})
    set(AMREX_DIM ${PELEC_DIM})
    build_amrex_library(${AMREX_DIM} ${AMREX_ENABLE_EB})
  endif()
endfunction(build_amrex)
