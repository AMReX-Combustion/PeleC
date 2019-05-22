function(build_fcompare)

  #Expose functions we want to be able to call
  include(${CMAKE_SOURCE_DIR}/CMake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/fcompare_sources.cmake)

  #Clear source file list from any previous executables
  set_property(GLOBAL PROPERTY GlobalFcompareSourceList "") 

  #Aggregate amrex and pelephysics source files
  get_fcompare_sources()
  
  #Put source list from global property into local list
  get_property(FCOMPARE_SOURCES GLOBAL PROPERTY GlobalFcompareSourceList)

  #Create an executable based on all the source files we aggregated
  add_executable(fcompare ${FCOMPARE_SOURCES})

  #Set definitions for our particular executable
  target_compile_definitions(fcompare PRIVATE BL_SPACEDIM=3)
  target_compile_definitions(fcompare PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(fcompare PRIVATE AMREX_SPACEDIM=3)
  target_compile_definitions(fcompare PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(fcompare PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(fcompare PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)
  # CMake BUILD_TYPE should already define this for us
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(fcompare PRIVATE NDEBUG)
  #endif()
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions(fcompare PRIVATE BL_Darwin)
    target_compile_definitions(fcompare PRIVATE AMREX_Darwin)
  endif()
  
  #AMReX include directories
  target_include_directories(fcompare SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(fcompare SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/Plotfile)
  
  #Link our executable to the MPI libraries, etc
  #if(PELEC_ENABLE_MPI)
  #  target_link_libraries(fcompare PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
  #  target_compile_definitions(fcompare PRIVATE BL_USE_MPI)
  #  target_compile_definitions(fcompare PRIVATE AMREX_USE_MPI)
  #endif()
 
  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(fcompare PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/fcompare_fortran_modules")

  #Define what we want to be installed during a make install 
  install(TARGETS fcompare
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
  
endfunction(build_fcompare)
