function(build_unit_test utest_exe_name utest_exe_options_file)

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

  get_filename_component(exe_directory ${utest_exe_options_file} DIRECTORY)
  include(${utest_exe_options_file})

  #Expose functions we would like to be able to call here
  include(${CMAKE_SOURCE_DIR}/CMake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/pelec_unit_test_sources.cmake)

  #Clear out list of source files for unit tests
  set_property(GLOBAL PROPERTY GlobalUnitSourceList "")

  #Aggregate source files for unit tests
  get_pelec_unit_test_sources()

  #Put source list from global property into local list
  get_property(PELE_UNIT_SOURCES GLOBAL PROPERTY GlobalUnitSourceList)

  #Create our unit test executable
  add_executable(${utest_exe_name} ${PELE_UNIT_SOURCES})

  #Link our unit test executable with Google Test, etc
  #target_link_libraries(${utest_exe_name} PRIVATE GTest::GTest GTest::Main)
  target_link_libraries(${utest_exe_name} PRIVATE gtest gtest_main)
  target_include_directories(${utest_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)

  #Set definitions for our particular executable
  target_compile_definitions(${utest_exe_name} PRIVATE BL_SPACEDIM=${PELEC_DIM})
  target_compile_definitions(${utest_exe_name} PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(${utest_exe_name} PRIVATE AMREX_SPACEDIM=${PELEC_DIM})
  target_compile_definitions(${utest_exe_name} PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(${utest_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(${utest_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)
  # CMake BUILD_TYPE should already define this for us
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(${utest_exe_name} PRIVATE NDEBUG)
  #endif()
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions(${utest_exe_name} PRIVATE BL_Darwin)
    target_compile_definitions(${utest_exe_name} PRIVATE AMREX_Darwin)
  endif()
  
  #Link our unit test executable with MPI, etc
  if(PELEC_ENABLE_MPI)
    target_link_libraries(${utest_exe_name} PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
    target_compile_definitions(${utest_exe_name} PRIVATE BL_USE_MPI)
    target_compile_definitions(${utest_exe_name} PRIVATE AMREX_USE_MPI)
  endif()

  #Choose what we want installed if we do a make install
  install(TARGETS ${utest_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_unit_test utest_exe_name)
