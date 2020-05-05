#=============================================================================
# PeleC Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

add_executable(${pelec_unit_test_exe_name})

if(PELEC_ENABLE_CUDA)
  get_target_property(UTEST_SOURCES ${pelec_unit_test_exe_name} SOURCES)
  set_source_files_properties(${UTEST_SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
