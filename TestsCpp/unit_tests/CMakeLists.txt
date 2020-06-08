add_executable(${pelec_unit_test_exe_name} "")

target_sources(${pelec_unit_test_exe_name}
  PRIVATE
  unit-tests-main.C
  test-config.C
  )

target_include_directories(${pelec_unit_test_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(AMReXBuildInfo)
generate_buildinfo(${pelec_unit_test_exe_name} ${CMAKE_SOURCE_DIR})
target_include_directories(${pelec_unit_test_exe_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

if(PELEC_ENABLE_CUDA)
  get_target_property(UTEST_SOURCES ${pelec_unit_test_exe_name} SOURCES)
  set_source_files_properties(${UTEST_SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

#Link our unit test executable with GoogleTest
#target_link_libraries(${pelec_unit_test_exe_name} PRIVATE GTest::GTest GTest::Main)
target_link_libraries(${pelec_unit_test_exe_name} PRIVATE gtest)
target_include_directories(${pelec_unit_test_exe_name} SYSTEM PRIVATE
  ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)

target_link_libraries(${pelec_unit_test_exe_name} PRIVATE amrex)

#Define what we want to be installed during a make install 
install(TARGETS ${pelec_unit_test_exe_name}
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
