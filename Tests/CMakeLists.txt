#=============================================================================
# PeleC Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

#Google Test Stuff (GoogleTest submodule must be checked out in this directory)
add_subdirectory(GoogleTest)
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)
set_target_properties(gtest PROPERTIES FOLDER Submodules)
set_target_properties(gtest_main PROPERTIES FOLDER Submodules)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
