#=============================================================================
# PeleC Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
