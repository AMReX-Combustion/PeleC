# Include files with functions we would like to call
include(${CMAKE_SOURCE_DIR}/CMake/build_unit_test.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/build_pelec.cmake)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

message("-- Test golds directory: ${CMAKE_CURRENT_SOURCE_DIR}/PeleCGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION}")

# Standard regression test
function(add_test_r TEST_NAME NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${CMAKE_CURRENT_SOURCE_DIR}/PeleCGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION}/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Get test options
    set(EXE_OPTIONS_FILE ${CURRENT_TEST_SOURCE_DIR}/exe_options.cmake)
    # Define our test options
    include(${EXE_OPTIONS_FILE})
    # Find fcompare
    set(FCOMPARE ${CMAKE_BINARY_DIR}/fcompare)
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "max_step=10 amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.probin_file=${TEST_NAME}.probin")
    # Build the exe for the test
    build_pelec(PeleC-${TEST_NAME} ${EXE_OPTIONS_FILE})
    # Place the exe in the correct working directory
    set_target_properties(PeleC-${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/")
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} ${CURRENT_TEST_BINARY_DIR}/PeleC-${TEST_NAME} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} && ${FCOMPARE} ${PLOT_GOLD} ${PLOT_TEST}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression")
endfunction(add_test_r)

# Standard unit test
function(add_test_u TEST_NAME NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Build exe for test
    build_unit_test(${TEST_NAME} ${CURRENT_TEST_SOURCE_DIR}/exe_options.cmake)
    # Place the exe in the correct working directory
    set_target_properties(${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/")
    # Add test and commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(fiab-2d 4)
add_test_r(fiab-3d 4)
add_test_r(hit-3d-1 4)
add_test_r(hit-3d-2 4)
add_test_r(hit-3d-3 4)
# add_test_r(mms-1d-1 4)
# add_test_r(mms-2d-1 4)
# add_test_r(mms-2d-2 4)
# add_test_r(mms-3d-1 4)
# add_test_r(mms-3d-2 4)
# add_test_r(mms-3d-3 4)
# add_test_r(mms-3d-4 1)
add_test_r(sod-3d-1 4)
add_test_r(tg-2d-1 4)
add_test_r(tg-3d-1 4)
add_test_r(tg-3d-2 4)

#=============================================================================
# Verification tests
#=============================================================================

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit-tests-3d 1)

#=============================================================================
# Performance tests
#=============================================================================
