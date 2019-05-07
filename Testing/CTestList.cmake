# Include files with functions we would like to call
include(${CMAKE_SOURCE_DIR}/CMake/build_unit_test.cmake)
include(${CMAKE_SOURCE_DIR}/CMake/build_pelec.cmake)

# Set location of gold files according to system/compiler/compiler_version
set(GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PeleCGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})

if(TEST_WITH_FCOMPARE)
  message("-- Test golds directory: ${GOLD_FILES_DIRECTORY}")
endif()

# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r TEST_NAME NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
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
    # Either just run the tests, or also use fcompare to test diffs in plots against gold files
    if(TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${FCOMPARE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
    # Place the exe in the correct working directory
    set_target_properties(PeleC-${TEST_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/")
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} ${CURRENT_TEST_BINARY_DIR}/PeleC-${TEST_NAME} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression")
endfunction(add_test_r)

# Verification test with 1 resolution
function(add_test_v1 TEST_NAME TEST_DEPENDENCY NP)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(TEST_DEPENDENCY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_DEPENDENCY})
    set(TEST_DEPENDENCY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_DEPENDENCY})
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Get test options
    set(EXE_OPTIONS_FILE ${TEST_DEPENDENCY_SOURCE_DIR}/exe_options.cmake)
    # Define our test options
    include(${EXE_OPTIONS_FILE})
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.probin_file=${TEST_NAME}.probin")
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
endfunction(add_test_v1)

# Verification test with resolutions of 8, 16, 32 (each test runs on maximum number of processes on node)
function(add_test_v2 TEST_NAME TEST_DEPENDENCY)
    foreach(RESOLUTION IN ITEMS 8 16 32)
      #message(STATUS "RESOLUTION=${RESOLUTION}")
      # Set variables for respective binary and source directories for the test
      set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
      set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME}/${RESOLUTION})
      set(TEST_DEPENDENCY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_DEPENDENCY})
      set(TEST_DEPENDENCY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_DEPENDENCY})
      # Make working directory for test
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
      # Gather all files in source directory for test
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
      # Copy files to test working directory
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
      # Get test options
      set(EXE_OPTIONS_FILE ${TEST_DEPENDENCY_SOURCE_DIR}/exe_options.cmake)
      # Define our test options
      include(${EXE_OPTIONS_FILE})
      # Set number of cells at runtime according to dimension
      if(${PELEC_DIM} EQUAL 3)
        set(NCELLS "${RESOLUTION} ${RESOLUTION} ${RESOLUTION}")
      elseif(${PELEC_DIM} EQUAL 2)
        set(NCELLS "${RESOLUTION} ${RESOLUTION}")
      elseif(${PELEC_DIM} EQUAL 1)
        set(NCELLS "${RESOLUTION}")
      endif()
      # Set some default runtime options for all tests in this category
      set(RUNTIME_OPTIONS "amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.probin_file=${TEST_NAME}.probin amr.n_cell=${NCELLS}")
      # Add test and actual test commands to CTest database
      add_test(${TEST_NAME}_${RESOLUTION} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCESSES} ${MPIEXEC_PREFLAGS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS}")
      # Set properties for test
      set_tests_properties(${TEST_NAME}_${RESOLUTION} PROPERTIES TIMEOUT 1500 PROCESSORS ${PROCESSES} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
    endforeach()
endfunction(add_test_v2)

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
add_test_r(mms-1d-1 4)
add_test_r(mms-2d-1 4)
add_test_r(mms-2d-2 4)
add_test_r(mms-3d-1 4)
add_test_r(mms-3d-2 4)
add_test_r(mms-3d-3 4)
add_test_r(mms-3d-4 1)
add_test_r(sod-3d-1 4)
add_test_r(tg-2d-1 4)
add_test_r(tg-3d-1 4)
add_test_r(tg-3d-2 4)

#=============================================================================
# Verification tests
#=============================================================================
if(ENABLE_VERIFICATION)
  add_test_v1(symmetry_3d mms-3d-1 4)
  add_test_v2(cns_no_amr_1d mms-1d-1)
  add_test_v2(cns_no_amr_2d mms-2d-1)
  add_test_v2(cns_no_amr_3d mms-3d-1)
  add_test_v2(cns_no_amr_mol_2d mms-2d-2)
  add_test_v2(cns_no_amr_mol_3d mms-3d-3)
  #add_test_v2(cns_les_noamr_3d mms-3d-3) # What MMS does this depend on?
  #add_test_v2(cns_amr_3d mms-3d-1) # This one takes a while
endif()

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit-tests-3d 1)

#=============================================================================
# Performance tests
#=============================================================================
