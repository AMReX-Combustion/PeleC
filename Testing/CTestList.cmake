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
    # Define our main run command
    set(RUN_COMMAND "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.probin_file=${TEST_NAME}.probin")
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${RUN_COMMAND} ${RUNTIME_OPTIONS} && nosetests ${TEST_NAME}.py")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
endfunction(add_test_v1)

# Verification test with multiple resolutions (each test runs on maximum number of processes on node)
function(add_test_v2 TEST_NAME TEST_DEPENDENCY)
    # Make sure run command is cleared before we construct it
    unset(MASTER_RUN_COMMAND)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(TEST_DEPENDENCY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_DEPENDENCY})
    set(TEST_DEPENDENCY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_DEPENDENCY})
    # Copy python file to test directory for running nosetests
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/test_second_order.py DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Get test dependency options (mainly just need the dimension value)
    set(EXE_OPTIONS_FILE ${TEST_DEPENDENCY_SOURCE_DIR}/exe_options.cmake)
    # Define our test options
    include(${EXE_OPTIONS_FILE})
    # Create list of resolutions we want to test with
    set(RESOLUTION_LIST 8 12 16 20)
    if(${PELEC_DIM} EQUAL 1)
      set(RESOLUTION_LIST 8 16 32 64)
    endif()
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET RESOLUTION_LIST -1 LAST_RESOLUTION_IN_LIST)
    # Create the commands to run for each resolution
    foreach(RESOLUTION IN LISTS RESOLUTION_LIST)
      # Make working directory for test
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${RESOLUTION})
      # Gather all files in source directory for test
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
      # Copy files to test working directory
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${RESOLUTION}/")
      # Set number of cells at runtime according to dimension
      if(${PELEC_DIM} EQUAL 3)
        set(NCELLS "${RESOLUTION} ${RESOLUTION} ${RESOLUTION}")
      elseif(${PELEC_DIM} EQUAL 2)
        set(NCELLS "${RESOLUTION} ${RESOLUTION}")
      elseif(${PELEC_DIM} EQUAL 1)
        set(NCELLS "${RESOLUTION}")
      endif()
      # Set the run command for this resolution
      set(RUN_COMMAND_${RESOLUTION} "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCESSES} ${MPIEXEC_PREFLAGS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${RESOLUTION}/${TEST_NAME}.i")
      # Set some runtime options for each resolution
      set(RUNTIME_OPTIONS_${RESOLUTION} "amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.probin_file=${TEST_NAME}.probin amr.n_cell=${NCELLS}")
      # Construct our large run command with everything &&'d together
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${RESOLUTION}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${RESOLUTION}} ${RUNTIME_OPTIONS_${RESOLUTION}}")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${RESOLUTION} EQUAL ${LAST_RESOLUTION_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    # Add test and actual test commands to CTest database (need to convert this to arrays for resolutions)
    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && nosetests test_second_order.py")
    # Set properties for test and make sure test dependencies have run before this test will run
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 7200 PROCESSORS ${PROCESSES} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
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
  #add_test_v3(cns_amr_3d mms-3d-1) # This one takes a while with AMR
endif()

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit-tests-3d 1)

#=============================================================================
# Performance tests
#=============================================================================
