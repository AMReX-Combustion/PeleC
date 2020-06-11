#=============================================================================
# PeleC Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

# Set location of gold files according to system/compiler/compiler_version
set(FCOMPARE_GOLD_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/PeleCGoldFiles/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})

if(TEST_WITH_FCOMPARE)
  message(STATUS "Test golds directory for fcompare: ${FCOMPARE_GOLD_FILES_DIRECTORY}")
endif()

# Have CMake discover the number of cores on the node
#include(ProcessorCount)
#ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/ExecCpp/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/ExecCpp/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/ExecCpp/RegTests/${TEST_EXE_DIR}/pelec_${TEST_EXE_DIR})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${FCOMPARE_GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(TEST_WITH_FCOMPARE)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/${AMREX_SUBMOD_LOCATION}/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    # Copy files to test working directory
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests in this category
    set(RUNTIME_OPTIONS "max_step=10 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.signal_handling=0")
    # Use fcompare to test diffs in plots against gold files
    if(TEST_WITH_FCOMPARE)
      set(FCOMPARE_COMMAND "&& ${FCOMPARE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
    if(PELEC_ENABLE_MPI)
      set(NP 4)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Place the exe in the correct working directory
    #set_target_properties(${pelec_exe_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/")
    # Add test and actual test commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i ${RUNTIME_OPTIONS} ${FEXTREMA_COMMAND} ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 1500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression")
endfunction(add_test_r)

# Verification test with 1 resolution
#function(add_test_v1 TEST_NAME TEST_DEPENDENCY NP)
#    # Set variables for respective binary and source directories for the test
#    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
#    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
#    set(TEST_DEPENDENCY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_DEPENDENCY})
#    set(TEST_DEPENDENCY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_DEPENDENCY})
#    # Make working directory for test
#    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
#    # Gather all files in source directory for test
#    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
#    # Copy files to test working directory
#    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
#    # Get test options
#    set(EXE_OPTIONS_FILE ${TEST_DEPENDENCY_SOURCE_DIR}/exe_options.cmake)
#    # Define our test options
#    include(${EXE_OPTIONS_FILE})
#    if(PELEC_ENABLE_MPI)
#      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
#    else()
#      unset(MPI_COMMANDS)
#    endif()
#    # Define our main run command
#    set(RUN_COMMAND "rm mmslog datlog || true && ${MPI_COMMANDS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.i")
#    # Set some default runtime options for all tests in this category
#    set(RUNTIME_OPTIONS "amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amrex.signal_handling=0")
#    # Add test and actual test commands to CTest database
#    add_test(${TEST_NAME} sh -c "${RUN_COMMAND} ${RUNTIME_OPTIONS} && nosetests ${TEST_NAME}.py")
#    # Set properties for test
#    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 1500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
#endfunction(add_test_v1)

# Verification test with multiple resolutions (each test runs on maximum number of processes on node)
#function(add_test_v2 TEST_NAME TEST_DEPENDENCY LIST_OF_GRID_SIZES)
#    # Make sure run command is cleared before we construct it
#    unset(MASTER_RUN_COMMAND)
#    # Set variables for respective binary and source directories for the test
#    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
#    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
#    set(TEST_DEPENDENCY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_DEPENDENCY})
#    set(TEST_DEPENDENCY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_DEPENDENCY})
#    # Copy python file to test directory for running nosetests
#    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/test_second_order.py DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
#    # Get test dependency options (mainly just need the dimension value)
#    set(EXE_OPTIONS_FILE ${TEST_DEPENDENCY_SOURCE_DIR}/exe_options.cmake)
#    # Define our test options
#    include(${EXE_OPTIONS_FILE})
#    # Get last item in resolution list so we can find out when we are on the last item in our loop
#    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
#    # Create the commands to run for each resolution
#    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
#      # Make working directory for test
#      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
#      # Gather all files in source directory for test
#      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
#      # Copy files to test working directory
#      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
#      # Set number of cells at runtime according to dimension
#      if(${PELEC_DIM} EQUAL 3)
#        set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
#      elseif(${PELEC_DIM} EQUAL 2)
#        set(NCELLS "${GRID_SIZE} ${GRID_SIZE}")
#      elseif(${PELEC_DIM} EQUAL 1)
#        set(NCELLS "${GRID_SIZE}")
#      endif()
#      # Set the command to delete files from previous test runs in each resolution
#      set(DELETE_PREVIOUS_FILES_COMMAND "rm mmslog datlog || true")
#      if(PELEC_ENABLE_MPI)
#        set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCESSES} ${MPIEXEC_PREFLAGS}")
#      else()
#        unset(MPI_COMMANDS)
#      endif()
#      # Set the run command for this resolution
#      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${TEST_DEPENDENCY_BINARY_DIR}/PeleC-${TEST_DEPENDENCY} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.i")
#      # Set some runtime options for each resolution
#      set(RUNTIME_OPTIONS_${GRID_SIZE} "amrex.signal_handling=0 amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1 amr.n_cell=${NCELLS}")
#      # Construct our large run command with everything &&'d together
#      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
#      string(APPEND MASTER_RUN_COMMAND " && ")
#      string(APPEND MASTER_RUN_COMMAND "${DELETE_PREVIOUS_FILES_COMMAND}")
#      string(APPEND MASTER_RUN_COMMAND " && ")
#      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}}")
#      # Add another " && " unless we are on the last resolution in the list
#      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
#        string(APPEND MASTER_RUN_COMMAND " && ")
#      endif()
#    endforeach()
#    # Set list of images to be uploaded for verification
#    set(IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/p_error.png ${CURRENT_TEST_BINARY_DIR}/rho_error.png ${CURRENT_TEST_BINARY_DIR}/u_error.png)
#    if(${PELEC_DIM} EQUAL 3)
#      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png ${CURRENT_TEST_BINARY_DIR}/w_error.png)
#    elseif(${PELEC_DIM} EQUAL 2)
#      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png)
#    endif()
#    # Add test and actual test commands to CTest database (need to convert this to arrays for resolutions)
#    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && nosetests test_second_order.py")
#    # Set properties for test and make sure test dependencies have run before this test will run
#    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 7200 PROCESSORS ${PROCESSES} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification" ATTACHED_FILES "${IMAGES_TO_UPLOAD}" FIXTURES_REQUIRED ${TEST_DEPENDENCY})
#endfunction(add_test_v2)

# Standard unit test
function(add_test_u TEST_NAME)
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/ExecCpp/UnitTests/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/ExecCpp/UnitTests/${TEST_NAME})
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    set(NP 1)
    if(PELEC_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    # Add test and commands to CTest database
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 500 PROCESSORS ${NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(pmf-1 PMF)
add_test_r(tg-1 TG)
add_test_r(tg-2 TG)
add_test_r(tg-3 TG)
add_test_r(tg-4 TG)
add_test_r(hit-1 HIT)
add_test_r(hit-2 HIT)
add_test_r(hit-3 HIT)
add_test_r(sod-1 Sod)

if(PELEC_ENABLE_MASA)
  add_test_r(mms-1 MMS)
  add_test_r(mms-2 MMS)
  #add_test_r(mms-3 MMS)
  add_test_r(mms-4 MMS)
  add_test_r(mms-5 MMS)
  #add_test_r(ebmms-1 EB_MMS)
endif()

#=============================================================================
# Verification tests
#=============================================================================
#if(PELEC_ENABLE_MASA)
#  add_test_v1(symmetry mms-1 4)
#  add_test_v1(eb_symmetry ebmms-1 4)
#
#  # Create list of resolutions we want to test with; make sure to pass it as a string in quotes
#  set(LIST_OF_GRID_SIZES 8 12 16 20)
#  add_test_v2(cns_no_amr mms-1 "${LIST_OF_GRID_SIZES}")
#  add_test_v2(cns_no_amr_mol mms-3 "${LIST_OF_GRID_SIZES}")
#  #add_test_v3(cns_amr mms-1 "${LIST_OF_GRID_SIZES}") # This one takes a while with AMR
#
#  set(LIST_OF_GRID_SIZES 8 12 16 24)
#  add_test_v2(cns_les_no_amr mms-5 "${LIST_OF_GRID_SIZES}")
#endif()

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(pelec_unit_tests)

#=============================================================================
# Performance tests
#=============================================================================
