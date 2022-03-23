#=============================================================================
# PeleC Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

if(PELEC_ENABLE_FCOMPARE_FOR_TESTS)
  if("${PELEC_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, PELEC_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELEC_REFERENCE_GOLDS_DIRECTORY})
      set(GOLD_FILES_DIRECTORY ${PELEC_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for reference gold files does not exist: ${PELEC_REFERENCE_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

if(PELEC_SAVE_GOLDS)
  if("${PELEC_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, PELEC_SAVED_GOLDS_DIRECTORY must be set and exist")
  else()
    if(EXISTS ${PELEC_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${PELEC_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${PELEC_SAVED_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

macro(setup_test)
    # Set variables for respective binary and source directories for the test
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/tests/${TEST_NAME})
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/RegTests/${TEST_EXE_DIR}/PeleC-${TEST_EXE_DIR})
    # Gold files should be submodule organized by machine and compiler (these are output during configure)
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_EXE_DIR}/tests/${TEST_NAME}/plt00010)
    # Test plot is currently expected to be after 10 steps
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    # Find fcompare
    if(PELEC_ENABLE_FCOMPARE_FOR_TESTS)
      set(FCOMPARE ${CMAKE_BINARY_DIR}/Submodules/AMReX/Tools/Plotfile/fcompare)
    endif()
    # Make working directory for test
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    # Gather all files in source directory for test
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
    # Copy files to test working directory
    file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    # Set some default runtime options for all tests
    set(RUNTIME_OPTIONS "amr.plot_file=plt amr.checkpoint_files_output=0 amr.plot_files_output=1")
    if(PELEC_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(PELEC_ENABLE_MPI)
      if(PELEC_ENABLE_CUDA)
        set(PELEC_NP 2) # 1 rank per GPU on Eagle
      else()
        set(PELEC_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELEC_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(PELEC_NP 1)
      unset(MPI_COMMANDS)
    endif()
    # Use fcompare to test diffs in plots against gold files
    if(PELEC_ENABLE_FCOMPARE_FOR_TESTS)
      if(PELEC_ENABLE_CUDA)
        set(FCOMPARE_TOLERANCE "-r 1e-12 --abs_tol 1.0e-12")
      endif()
      set(FCOMPARE_COMMAND "&& ${MPI_COMMANDS} ${FCOMPARE} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
    if(PELEC_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_EXE_DIR}/${TEST_NAME}/")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUNTIME_OPTIONS "max_step=10 ${RUNTIME_OPTIONS}")
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} > ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELEC_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "regression" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no-ci")
endfunction(add_test_re)

# Regression tests expected to fail
function(add_test_rf TEST_NAME TEST_EXE_DIR)
    add_test_r(${TEST_NAME} ${TEST_EXE_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES WILL_FAIL TRUE)
endfunction(add_test_rf)

# Verification test with 1 resolution
function(add_test_v1 TEST_NAME TEST_EXE_DIR)
    setup_test()
    set(RUN_COMMAND "rm -f mmslog datlog && ${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp")
    add_test(${TEST_NAME} sh -c "${RUN_COMMAND} ${RUNTIME_OPTIONS} > ${TEST_NAME}.log && nosetests ${TEST_NAME}.py")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELEC_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "verification" ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log")
endfunction(add_test_v1)

# Verification test with multiple resolutions (each test runs on maximum number of processes on node)
function(add_test_v2 TEST_NAME TEST_EXE_DIR LIST_OF_GRID_SIZES)
    setup_test()
    unset(MASTER_RUN_COMMAND)
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
    # Create the commands to run for each resolution
    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*.dat" "${CURRENT_TEST_SOURCE_DIR}/*.py")
      file(COPY ${CURRENT_TEST_SOURCE_DIR}/${TEST_NAME}.inp DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      if(${PELEC_DIM} EQUAL 3)
        set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
      elseif(${PELEC_DIM} EQUAL 2)
        set(NCELLS "${GRID_SIZE} ${GRID_SIZE}")
      elseif(${PELEC_DIM} EQUAL 1)
        set(NCELLS "${GRID_SIZE}")
      endif()
      set(DELETE_PREVIOUS_FILES_COMMAND "rm -f mmslog datlog")
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.inp")
      set(RUNTIME_OPTIONS_${GRID_SIZE} "${RUNTIME_OPTIONS} amr.n_cell=${NCELLS}")
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${DELETE_PREVIOUS_FILES_COMMAND}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}} > ${TEST_NAME}-${GRID_SIZE}.log")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    set(IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/p_error.png ${CURRENT_TEST_BINARY_DIR}/rho_error.png ${CURRENT_TEST_BINARY_DIR}/u_error.png)
    if(${PELEC_DIM} EQUAL 3)
      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png ${CURRENT_TEST_BINARY_DIR}/w_error.png)
    elseif(${PELEC_DIM} EQUAL 2)
      list(APPEND IMAGES_TO_UPLOAD ${CURRENT_TEST_BINARY_DIR}/v_error.png)
    endif()
    add_test(${TEST_NAME} sh -c "${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && nosetests ${CMAKE_CURRENT_SOURCE_DIR}/test_second_order.py")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 18000 PROCESSORS ${PELEC_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification;no-ci" ATTACHED_FILES "${IMAGES_TO_UPLOAD}")
endfunction(add_test_v2)

# Standard unit test
function(add_test_u TEST_NAME)
    set(CURRENT_TEST_EXE ${CMAKE_BINARY_DIR}/Exec/UnitTests/PeleC-UnitTests)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/Exec/UnitTests/${TEST_EXE_DIR})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_BINARY_DIR}/Exec/UnitTests/${TEST_EXE_DIR})
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    set(PELEC_NP 1)
    if(PELEC_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PELEC_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} sh -c "${MPI_COMMANDS} ${CURRENT_TEST_EXE} ${MPIEXEC_POSTFLAGS}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 1800 PROCESSORS ${PELEC_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/" LABELS "unit")
endfunction(add_test_u)

#=============================================================================
# Regression tests
#=============================================================================

# Run in CI
add_test_r(multispecsod-1 MultiSpecSod)
add_test_r(pmf-lidryer-arkode PMF)
add_test_r(pmf-srk-1 PMF-SRK)
add_test_r(tg-1 TG)
add_test_r(tg-2 TG)
add_test_r(hit-1 HIT)
add_test_r(hit-2 HIT)
add_test_r(hit-3 HIT)
add_test_r(sod-1 Sod)
add_test_r(sod-2 Sod)
add_test_r(eb-c4 EB-C4-5)
add_test_r(eb-c5 EB-C4-5)
add_test_r(eb-c9 EB-C9)
add_test_r(eb-c10 EB-C10)
add_test_r(eb-c11 EB-C11)
add_test_r(eb-c12 EB-C12)
# add_test_r(eb-c14 EB-C14) # disable due to FPE in ghost cells
add_test_r(eb-converging-nozzle EB-ConvergingNozzle)
if(PELEC_ENABLE_MASA)
  add_test_r(mms-3 MMS)
  add_test_r(mms-1 MMS)
  add_test_r(mms-2 MMS)
  add_test_r(mms-4 MMS)
  add_test_r(ebmms-1 EB-C1)
  if(PELEC_DIM GREATER 2)
    add_test_r(mms-5 MMS)
  endif()
endif()

# Not run in CI
if(PELEC_DIM GREATER 2)
  add_test_re(pmf-lidryer-rk64 PMF)
  add_test_re(pmf-lidryer-cvode PMF)
  add_test_re(tg-3 TG)
  add_test_re(tg-4 TG)
  add_test_re(sedov-1 Sedov)
  add_test_re(shu-osher-1 Shu-Osher)
  add_test_re(zerod-1 zeroD)
  add_test_re(eb-c3 EB-C3)
  add_test_re(eb-c7 EB-C7)
  add_test_re(eb-c8 EB-C8)
  add_test_re(eb-bluffbody-1 EB-BluffBody)
  add_test_re(eb-plane-1 EB-Plane)
endif()

#=============================================================================
# Verification tests
#=============================================================================
if(PELEC_ENABLE_MASA AND (PELEC_DIM GREATER 2))
  add_test_v1(symmetry MMS)
  add_test_v1(eb-symmetry EB-C1)

  # Created an option to exclude these tests when using ASAN and UBSAN because they take too long
  if(NOT PELEC_ENABLE_SANITIZE_FOR_TESTS)
    # Create list of resolutions we want to test with; make sure to pass it as a string in quotes
    set(LIST_OF_GRID_SIZES 12 16 24)
    add_test_v2(cns-no-amr MMS "${LIST_OF_GRID_SIZES}")
    add_test_v2(cns-no-amr-mol MMS "${LIST_OF_GRID_SIZES}")
    add_test_v2(cns-les-no-amr MMS "${LIST_OF_GRID_SIZES}")
    #add_test_v3(cns-amr MMS "${LIST_OF_GRID_SIZES}") # This one takes a while with AMR
  endif()
endif()

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit-tests)

#=============================================================================
# Performance tests
#=============================================================================
