if(NOT "${TESTING_ROOT_DIR}" STREQUAL "")
  message("Testing root directory is ${TESTING_ROOT_DIR}")
else()
  message(FATAL_ERROR "You need to set the TESTING_ROOT_DIR variable. CMake will exit." )
endif()

if(NOT "${HOST_NAME}" STREQUAL "")
  message("Hostname is ${HOST_NAME}")
else()
  message(FATAL_ERROR "You need to set the HOST_NAME variable. CMake will exit." )
endif()

if(NOT "${PELEC_DIR}" STREQUAL "")
  message("PELEC_DIR is ${PELEC_DIR}")
else()
  message(FATAL_ERROR "You need to set the PELEC_DIR variable. CMake will exit." )
endif()

# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

# Set important configuration variables
set(CTEST_SITE "${HOST_NAME}")
set(CTEST_BUILD_NAME "${CMAKE_SYSTEM_NAME}${EXTRA_BUILD_NAME}")
set(CTEST_SOURCE_DIRECTORY "${PELEC_DIR}")
set(CTEST_BINARY_DIRECTORY "${PELEC_DIR}/build")
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
find_program(CTEST_GIT_COMMAND NAMES git)
find_program(MAKE NAMES make)

# Add parallelism capability to testing
include(ProcessorCount)
ProcessorCount(NP)
message(STATUS "\nNumber of processors detected: ${NP}")
set(CTEST_BUILD_FLAGS "-j${NP}")
set(CTEST_PARALLEL_LEVEL ${NP})

# Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_INIT_SUBMODULES TRUE)

# Configure Command
set(CTEST_CONFIGURE_COMMAND "cmake ${CMAKE_CONFIGURE_ARGS} -DENABLE_TESTS:BOOL=ON ${CTEST_SOURCE_DIRECTORY}")

# Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${CTEST_BUILD_FLAGS}")

# -----------------------------------------------------------
# -- Run CTest
# -----------------------------------------------------------

#ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
 
message("\n -- Start dashboard - ${CTEST_BUILD_NAME} --")
ctest_start("Nightly" TRACK "Nightly")

message("\n -- Update - ${CTEST_BUILD_NAME} --")
ctest_update(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE result)
message(" -- Update exit code = ${result} --")
if(result GREATER -1)
  message("\n -- Configure - ${CTEST_BUILD_NAME} --")
  ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE result)
  message(" -- Configure exit code = ${result} --")
  if(result EQUAL 0)
    message("\n -- Build - ${CTEST_BUILD_NAME} --")
    ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
    ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE result)
    message(" -- Build exit code = ${result} --")
    if(result EQUAL 0)
      # Need to have TMPDIR set to disk on certain NREL machines for building so builds
      # do not run out of space but unset when running to stop OpenMPI from complaining
      if(UNSET_TMPDIR_VAR)
        message("Clearing TMPDIR variable...")
        unset(ENV{TMPDIR})
      endif()
      message("\n -- Test - ${CTEST_BUILD_NAME} --")
      ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}"
                 PARALLEL_LEVEL ${CTEST_PARALLEL_LEVEL}
                 RETURN_VALUE result)
      message(" -- Test exit code = ${result} --")
    endif()
  endif()
endif()

message("\n -- Submit - ${CTEST_BUILD_NAME} --")
set(CTEST_NOTES_FILES "${TEST_LOG}")
if(HAVE_STATIC_ANALYSIS_OUTPUT)
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} "${STATIC_ANALYSIS_LOG}")
endif()
ctest_submit(RETRY_COUNT 20
             RETRY_DELAY 20
             RETURN_VALUE result)
message(" -- Submit exit code = ${result} --")

message("\n -- Finished - ${CTEST_BUILD_NAME} --")
