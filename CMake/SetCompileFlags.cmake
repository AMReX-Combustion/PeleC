# Logic for handling warnings
list(APPEND PELEC_CXX_FLAGS "-Wno-pass-failed") # Ignore loop not vectorized warnings
if(PELEC_ENABLE_ALL_WARNINGS)
  # GCC, Clang, and Intel seem to accept these
  list(APPEND PELEC_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # Intel always reports some diagnostics we don't necessarily care about
    list(APPEND PELEC_CXX_FLAGS "-diag-disable:11074,11076")
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    # Avoid notes about -faligned-new with GCC > 7
    list(APPEND PELEC_CXX_FLAGS "-faligned-new")
  endif()
endif()

# Add our extra flags according to language
separate_arguments(PELEC_CXX_FLAGS)
target_compile_options(${pelec_exe_name} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PELEC_CXX_FLAGS}>)
