if(PELEC_ENABLE_ALL_WARNINGS)
  list(APPEND PELEC_CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
      list(APPEND PELEC_CXX_FLAGS "-faligned-new"
                                  "-Wunreachable-code"
                                  "-Wnull-dereference"
                                  "-Wfloat-conversion"
                                  "-Wshadow"
                                  "-Woverloaded-virtual")
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    list(APPEND PELEC_CXX_FLAGS "-diag-disable:11074,11076,10397,15335")
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 9.0.0.9000037)
      list(APPEND PELEC_CXX_FLAGS "-Wno-missing-braces")
    endif()
  endif()
endif()

if(PELEC_ENABLE_WERROR)
  list(APPEND PELEC_CXX_FLAGS "-Werror")
endif()

separate_arguments(PELEC_CXX_FLAGS)
target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${PELEC_CXX_FLAGS}>)
