if(PELEC_ENABLE_ALL_WARNINGS)
  list(APPEND PELEC_CXX_FLAGS "-Wall" "-Wextra" "-pedantic" "-Wno-unused-function")
  if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
      list(APPEND PELEC_CXX_FLAGS "-faligned-new"
                                  "-Wunreachable-code"
                                  "-Wnull-dereference"
                                  "-Wfloat-conversion"
                                  "-Wshadow"
                                  "-Woverloaded-virtual")
      if(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")
        list(APPEND PELEC_CXX_FLAGS "-Wno-pass-failed" "-Wno-c++17-extensions")
      endif()
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    list(APPEND PELEC_CXX_FLAGS "-diag-disable:11074,11076,10397,15335")
  endif()
endif()

if(PELEC_ENABLE_WERROR)
  list(APPEND PELEC_CXX_FLAGS "-Werror")
endif()

separate_arguments(PELEC_CXX_FLAGS)
target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${PELEC_CXX_FLAGS}>)
