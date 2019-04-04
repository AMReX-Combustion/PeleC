function(set_compile_flags)
  # Clear our local compiler flag variables
  set(CXX_FLAGS "")
  set(C_FLAGS "")
  set(Fortran_FLAGS "")

  # Add any extra flags based on compiler and/or OS
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    list(APPEND CXX_FLAGS "")
    list(APPEND C_FLAGS "")
    list(APPEND Fortran_FLAGS "-ffree-line-length-none"
                              "-ffixed-line-length-none"
                              "-fno-range-check"
                              "-fno-second-underscore")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    list(APPEND CXX_FLAGS "")
    list(APPEND C_FLAGS "")
    list(APPEND Fortran_FLAGS "-ffree-line-length-none"
                              "-ffixed-line-length-none"
                              "-fno-range-check"
                              "-fno-second-underscore")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    list(APPEND CXX_FLAGS "")
    list(APPEND C_FLAGS "")
    list(APPEND Fortran_FLAGS "-ffree-line-length-none"
                              "-ffixed-line-length-none"
                              "-fno-range-check"
                              "-fno-second-underscore")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    list(APPEND CXX_FLAGS "-restrict")
    list(APPEND C_FLAGS "-restrict")
    list(APPEND Fortran_FLAGS "")
  endif()
  
  # Logic for managing warnings
  if(ENABLE_ALL_WARNINGS)
    # GCC, Clang, and Intel seem to accept these
    list(APPEND CXX_FLAGS "-Wall" "-Wextra" "-pedantic")
    list(APPEND C_FLAGS "-Wall" "-Wextra" "-pedantic")
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      # ifort doesn't like -Wall
      list(APPEND Fortran_FLAGS "-Wall")
    else()
      # Intel reports some diagnostics we don't necessarily care about
      list(APPEND CXX_FLAGS "-diag-disable:11074,11076")
    endif()
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
      # Avoid notes about -faligned-new with GCC > 7
      list(APPEND CXX_FLAGS "-faligned-new")
    endif()
  else()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
      # Avoid warning about explicit vectorization not happening in AMReX
      list(APPEND CXX_FLAGS "-Wno-pass-failed")
      list(APPEND C_FLAGS "-Wno-pass-failed")
    endif()
  endif()

  # Since we created the flags as lists, we replace the semicolons when it gets converted to a string 
  string(REPLACE ";" " " CXX_FLAGS "${CXX_FLAGS}")
  string(REPLACE ";" " " C_FLAGS "${C_FLAGS}")
  string(REPLACE ";" " " Fortran_FLAGS "${Fortran_FLAGS}")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS}" PARENT_SCOPE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_FLAGS}" PARENT_SCOPE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${Fortran_FLAGS}" PARENT_SCOPE)

  message("-- CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("-- CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
  message("-- CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")

endfunction(set_compile_flags)
