# target_link_libraries_system
#
# This function is similar to target_link_libraries but allows the includes
# determined from the library to be added as system includes to suppress
# warnings generated from those header files
#
# https://stackoverflow.com/questions/52135983/cmake-target-link-libraries-include-as-system-to-suppress-compiler-warnings/52136398#52136398
#
function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

macro(init_clang_tidy)
  if(PELEC_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
    if(CLANG_TIDY_EXE)
      message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    else()
      message(WARNING "clang-tidy not found.")
    endif()
  endif()
endmacro(init_clang_tidy)
