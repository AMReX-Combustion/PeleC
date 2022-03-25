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

macro(init_code_checks)
  if(PELEC_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
    if(CLANG_TIDY_EXE)
      message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    else()
      message(WARNING "clang-tidy not found.")
    endif()
  endif()

  if(PELEC_ENABLE_CPPCHECK)
    find_program(CPPCHECK_EXE NAMES "cppcheck")
    if(CPPCHECK_EXE)
      message(STATUS "cppcheck found: ${CPPCHECK_EXE}")
      include(ProcessorCount)
      ProcessorCount(NP)
      if(NP EQUAL 0)
        set(NP 1)
      endif()
      add_custom_target(cppcheck ALL
          COMMAND ${CMAKE_COMMAND} -E echo "Running cppcheck on project using ${NP} cores..."
          COMMAND ${CMAKE_COMMAND} -E make_directory cppcheck
          # cppcheck ignores -isystem directories, so we change them to regular -I include directories (with no spaces either)
          COMMAND sed "s/isystem /I/g" compile_commands.json > cppcheck/cppcheck_compile_commands.json
          COMMAND ${CPPCHECK_EXE} --version
          COMMAND ${CPPCHECK_EXE} --template=gcc --inline-suppr --suppress=unusedFunction --std=c++14 --language=c++ --enable=all --project=cppcheck/cppcheck_compile_commands.json -i ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src -i ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest --output-file=cppcheck/cppcheck-full-report.txt -j ${NP}
          # Currently we filter analysis from submodules after cppcheck has run
          #COMMAND awk -v nlines=2 "/Submodules\/AMReX/ || /Submodules\/GoogleTest/ {for (i=0; i<nlines; i++) {getline}; next} 1" < cppcheck/cppcheck-full-report.txt > cppcheck/cppcheck-short-report.txt
          COMMAND awk -v nlines=2 "/Submodules/ {for (i=0; i<nlines; i++) {getline}; next} 1" < cppcheck/cppcheck-full-report.txt > cppcheck/cppcheck-short-report.txt
          COMMAND cat cppcheck/cppcheck-short-report.txt | egrep "information:|error:|performance:|portability:|style:|warning:" | sort | uniq > cppcheck-warnings.txt
          COMMAND printf "Warnings: " >> cppcheck-warnings.txt
          COMMAND cat cppcheck-warnings.txt | awk "END{print NR-1}" >> cppcheck-warnings.txt
          COMMENT "Run cppcheck on project compile_commands.json"
          BYPRODUCTS cppcheck-warnings.txt
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          VERBATIM USES_TERMINAL
      )
    else()
      message(WARNING "cppcheck not found.")
    endif()
  endif()
endmacro(init_code_checks)
