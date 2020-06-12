set(pelec_exe_name pelec_unit_tests)

#Compile-time options for executable
set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL GammaLaw)
set(PELEC_REACTIONS_MODEL Null)
set(PELEC_CHEMISTRY_MODEL Null)
set(PELEC_TRANSPORT_MODEL Constant)

add_executable(${pelec_exe_name} "")
target_sources(${pelec_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  prob.cpp
  prob.H
  prob_parm.H
  )

target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildPeleCExe.cmake)
build_pelec_exe(${pelec_exe_name})

target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)
target_link_libraries(${pelec_exe_name} PRIVATE gtest)
