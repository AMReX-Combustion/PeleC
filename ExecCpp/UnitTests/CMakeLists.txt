set(pelec_exe_name pelec_unit_tests)
set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL GammaLaw)
set(PELEC_CHEMISTRY_MODEL Null)
set(PELEC_TRANSPORT_MODEL Constant)

include(${CMAKE_SOURCE_DIR}/CMake/BuildPeleCExe.cmake)
build_pelec_exe(${pelec_exe_name})
target_sources(${pelec_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  )

target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)
target_link_libraries(${pelec_exe_name} PRIVATE gtest)
