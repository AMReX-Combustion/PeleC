set(PELE_PHYSICS_EOS_MODEL GammaLaw)
set(PELE_PHYSICS_CHEMISTRY_MODEL Null)
set(PELE_PHYSICS_TRANSPORT_MODEL Constant)
set(PELE_PHYSICS_ENABLE_SOOT OFF)
set(PELE_PHYSICS_ENABLE_SPRAY OFF)
set(PELE_PHYSICS_SPRAY_FUEL_NUM 0)
include(BuildExeAndLib)

target_sources(${pelec_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  )

if(PELE_ENABLE_CUDA)
  set_source_files_properties(unit-tests-main.cpp test-config.cpp PROPERTIES LANGUAGE CUDA)
endif()

target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)
target_link_libraries(${pelec_exe_name} PRIVATE gtest)
