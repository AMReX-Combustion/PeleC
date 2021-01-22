set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL GammaLaw)
set(PELEC_CHEMISTRY_MODEL Null)
set(PELEC_TRANSPORT_MODEL Constant)

get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(pelec_exe_name PeleC-${DIR_NAME})
include(BuildPeleCExe)
build_pelec_exe(${pelec_exe_name})
target_sources(${pelec_exe_name}
  PUBLIC
  unit-tests-main.cpp
  test-config.cpp
  )

if(PELEC_ENABLE_CUDA)
  set_source_files_properties(unit-tests-main.cpp test-config.cpp PROPERTIES LANGUAGE CUDA)
endif()

target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/GoogleTest/googletest/include)
target_link_libraries(${pelec_exe_name} PRIVATE gtest)
