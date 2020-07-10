set(pelec_exe_name pelec_EB_MMS)

#Compile-time options for executable
set(PELEC_ENABLE_EB ON)
set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL GammaLaw)
set(PELEC_REACTIONS_MODEL Null)
set(PELEC_CHEMISTRY_MODEL Null)
set(PELEC_TRANSPORT_MODEL Constant)

add_executable(${pelec_exe_name} "")
target_sources(${pelec_exe_name}
   PRIVATE
     prob_parm.H
     prob.H
     prob.cpp
)

target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

include(${CMAKE_SOURCE_DIR}/CMake/BuildPeleCExe.cmake)
build_pelec_exe(${pelec_exe_name})
