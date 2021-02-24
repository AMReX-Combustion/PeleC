set(PELEC_ENABLE_EB ON)
set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL GammaLaw)
set(PELEC_CHEMISTRY_MODEL air)
set(PELEC_TRANSPORT_MODEL Constant)

#Generic stuff for every executable
get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(pelec_exe_name PeleC-${DIR_NAME})
include(BuildPeleCExe)
build_pelec_exe(${pelec_exe_name})
