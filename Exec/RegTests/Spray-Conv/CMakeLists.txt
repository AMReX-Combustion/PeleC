set(PELEC_ENABLE_EB OFF)
set(PELEC_ENABLE_PARTICLES ON)
set(PELEC_EOS_MODEL Fuego)
set(PELEC_CHEMISTRY_MODEL decane_3sp)
set(PELEC_TRANSPORT_MODEL Simple)
set(PELEMP_SPRAY_FUEL_NUM 1)

#Generic stuff for every executable
get_filename_component(DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(pelec_exe_name PeleC-${DIR_NAME})
include(BuildPeleCExe)
build_pelec_exe(${pelec_exe_name})
