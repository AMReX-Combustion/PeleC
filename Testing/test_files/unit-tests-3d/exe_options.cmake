#User-specific source files
#list(APPEND PELEC_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Testing/test_files/unit-tests-3d/probdata.f90)
#list(APPEND PELEC_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Testing/test_files/unit-tests-3d/pmf_generic.f90)
#list(APPEND PELEC_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Testing/test_files/unit-tests-3d/Prob_nd.F90)
#list(APPEND PELEC_EXTRA_SOURCES ${CMAKE_SOURCE_DIR}/Testing/test_files/unit-tests-3d/bc_fill_nd.F90)

#Compile-time options for executable
set(PELEC_DIM 3)
set(PELEC_ENABLE_EB ON)
set(PELEC_ENABLE_MASA OFF)
set(PELEC_ENABLE_REACTIONS OFF)
set(PELEC_ENABLE_MOL OFF)
set(PELEC_ENABLE_PARTICLES OFF)
set(PELEC_EOS_MODEL Fuego)
set(PELEC_REACTIONS_MODEL Fuego)
set(PELEC_CHEMISTRY_MODEL LiDryer)
set(PELEC_TRANSPORT_MODEL EGLib)
