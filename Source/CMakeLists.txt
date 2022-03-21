set(pelec_lib_name "PeleCLib")

add_library(${pelec_lib_name} OBJECT)

if(CLANG_TIDY_EXE)
  set_target_properties(${pelec_lib_name} PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
endif()

target_include_directories(${pelec_lib_name} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
set(PELE_PHYSICS_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/PelePhysics/${pelec_lib_name})

set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${pelec_lib_name})

include(SetPeleCCompileFlags)

add_subdirectory(${SRC_DIR}/Params ${BIN_DIR}/Params/${pelec_lib_name})

target_include_directories(${pelec_lib_name} PUBLIC "${PELE_PHYSICS_SRC_DIR}/Source")

target_sources(${pelec_lib_name}
  PUBLIC
  ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.cpp
  ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.H)
target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow)

target_sources(${pelec_lib_name}
  PUBLIC
  ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.cpp
  ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.H
  ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManagerBCFill.H)
  target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager)

target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Reactions)
target_link_libraries(${pelec_lib_name} PUBLIC sundials_arkode sundials_cvode)

target_sources(${pelec_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_Sundials.H
                                        ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_Sundials.cpp
                                        ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.cpp
                                        ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.H
                                        ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_SUNMemory.cpp
                                        ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_SUNMemory.H)
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_Sundials.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_Sundials.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_NVector_MultiFab.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_SUNMemory.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS/AMReX_SUNMemory.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")

target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS)

if(PELEC_ENABLE_CUDA)
  target_link_libraries(${pelec_lib_name} PUBLIC sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
elseif(PELEC_ENABLE_HIP)
  target_link_libraries(${pelec_lib_name} PUBLIC sundials_nvechip)
elseif(PELEC_ENABLE_DPCPP)
  target_link_libraries(${pelec_lib_name} PUBLIC sundials_nvecsycl)
endif()

#target_sources(${pelec_lib_name}
#   PUBLIC
#     ${SRC_DIR}/Advance.cpp
#     ${SRC_DIR}/BCfill.cpp
#     ${SRC_DIR}/Bld.cpp
#     ${SRC_DIR}/Constants.H
#     ${SRC_DIR}/Derive.H
#     ${SRC_DIR}/Derive.cpp
#     ${SRC_DIR}/Diffterm.H
#     ${SRC_DIR}/Diffterm.cpp
#     ${SRC_DIR}/Diffusion.H
#     ${SRC_DIR}/Diffusion.cpp
#     ${SRC_DIR}/EB.H
#     ${SRC_DIR}/EB.cpp
#     ${SRC_DIR}/EBStencilTypes.H
#     ${SRC_DIR}/External.cpp
#     ${SRC_DIR}/Filter.H
#     ${SRC_DIR}/Filter.cpp
#     ${SRC_DIR}/Forcing.cpp
#     ${SRC_DIR}/GradUtil.H
#     ${SRC_DIR}/GradUtil.cpp
#     ${SRC_DIR}/Hydro.H
#     ${SRC_DIR}/Hydro.cpp
#     ${SRC_DIR}/Geometry.H
#     ${SRC_DIR}/Geometry.cpp
#     ${SRC_DIR}/Godunov.H
#     ${SRC_DIR}/Godunov.cpp
#     ${SRC_DIR}/PLM.H
#     ${SRC_DIR}/PPM.H
#     ${SRC_DIR}/PPM.cpp
#     ${SRC_DIR}/InitEB.cpp
#     ${SRC_DIR}/IndexDefines.H
#     ${SRC_DIR}/IndexDefines.cpp
#     ${SRC_DIR}/IO.H
#     ${SRC_DIR}/IO.cpp
#     ${SRC_DIR}/LES.H
#     ${SRC_DIR}/LES.cpp
#     ${SRC_DIR}/MOL.H
#     ${SRC_DIR}/MOL.cpp
#     ${SRC_DIR}/Particle.cpp
#     ${SRC_DIR}/PeleC.H
#     ${SRC_DIR}/PeleC.cpp
#     ${SRC_DIR}/Problem.H
#     ${SRC_DIR}/ProblemDerive.H
#     ${SRC_DIR}/React.cpp
#     ${SRC_DIR}/Riemann.H
#     ${SRC_DIR}/Setup.cpp
#     ${SRC_DIR}/Sources.cpp
#     ${SRC_DIR}/SparseData.H
#     ${SRC_DIR}/SumIQ.cpp
#     ${SRC_DIR}/SumUtils.cpp
#     ${SRC_DIR}/Tagging.H
#     ${SRC_DIR}/Tagging.cpp
#     ${SRC_DIR}/Timestep.H
#     ${SRC_DIR}/Utilities.H
#     ${SRC_DIR}/Utilities.cpp
#     ${SRC_DIR}/WENO.H
#)

include(AMReXBuildInfo)
generate_buildinfo(${pelec_lib_name} ${CMAKE_SOURCE_DIR})
target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

if(PELEC_ENABLE_MPI)
  target_link_libraries(${pelec_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
endif()

#PeleC include directories
target_include_directories(${pelec_lib_name} PUBLIC ${SRC_DIR})
target_include_directories(${pelec_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

#Link to amrex libraries
target_link_libraries(${pelec_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${pelec_lib_name} PUBLIC AMReX::amrex)

if(PELEC_ENABLE_CUDA)
  set(pctargets "${pelec_lib_name}")
  foreach(tgt IN LISTS pctargets)
    get_target_property(PELEC_SOURCES ${tgt} SOURCES)
    list(FILTER PELEC_SOURCES INCLUDE REGEX "\\.cpp")
    set_source_files_properties(${PELEC_SOURCES} PROPERTIES LANGUAGE CUDA)
  endforeach()
  set_target_properties(${pelec_lib_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_compile_options(${pelec_lib_name} PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
endif()

#Define what we want to be installed during a make install 
#install(TARGETS ${pelec_lib_name}
#        RUNTIME DESTINATION bin
#        ARCHIVE DESTINATION lib
#        LIBRARY DESTINATION lib)
