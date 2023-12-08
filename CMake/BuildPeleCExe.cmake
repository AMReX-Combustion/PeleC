function(build_pelec_exe pelec_exe_name pelec_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${pelec_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pelec_exe_name} PROPERTIES CXX_CLANG_TIDY
                          "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
  endif()

  target_sources(${pelec_exe_name}
     PRIVATE
       prob_parm.H
       prob.H
       prob.cpp
  )

  #PeleC include directories
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_BINARY_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/Params/param_includes)

  # Spray
  set(PELE_PHYSICS_SPRAY_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Source/Spray)
  if(PELEC_ENABLE_AMREX_PARTICLES AND PELE_SPRAY_FUEL_NUM GREATER 0)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_SPRAY)
    target_compile_definitions(${pelec_exe_name} PRIVATE SPRAY_FUEL_NUM=${PELE_SPRAY_FUEL_NUM})
    target_sources(${pelec_exe_name} PRIVATE
                   SprayParticlesInitInsert.cpp
                   ${SRC_DIR}/Particle.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/Drag.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayDerive.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayFuelData.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayIO.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayInjection.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayInterpolation.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.H
                   ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/SpraySB.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/SpraySetup.cpp
                   ${PELE_PHYSICS_SPRAY_DIR}/WallFunctions.H
                   ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/AhamedSplash.H
                   ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/ReitzKHRT.H
                   ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/SBData.H
                   ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/TABBreakup.H
                   ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/WallFilm.H
                   ${PELE_PHYSICS_SPRAY_DIR}/Distribution/DistBase.H
                   ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.H
                   ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.cpp)
    target_include_directories(${pelec_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR})
    target_include_directories(${pelec_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/Distribution)
    target_include_directories(${pelec_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash)
  endif()

  target_sources(${pelec_exe_name}
     PRIVATE
       ${SRC_DIR}/Advance.cpp
       ${SRC_DIR}/BCfill.cpp
       ${SRC_DIR}/Bld.cpp
       ${SRC_DIR}/Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/Diffterm.H
       ${SRC_DIR}/Diffusion.H
       ${SRC_DIR}/Diffusion.cpp
       ${SRC_DIR}/EB.H
       ${SRC_DIR}/EB.cpp
       ${SRC_DIR}/EBStencilTypes.H
       ${SRC_DIR}/External.cpp
       ${SRC_DIR}/Forcing.cpp
       ${SRC_DIR}/GradUtil.H
       ${SRC_DIR}/Hydro.H
       ${SRC_DIR}/Hydro.cpp
       ${SRC_DIR}/Geometry.H
       ${SRC_DIR}/Geometry.cpp
       ${SRC_DIR}/Godunov.H
       ${SRC_DIR}/Godunov.cpp
       ${SRC_DIR}/PLM.H
       ${SRC_DIR}/PPM.H
       ${SRC_DIR}/PPM.cpp
       ${SRC_DIR}/InitEB.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/IO.H
       ${SRC_DIR}/IO.cpp
       ${SRC_DIR}/LES.H
       ${SRC_DIR}/LES.cpp
       ${SRC_DIR}/MOL.H
       ${SRC_DIR}/MOL.cpp
       ${SRC_DIR}/PeleC.H
       ${SRC_DIR}/PeleC.cpp
       ${SRC_DIR}/PeleCAmr.H
       ${SRC_DIR}/PeleCAmr.cpp
       ${SRC_DIR}/ProblemDerive.H
       ${SRC_DIR}/React.cpp
       ${SRC_DIR}/Riemann.H
       ${SRC_DIR}/Setup.cpp
       ${SRC_DIR}/Sources.cpp
       ${SRC_DIR}/SparseData.H
       ${SRC_DIR}/SumIQ.cpp
       ${SRC_DIR}/SumUtils.cpp
       ${SRC_DIR}/Tagging.H
       ${SRC_DIR}/Tagging.cpp
       ${SRC_DIR}/Timestep.H
       ${SRC_DIR}/Utilities.H
       ${SRC_DIR}/Utilities.cpp
       ${SRC_DIR}/WENO.H
  )

  # Soot
  if(PELEC_ENABLE_SOOT)
      target_sources(${pelec_exe_name} PRIVATE
                   ${SRC_DIR}/Soot.cpp)
  endif()

  if(PELEC_ENABLE_ASCENT)
    target_sources(${pelec_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleAscent.H
        ${SRC_DIR}/PeleAscent.cpp
    )
  endif()

  if(PELEC_ENABLE_MASA)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_MASA)
    target_sources(${pelec_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
    target_link_libraries(${pelec_exe_name} PRIVATE MASA::MASA)
    if(PELEC_ENABLE_FPE_TRAP_FOR_TESTS)
      set_source_files_properties(${SRC_DIR}/PeleC.cpp PROPERTIES COMPILE_DEFINITIONS PELEC_ENABLE_FPE_TRAP)
    endif()
  endif()

  if(NOT "${pelec_exe_name}" STREQUAL "PeleC-UnitTests")
    target_sources(${pelec_exe_name}
       PRIVATE
         ${CMAKE_SOURCE_DIR}/Source/main.cpp
    )
  endif()

  if(PELEC_ENABLE_CUDA)
    set(pctargets "${pelec_exe_name};${pelec_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(PELEC_SOURCES ${tgt} SOURCES)
      list(FILTER PELEC_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELEC_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pelec_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()

  target_link_libraries(${pelec_exe_name} PRIVATE ${pelec_lib_name} AMReX::amrex)

  #Define what we want to be installed during a make install
  install(TARGETS ${pelec_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
