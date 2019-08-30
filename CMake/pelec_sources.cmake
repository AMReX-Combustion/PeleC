#Aggregate PeleC source files
function(get_pelec_sources pelec_exe_name)

  if(${PELEC_DIM} EQUAL 1)
    set(PELEC_ENABLE_EB OFF)
  endif()
  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Source")
  add_sources(GlobalSourceList
     ${PELEC_SOURCE_DIR}/PeleC.cpp
     ${PELEC_SOURCE_DIR}/PeleCBld.cpp
     ${PELEC_SOURCE_DIR}/PeleC_MOL.cpp
     ${PELEC_SOURCE_DIR}/PeleC_advance.cpp
     ${PELEC_SOURCE_DIR}/PeleC_bcfill.cpp
     ${PELEC_SOURCE_DIR}/PeleC_external.cpp
     ${PELEC_SOURCE_DIR}/PeleC_forcing.cpp
     ${PELEC_SOURCE_DIR}/PeleC_hydro.cpp
     ${PELEC_SOURCE_DIR}/PeleC_init_eb.cpp
     ${PELEC_SOURCE_DIR}/PeleC_io.cpp
     ${PELEC_SOURCE_DIR}/PeleC_setup.cpp
     ${PELEC_SOURCE_DIR}/PeleC_sources.cpp
     ${PELEC_SOURCE_DIR}/Prob.cpp
     ${PELEC_SOURCE_DIR}/main.cpp
     ${PELEC_SOURCE_DIR}/sum_integrated_quantities.cpp
     ${PELEC_SOURCE_DIR}/sum_utils.cpp
     ${PELEC_SOURCE_DIR}/Filter.cpp
     ${PELEC_SOURCE_DIR}/PeleC_les.cpp
     #${PELEC_SOURCE_DIR}/PeleC_error.cpp
  )
  if(PELEC_ENABLE_REACTIONS)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/PeleC_react.cpp
     )
  endif()
  if(PELEC_ENABLE_PARTICLES)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/PeleCParticles.cpp
     )
  endif()
  if(PELEC_ENABLE_MASA)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/PeleC_mms.cpp
     )
  endif()
  if(PELEC_ENABLE_EB)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/EBStencilTypes_mod.F90
     )
  endif()

  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Source/Src_${PELEC_DIM}d")
  add_sources(GlobalSourceList
     ${PELEC_SOURCE_DIR}/advection_util_${PELEC_DIM}d.f90
     ${PELEC_SOURCE_DIR}/impose_NSCBC_${PELEC_DIM}d.f90
     ${PELEC_SOURCE_DIR}/riemann_${PELEC_DIM}d.F90
     ${PELEC_SOURCE_DIR}/set_bc_mask_${PELEC_DIM}d.f90
     ${PELEC_SOURCE_DIR}/filter_${PELEC_DIM}d.f90
     ${PELEC_SOURCE_DIR}/lesterm_${PELEC_DIM}d.f90
     #${PELEC_SOURCE_DIR}/slope_mol_${PELEC_DIM}d.f90
     #${PELEC_SOURCE_DIR}/PeleC_mol_${PELEC_DIM}d.F90
  )
  if(${PELEC_DIM} GREATER 1)
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/grad_utils_${PELEC_DIM}d.F90
     )
  endif()
  if(${PELEC_DIM} GREATER 1 AND PELEC_ENABLE_EB)
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/PeleC_init_eb_${PELEC_DIM}d.F90
     )
  endif()
  if("${PELEC_TRANSPORT_TYPE}" STREQUAL "IDEAL_GAS")
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/diffterm_${PELEC_DIM}d.f90
     )
  elseif("${PELEC_TRANSPORT_TYPE}" STREQUAL "REAL_GAS")
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/diffterm_nonideal_${PELEC_DIM}d.f90
     )
  endif()
  if(PELEC_ENABLE_MOL AND PELEC_ENABLE_EB)
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/Hyp_pele_MOL_${PELEC_DIM}d.F90
       ${PELEC_SOURCE_DIR}/slope_mol_${PELEC_DIM}d_EB.f90
     )
  else()
     add_sources(GlobalSourceList
       ${PELEC_SOURCE_DIR}/trace_${PELEC_DIM}d.f90
       ${PELEC_SOURCE_DIR}/trace_ppm_${PELEC_DIM}d.f90
       ${PELEC_SOURCE_DIR}/ppm_${PELEC_DIM}d.f90
       ${PELEC_SOURCE_DIR}/slope_${PELEC_DIM}d.f90
       ${PELEC_SOURCE_DIR}/PeleC_advection_${PELEC_DIM}d.F90
       ${PELEC_SOURCE_DIR}/PeleC_${PELEC_DIM}d.F90
     )
     if(${PELEC_DIM} GREATER 1)
       add_sources(GlobalSourceList
         ${PELEC_SOURCE_DIR}/trans_${PELEC_DIM}d.F90
       )
     endif()
  endif()

  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Source/Src_nd")
  add_sources(GlobalSourceList
     ${PELEC_SOURCE_DIR}/Derive_nd.F90
     ${PELEC_SOURCE_DIR}/Diffusion_nd.f90
     ${PELEC_SOURCE_DIR}/Make.package
     ${PELEC_SOURCE_DIR}/PeleC_nd.F90
     ${PELEC_SOURCE_DIR}/PeleC_util.F90
     ${PELEC_SOURCE_DIR}/Problem.f90
     ${PELEC_SOURCE_DIR}/Tagging_nd.f90
     ${PELEC_SOURCE_DIR}/advection_util_nd.F90
     ${PELEC_SOURCE_DIR}/amrinfo.f90
     ${PELEC_SOURCE_DIR}/ext_src_nd.f90
     ${PELEC_SOURCE_DIR}/filcc_nd.F90
     ${PELEC_SOURCE_DIR}/flatten_nd.F90
     ${PELEC_SOURCE_DIR}/forcing_src_nd.F90
     ${PELEC_SOURCE_DIR}/interpolate.f90
     ${PELEC_SOURCE_DIR}/io.f90
     ${PELEC_SOURCE_DIR}/math.f90
     ${PELEC_SOURCE_DIR}/meth_params.F90
     ${PELEC_SOURCE_DIR}/parmparse_fi.cpp
     ${PELEC_SOURCE_DIR}/parmparse_mod.f90
     ${PELEC_SOURCE_DIR}/prob_params.f90
     ${PELEC_SOURCE_DIR}/problem_derive_nd.F90
     ${PELEC_SOURCE_DIR}/problem_tagging_nd.F90
     ${PELEC_SOURCE_DIR}/rk_params.f90
     ${PELEC_SOURCE_DIR}/riemann_util.f90
     ${PELEC_SOURCE_DIR}/string_mod.f90
     ${PELEC_SOURCE_DIR}/sums_nd.f90
     ${PELEC_SOURCE_DIR}/timestep.F90
     ${PELEC_SOURCE_DIR}/weno.f90
     #${PELEC_SOURCE_DIR}/bc_fill_nd.F90
     #${PELEC_SOURCE_DIR}/Prob_nd.F90
  )
  if(PELEC_ENABLE_REACTIONS)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/React_nd.F90
     )
  endif()
  if(PELEC_ENABLE_MASA)
     add_sources(GlobalSourceList
        ${PELEC_SOURCE_DIR}/mms_src_nd.F90
     )
  endif()

  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/constants")
  add_sources(GlobalSourceList
     ${PELEC_SOURCE_DIR}/constants_cgs.f90
  )

  #Add generated source files
  add_sources(GlobalSourceList
     ${CMAKE_BINARY_DIR}/generated_files/${pelec_exe_name}_generated_files/extern.f90
     ${CMAKE_BINARY_DIR}/generated_files/${pelec_exe_name}_generated_files/AMReX_buildInfo.cpp
  )
endfunction(get_pelec_sources pelec_exe_name)
