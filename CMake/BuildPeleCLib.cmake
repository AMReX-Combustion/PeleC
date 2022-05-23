function(build_pelec_lib pelec_lib_name)
  if (NOT (TARGET ${pelec_lib_name}))
    add_library(${pelec_lib_name} OBJECT)

    set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
    set(PELE_PHYSICS_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport")
    set(PELE_PHYSICS_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos")
    set(PELE_PHYSICS_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Mechanism/Models/${PELEC_CHEMISTRY_MODEL}")
    set(AMREX_SUNDIALS_DIR ${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS)

    if(CLANG_TIDY_EXE)
      set_target_properties(${pelec_lib_name} PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
    endif()

    include(SetPeleCCompileFlags)
    
    target_sources(${pelec_lib_name}
      PRIVATE
      ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.cpp
      ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow/turbinflow.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/TurbInflow)
    
    target_sources(${pelec_lib_name}
      PRIVATE
      ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.cpp
      ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManager.H
      ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager/PltFileManagerBCFill.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Utility/PltFileManager)
    
    target_sources(${pelec_lib_name} PRIVATE ${AMREX_SUNDIALS_DIR}/AMReX_Sundials.H
                                             ${AMREX_SUNDIALS_DIR}/AMReX_Sundials.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.H
                                             ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.H)
    target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${AMREX_SUNDIALS_DIR})

    target_include_directories(${pelec_lib_name} PUBLIC "${PELE_PHYSICS_SRC_DIR}/Source")

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.cpp
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportParams.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportTypes.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Constant.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Simple.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Sutherland.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_TRANSPORT_DIR})
    if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_CONSTANT_TRANSPORT)
    endif()
    if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_SIMPLE_TRANSPORT)
    endif()
    if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Sutherland")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_SUTHERLAND_TRANSPORT)
    endif()

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_EOS_DIR}/EOS.cpp
                   ${PELE_PHYSICS_EOS_DIR}/EOS.H
                   ${PELE_PHYSICS_EOS_DIR}/GammaLaw.H
                   ${PELE_PHYSICS_EOS_DIR}/Fuego.H
                   ${PELE_PHYSICS_EOS_DIR}/SRK.H)
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_EOS_DIR})
    if("${PELEC_EOS_MODEL}" STREQUAL "GammaLaw")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_GAMMALAW_EOS)
    endif()
    if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_FUEGO_EOS)
    endif()
    if("${PELEC_EOS_MODEL}" STREQUAL "Soave-Redlich-Kwong")
      target_compile_definitions(${pelec_lib_name} PUBLIC USE_SRK_EOS)
    endif()

    target_sources(${pelec_lib_name} PRIVATE
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.cpp
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.H)
    target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${PELE_PHYSICS_MECHANISM_DIR})

    target_sources(${pelec_lib_name}
      PRIVATE
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorArkode.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorBase.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvode.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeCustomLinSolver.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeJacobian.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodePreconditioner.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorCvodeUtils.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorNull.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorRK64.cpp
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorTypes.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.H
        ${PELE_PHYSICS_SRC_DIR}/Reactions/ReactorUtils.cpp
    )
    target_include_directories(${pelec_lib_name} PUBLIC ${PELE_PHYSICS_SRC_DIR}/Reactions)

    include(AMReXBuildInfo)
    generate_buildinfo(${pelec_lib_name} ${CMAKE_SOURCE_DIR})
    target_include_directories(${pelec_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)
    
    target_link_libraries(${pelec_lib_name} PUBLIC sundials_arkode sundials_cvode)
    
    if(PELEC_ENABLE_CUDA)
      target_link_libraries(${pelec_lib_name} PUBLIC sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
    elseif(PELEC_ENABLE_HIP)
      target_link_libraries(${pelec_lib_name} PUBLIC sundials_nvechip)
    elseif(PELEC_ENABLE_DPCPP)
      target_link_libraries(${pelec_lib_name} PUBLIC sundials_nvecsycl)
    endif()
    
    if(PELEC_ENABLE_MPI)
      target_link_libraries(${pelec_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    endif()
    
    #Link to amrex libraries
    target_link_libraries(${pelec_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api AMReX::amrex)

  endif()
endfunction()
