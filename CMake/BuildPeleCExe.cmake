function(build_pelec_exe pelec_exe_name)

  add_executable(${pelec_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pelec_exe_name} PROPERTIES CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
  endif()

  target_sources(${pelec_exe_name}
     PRIVATE
       prob_parm.H
       prob.H
       prob.cpp
  )
  
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  set(PELE_PHYSICS_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics)
  set(PELE_PHYSICS_BIN_DIR ${CMAKE_BINARY_DIR}/Submodules/PelePhysics/${pelec_exe_name})

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${pelec_exe_name})

  include(SetPeleCCompileFlags)

  add_subdirectory(${SRC_DIR}/Params ${BIN_DIR}/Params/${pelec_exe_name})

  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE "${PELE_PHYSICS_SRC_DIR}/Source")

  set(PELEC_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport/${PELEC_TRANSPORT_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_TRANSPORT_DIR}/Transport.H
                 ${PELEC_TRANSPORT_DIR}/Transport.cpp
                 ${PELEC_TRANSPORT_DIR}/TransportParams.H)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_TRANSPORT_DIR})
  if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_SIMPLE)
  endif()

  set(PELEC_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos/${PELEC_EOS_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_EOS_DIR}/EOS.cpp
                 ${PELEC_EOS_DIR}/EOS.H)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_EOS_DIR})
  if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_FUEGO)
  endif()
  if("${PELEC_EOS_MODEL}" STREQUAL "Soave-Redlich-Kwong")
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_SRK)
  endif()

  set(PELEC_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models/${PELEC_CHEMISTRY_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_MECHANISM_DIR}/mechanism.cpp
                 ${PELEC_MECHANISM_DIR}/mechanism.H)
  # Avoid warnings from certain files
  if((NOT PELEC_ENABLE_CUDA) AND (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$"))
    list(APPEND MY_CXX_FLAGS "-w")
  endif()
  separate_arguments(MY_CXX_FLAGS)
  set_source_files_properties(${PELEC_MECHANISM_DIR}/mechanism.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  set_source_files_properties(${PELEC_MECHANISM_DIR}/mechanism.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_MECHANISM_DIR})
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Evaluation)

  if(PELEC_ENABLE_EB)
    set_source_files_properties(${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_create_itracker_${PELEC_DIM}d.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
    set_source_files_properties(${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_redistribution.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
    set_source_files_properties(${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_redistribution.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
    set_source_files_properties(${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_state_redistribute.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  endif()
  
  if(PELEC_ENABLE_REACTIONS)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_REACTIONS)
    target_sources(${pelec_exe_name} PRIVATE
                   ${SRC_DIR}/React.H
                   ${SRC_DIR}/React.cpp)
    if(PELEC_ENABLE_SUNDIALS)
      target_compile_definitions(${pelec_exe_name} PRIVATE USE_SUNDIALS_PP)
      target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/reactor.cpp
                                               ${PELE_PHYSICS_SRC_DIR}/Reactions/reactor.H
                                               ${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_utils.H)
      if("${PELEC_SUNDIALS_INTEGRATOR}" STREQUAL "arkode")
          target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}.cpp)
          set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      elseif("${PELEC_SUNDIALS_INTEGRATOR}" STREQUAL "cvode")
        target_compile_definitions(${pelec_exe_name} PRIVATE COMPILE_JACOBIAN)
        if(PELEC_ENABLE_CUDA)
          target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}_GPU.cpp)
          set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}_GPU.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
        else()
          target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}_CPU.cpp)
          set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_${PELEC_SUNDIALS_INTEGRATOR}_CPU.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
        endif()
      endif()
      set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/reactor.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/reactor.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR}/reactor_utils.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      target_include_directories(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions)
      target_include_directories(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/${PELEC_SUNDIALS_INTEGRATOR})
      target_link_libraries(${pelec_exe_name} PRIVATE sundials_${PELEC_SUNDIALS_INTEGRATOR})
      if(PELEC_ENABLE_CUDA)
        target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/AMReX_SUNMemory.cpp
                                                 ${PELE_PHYSICS_SRC_DIR}/Reactions/AMReX_SUNMemory.H)
        target_link_libraries(${pelec_exe_name} PRIVATE sundials_nveccuda)
        if(PELEC_SUNDIALS_INTEGRATOR STREQUAL "cvode")
          target_link_libraries(${pelec_exe_name} PRIVATE sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
        endif()
      endif()
    endif()
  endif()

  if(PELEC_ENABLE_FORCING)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_FORCING)
  endif()
  
  if(PELEC_ENABLE_EB)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_EB)
    target_sources(${pelec_exe_name}
                   PRIVATE
                   ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_create_itracker_${PELEC_DIM}d.cpp
                   ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_redistribution.H
                   ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_redistribution.cpp
                   ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_state_utils.cpp
                   ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution/hydro_state_redistribute.cpp
                   ${SRC_DIR}/EB.H
                   ${SRC_DIR}/EB.cpp
                   ${SRC_DIR}/InitEB.cpp
                   ${SRC_DIR}/SparseData.H
                   ${SRC_DIR}/EBStencilTypes.H)
     target_include_directories(${pelec_exe_name} PRIVATE ${AMREX_HYDRO_SUBMOD_LOCATION}/Redistribution)
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
       ${SRC_DIR}/Diffterm.cpp
       ${SRC_DIR}/Diffusion.H
       ${SRC_DIR}/Diffusion.cpp
       ${SRC_DIR}/External.cpp
       ${SRC_DIR}/Filter.H
       ${SRC_DIR}/Filter.cpp
       ${SRC_DIR}/Forcing.H
       ${SRC_DIR}/Forcing.cpp
       ${SRC_DIR}/GradUtil.H
       ${SRC_DIR}/GradUtil.cpp
       ${SRC_DIR}/Hydro.H
       ${SRC_DIR}/Hydro.cpp
       ${SRC_DIR}/Godunov.H
       ${SRC_DIR}/Godunov.cpp
       ${SRC_DIR}/PLM.H
       ${SRC_DIR}/PPM.H
       ${SRC_DIR}/PPM.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/IndexDefines.cpp
       ${SRC_DIR}/IO.H
       ${SRC_DIR}/IO.cpp
       ${SRC_DIR}/LES.H
       ${SRC_DIR}/LES.cpp
       ${SRC_DIR}/MOL.H
       ${SRC_DIR}/MOL.cpp
       ${SRC_DIR}/Particle.cpp
       ${SRC_DIR}/PeleC.H
       ${SRC_DIR}/PeleC.cpp
       ${SRC_DIR}/Problem.H
       ${SRC_DIR}/ProblemDerive.H
       ${SRC_DIR}/Riemann.H
       ${SRC_DIR}/Setup.cpp
       ${SRC_DIR}/Sources.cpp
       ${SRC_DIR}/SumIQ.cpp
       ${SRC_DIR}/SumUtils.cpp
       ${SRC_DIR}/Tagging.H
       ${SRC_DIR}/Tagging.cpp
       ${SRC_DIR}/Timestep.H
       ${SRC_DIR}/Utilities.H
       ${SRC_DIR}/Utilities.cpp
       ${SRC_DIR}/WENO.H
  )

  if(NOT "${pelec_exe_name}" STREQUAL "PeleC-UnitTests")
    target_sources(${pelec_exe_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()
  
  include(AMReXBuildInfo)
  generate_buildinfo(${pelec_exe_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(PELEC_ENABLE_MASA)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_MASA)
    target_sources(${pelec_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
    target_link_libraries(${pelec_exe_name} PRIVATE MASA::MASA)
    if(PELEC_ENABLE_FPE_TRAP_FOR_TESTS)
      set_source_files_properties(${SRC_DIR}/PeleC.cpp PROPERTIES COMPILE_DEFINITIONS PELEC_ENABLE_FPE_TRAP)
    endif()
  endif()

  if(PELEC_ENABLE_MPI)
    target_link_libraries(${pelec_exe_name} PRIVATE $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #PeleC include directories
  target_include_directories(${pelec_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

  #Link to amrex library
  target_link_libraries(${pelec_exe_name} PRIVATE AMReX::amrex)

  if(PELEC_ENABLE_CUDA)
    set(pctargets "${pelec_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(PELEC_SOURCES ${tgt} SOURCES)
      list(FILTER PELEC_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELEC_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pelec_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pelec_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${pelec_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
