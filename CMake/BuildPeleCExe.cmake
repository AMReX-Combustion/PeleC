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

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/SourceCpp)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/SourceCpp/${pelec_exe_name})

  include(SetPeleCCompileFlags)

  add_subdirectory(${SRC_DIR}/Params ${BIN_DIR}/Params/${pelec_exe_name})

  set(PELEC_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Transport/${PELEC_TRANSPORT_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_TRANSPORT_DIR}/Transport.H
                 ${PELEC_TRANSPORT_DIR}/Transport.cpp
                 ${PELEC_TRANSPORT_DIR}/TransportParams.H)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_TRANSPORT_DIR})

  set(PELEC_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Eos/${PELEC_EOS_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_EOS_DIR}/EOS.cpp
                 ${PELEC_EOS_DIR}/EOS.H)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_EOS_DIR})

  set(PELEC_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Mechanism/Models/${PELEC_CHEMISTRY_MODEL}")
  target_sources(${pelec_exe_name} PRIVATE
                 ${PELEC_MECHANISM_DIR}/chemistry_file.H
                 ${PELEC_MECHANISM_DIR}/mechanism.cpp
                 ${PELEC_MECHANISM_DIR}/mechanism.h)
  # Avoid warnings from mechanism.cpp for now
  if(NOT PELEC_ENABLE_CUDA)
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
      if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
        list(APPEND MY_CXX_FLAGS "-Wno-unreachable-code"
                                 "-Wno-null-dereference"
                                 "-Wno-float-conversion"
                                 "-Wno-shadow"
                                 "-Wno-overloaded-virtual")
      endif()
      if(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-variable")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-parameter")
        list(APPEND MY_CXX_FLAGS "-Wno-vla-extension")
        list(APPEND MY_CXX_FLAGS "-Wno-zero-length-array")
      elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-variable")
        list(APPEND MY_CXX_FLAGS "-Wno-unused-parameter")
        list(APPEND MY_CXX_FLAGS "-Wno-vla")
        list(APPEND MY_CXX_FLAGS "-Wno-pedantic")
      endif()
    endif()
  endif()
  separate_arguments(MY_CXX_FLAGS)
  set_source_files_properties(${PELEC_MECHANISM_DIR}/mechanism.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  set_source_files_properties(${PELEC_MECHANISM_DIR}/chemistry_file.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  set_source_files_properties(${PELEC_MECHANISM_DIR}/mechanism.H PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELEC_MECHANISM_DIR})
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${PELE_PHYSICS_SRC_DIR}/Support/Fuego/Evaluation)
  
  if(PELEC_ENABLE_REACTIONS)
    if(PELEC_ENABLE_SUNDIALS)
      if(PELEC_ENABLE_CUDA)
        set(DEVICE GPU)
      else()
        set(DEVICE CPU)
      endif()
      target_compile_definitions(${pelec_exe_name} PRIVATE USE_SUNDIALS_PP)
      target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/arkode/reactor.cpp
                                               ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/arkode/reactor.h
                                               ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/AMREX_misc.H)
      set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/arkode/reactor.cpp PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      set_source_files_properties(${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/arkode/reactor.h PROPERTIES COMPILE_OPTIONS "${MY_CXX_FLAGS}")
      target_link_libraries(${pelec_exe_name} PRIVATE sundials_arkode)
      target_include_directories(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE})
      target_include_directories(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/arkode)
      if(PELEC_ENABLE_CUDA)
        target_sources(${pelec_exe_name} PRIVATE ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/AMReX_SUNMemory.cpp
                                                 ${PELE_PHYSICS_SRC_DIR}/Reactions/Fuego/${DEVICE}/AMReX_SUNMemory.H)
        target_link_libraries(${pelec_exe_name} PRIVATE sundials_nveccuda)
      endif()
    endif()
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_REACTIONS)
    target_sources(${pelec_exe_name} PRIVATE
                   ${SRC_DIR}/React.H
                   ${SRC_DIR}/React.cpp)
  endif()

  if(PELEC_ENABLE_FORCING)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_FORCING)
  endif()
  
  if(PELEC_ENABLE_EB)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_EB)
    target_sources(${pelec_exe_name}
                   PRIVATE
                   ${SRC_DIR}/EB.H
                   ${SRC_DIR}/EB.cpp
                   ${SRC_DIR}/InitEB.cpp
                   ${SRC_DIR}/SparseData.H
                   ${SRC_DIR}/EBStencilTypes.H)
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
