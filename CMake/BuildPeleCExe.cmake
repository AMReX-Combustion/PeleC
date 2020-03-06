function(build_pelec_exe pelec_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/SourceCpp)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/SourceCpp)

  include(${CMAKE_SOURCE_DIR}/CMake/SetCompileFlags.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/SetRpath.cmake)

  #Gather all other source files  
  add_subdirectory(${SRC_DIR}/Chemistry ${BIN_DIR}/Chemistry)
  add_subdirectory(${SRC_DIR}/Params ${BIN_DIR}/Params)
  
  if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
    add_subdirectory(${SRC_DIR}/EOS/Fuego ${BIN_DIR}/EOS/Fuego)
  else()
    add_subdirectory(${SRC_DIR}/EOS/Gamma ${BIN_DIR}/EOS/Gamma)
  endif()
  
  if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
    add_subdirectory(${SRC_DIR}/Transport/Simple ${BIN_DIR}/Transport/Simple)
  elseif("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
    add_subdirectory(${SRC_DIR}/Transport/Constant ${BIN_DIR}/Transport/Constant)
  endif()
  
  if(PELEC_ENABLE_REACTIONS)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_REACTIONS)
    target_sources(${pelec_exe_name} PRIVATE ${SRC_DIR}/React.H ${SRC_DIR}/React.cpp)
    target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Evaluation)
  endif()
  
  if(PELEC_ENABLE_EXPLICIT_REACT)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_EXPLICIT_REACT)
  endif()
  
  if(PELEC_ENABLE_MASA)
    target_sources(${pelec_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
    target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_MASA)
  endif()
  
  if(PELEC_ENABLE_EB)
    target_sources(${pelec_exe_name} PRIVATE ${SRC_DIR}/EB.H ${SRC_DIR}/EB.cpp ${SRC_DIR}/InitEB.cpp ${SRC_DIR}/SparseData.H ${SRC_DIR}/EBStencilTypes.H)
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
       ${SRC_DIR}/IO.H
       ${SRC_DIR}/IO.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/LES.H
       ${SRC_DIR}/LES.cpp
       ${SRC_DIR}/MOL.H
       ${SRC_DIR}/MOL.cpp
       ${SRC_DIR}/Particle.cpp
       ${SRC_DIR}/PeleC.H
       ${SRC_DIR}/PeleC.cpp
       ${SRC_DIR}/PLM.H
       ${SRC_DIR}/PLM.cpp
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
       ${SRC_DIR}/Timestep.cpp
       ${SRC_DIR}/Utilities.H
       ${SRC_DIR}/Utilities.cpp
       ${SRC_DIR}/main.cpp
     )
  
  #Add generated source files
  set_property(SOURCE ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp PROPERTY GENERATED 1)
  target_sources(${pelec_exe_name}
     PRIVATE
        ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
  )

  if(PELEC_ENABLE_MASA)
   if(MASA_FOUND)
     #Link our executable to the MASA libraries, etc
     target_link_libraries(${pelec_exe_name} PRIVATE ${MASA_LIBRARY})
     target_compile_definitions(${pelec_exe_name} PRIVATE USE_MASA DO_PROBLEM_POST_TIMESTEP DO_PROBLEM_POST_INIT)
     target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_INCLUDE_DIRS})
     target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${MASA_MOD_DIRS})
   endif()
  endif()

  if(PELEC_ENABLE_MPI)
    target_link_libraries(${pelec_exe_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #PeleC include directories
  target_include_directories(${pelec_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_BINARY_DIR})
  
  #Needed for AMReX_buildInfo.H
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/C_scripts)
  
  #Link to amrex library
  target_link_libraries(${pelec_exe_name} PRIVATE amrex)

  #Set the dependencies on targets so the generated source code files are there before we try to build the executable 
  add_dependencies(${pelec_exe_name} generate_build_info)
 
  #Define what we want to be installed during a make install 
  install(TARGETS ${pelec_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
