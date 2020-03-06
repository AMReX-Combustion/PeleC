if("${PELEC_EOS_MODEL}" STREQUAL "Fuego")
  target_sources(${pelec_exe_name} PRIVATE EOS_FUEGO.H EOS_FUEGO.cpp) 
else()
  target_sources(${pelec_exe_name} PRIVATE EOS_GAMMA.H EOS_GAMMA.cpp) 
endif()

#Create the full path to the extra case-specific source files
#Each PELEC_EXTRA_SOURCE must be an explicit path to each source file at the moment
include(${CMAKE_BINARY_DIR}/extra_sources.cmake)
foreach(PELEC_EXTRA_SOURCE ${PELEC_EXTRA_SOURCES})
  list(APPEND MY_EXTRA_SOURCES ${PELEC_EXTRA_SOURCE})
endforeach()

if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
  target_compile_definitions(${pelec_exe_name} PRIVATE SIMPLE_TRANSPORT)
elseif("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
  target_compile_definitions(${pelec_exe_name} PRIVATE CONSTANT_TRANSPORT)
endif()

if(PELEC_ENABLE_REACTIONS)
  target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_REACTIONS)
  target_sources(${pelec_exe_name} PRIVATE React.H React.cpp)
  target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Evaluation)
endif()

if(PELEC_ENABLE_EXPLICIT_REACT)
  target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_EXPLICIT_REACT)
endif()

if(PELEC_ENABLE_MASA)
  target_sources(${pelec_exe_name} PRIVATE MMS.cpp)
  target_compile_definitions(${pelec_exe_name} PRIVATE PELEC_USE_MASA)
endif()

if(PELEC_ENABLE_EB)
  target_sources(${pelec_exe_name} PRIVATE EB.H EB.cpp InitEB.cpp SparseData.H EBStencilTypes.H)
endif()

if("${PELEC_TRANSPORT_MODEL}" STREQUAL "Simple")
  target_sources(${pelec_exe_name} PRIVATE SimpleTransport.H SimpleTransportParams.H SimpleTransportParams.cpp)
elseif("${PELEC_TRANSPORT_MODEL}" STREQUAL "Constant")
  target_sources(${pelec_exe_name} PRIVATE ConstantTransport.H ConstantTransportParams.H ConstantTransportParams.cpp)
endif()

target_sources(${pelec_exe_name}
   PRIVATE
      Advance.cpp
      BCfill.cpp
      Bld.cpp
      Constants.H
      Derive.H
      Derive.cpp
      Diffterm.H
      Diffterm.cpp
      Diffusion.H
      Diffusion.cpp
      EOS.H
      EOS.cpp
      External.cpp
      Filter.H
      Filter.cpp
      Forcing.H
      Forcing.cpp
      GradUtil.H
      GradUtil.cpp
      Hydro.H
      Hydro.cpp
      IO.H
      IO.cpp
      IndexDefines.H
      LES.H
      LES.cpp
      MOL.H
      MOL.cpp
      Particle.cpp
      PeleC.H
      PeleC.cpp
      PLM.H
      PLM.cpp
      Problem.H
      ProblemDerive.H
      Riemann.H
      Setup.cpp
      Sources.cpp
      SumIQ.cpp
      SumUtils.cpp
      Tagging.H
      Tagging.cpp
      Timestep.H
      Timestep.cpp
      Transport.H
      Transport.cpp
      Utilities.H
      Utilities.cpp
      main.cpp
      ${MY_EXTRA_SOURCES}
   )

#Add generated source files
set_property(SOURCE ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp PROPERTY GENERATED 1)
target_sources(${pelec_exe_name}
   PRIVATE
      ${GENERATED_FILES_DIR}/AMReX_buildInfo.cpp
)

#PeleC include directories
target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

#Needed for AMReX_buildInfo.H
target_include_directories(${pelec_exe_name} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/C_scripts)

#Gather all other source files  
add_subdirectory(Chemistry)
add_subdirectory(Params)

#Link to amrex library
target_link_libraries(${pelec_exe_name} PRIVATE amrex)
