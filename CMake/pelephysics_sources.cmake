function(get_pelephysics_sources)

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Eos")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/eos_type.f90
   )

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Eos/${PELEC_EOS_MODEL}")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/eos.F90
   )

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Evaluation")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/LinAlg.f
      ${PELEPHYSICS_SOURCE_DIR}/bdf.f90
      ${PELEPHYSICS_SOURCE_DIR}/bdf_data.f90
      ${PELEPHYSICS_SOURCE_DIR}/math_d.f
      ${PELEPHYSICS_SOURCE_DIR}/vode.f
      ${PELEPHYSICS_SOURCE_DIR}/vode_module.f90
   )
   if(USE_FUEGO)
     add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/egz_module.f90
     )
   endif()

   # CMake has the worst possible ways to check whether a string is defined or not
   if(NOT DEFINED PELEC_CHEMISTRY_MODEL)
     set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism")
     add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/chemistry_module_null.f90
     )
   else()
     set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism")
     add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/chemistry_module.F90
     )
     set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism/Models")
     add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/mod_fuego.f90
     )
     set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Support/Fuego/Mechanism/Models/${PELEC_CHEMISTRY_MODEL}")
     add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/${PELEC_CHEMISTRY_MODEL}.cpp
     )
   endif()

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/transport_type.f90
      ${PELEPHYSICS_SOURCE_DIR}/transport.F90
   )

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Transport/${PELEC_TRANSPORT_MODEL}")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/actual_transport.f90
   )

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Reactions")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/network.f90
   )

   set(PELEPHYSICS_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Reactions/${PELEC_REACTIONS_MODEL}")
   add_sources(GlobalSourceList
      ${PELEPHYSICS_SOURCE_DIR}/actual_network.f90
      ${PELEPHYSICS_SOURCE_DIR}/actual_reactor.F90
   )
   if("${PELEC_REACTIONS_MODEL}" STREQUAL "Fuego")
      add_sources(GlobalSourceList
        ${PELEPHYSICS_SOURCE_DIR}/react_type.F90
      )
   endif()
endfunction(get_pelephysics_sources)
