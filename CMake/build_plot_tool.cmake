function(build_plot_tool PLOT_TOOL_NAME)

  #Expose functions we want to be able to call
  include(${CMAKE_SOURCE_DIR}/CMake/add_source_function.cmake)
  include(${CMAKE_SOURCE_DIR}/CMake/plot_tool_sources.cmake)

  #Clear source file list from any previous executables
  set_property(GLOBAL PROPERTY GlobalPlotToolSourceList "") 

  #Aggregate amrex and pelephysics source files
  get_plot_tool_sources(${PLOT_TOOL_NAME})
  
  #Put source list from global property into local list
  get_property(PLOT_TOOL_SOURCES GLOBAL PROPERTY GlobalPlotToolSourceList)

  #Create an executable based on all the source files we aggregated
  add_executable(${PLOT_TOOL_NAME} ${PLOT_TOOL_SOURCES})

  #Set definitions for our particular executable
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE BL_SPACEDIM=3)
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE BL_FORT_USE_UNDERSCORE)
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE AMREX_SPACEDIM=3)
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE AMREX_FORT_USE_UNDERSCORE)
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:BL_LANG_FORT>)
  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:AMREX_LANG_FORT>)
  # CMake BUILD_TYPE should already define this for us
  #if(${CMAKE_BUILD_TYPE} MATCHES "Release")
  #  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE NDEBUG)
  #endif()
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE BL_Darwin)
    target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE AMREX_Darwin)
  endif()
  
  #AMReX include directories
  target_include_directories(${PLOT_TOOL_NAME} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Src/Base)
  target_include_directories(${PLOT_TOOL_NAME} SYSTEM PRIVATE ${CMAKE_SOURCE_DIR}/Submodules/AMReX/Tools/Plotfile)
  
  target_link_libraries(${PLOT_TOOL_NAME} Threads::Threads)

  #Link our executable to the MPI libraries, etc
  #if(PELEC_ENABLE_MPI)
  #  target_link_libraries(${PLOT_TOOL_NAME} PRIVATE MPI::MPI_CXX MPI::MPI_C MPI::MPI_Fortran)
  #  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE BL_USE_MPI)
  #  target_compile_definitions(${PLOT_TOOL_NAME} PRIVATE AMREX_USE_MPI)
  #endif()
 
  #Keep our Fortran module files confined to a unique directory for each executable 
  set_target_properties(${PLOT_TOOL_NAME} PROPERTIES Fortran_MODULE_DIRECTORY
                       "${CMAKE_BINARY_DIR}/fortran_modules/${PLOT_TOOL_NAME}_fortran_modules")

  #Define what we want to be installed during a make install 
  install(TARGETS ${PLOT_TOOL_NAME}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)
  
endfunction(build_plot_tool PLOT_TOOL_NAME)
