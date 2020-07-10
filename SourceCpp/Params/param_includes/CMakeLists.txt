target_sources(${pelec_exe_name}
   PRIVATE
     pelec_defaults.H
     pelec_params.H
     pelec_queries.H 
   )

target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
