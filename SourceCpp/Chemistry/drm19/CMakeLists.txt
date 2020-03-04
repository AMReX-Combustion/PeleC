target_sources(${pelec_exe_name}
   PRIVATE
      chemistry_file.H
      mechanism.h
      mechanism.cpp
   )

target_include_directories(${pelec_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
