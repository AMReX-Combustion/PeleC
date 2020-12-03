find_package(Sphinx)

if(SPHINX_FOUND)
  add_custom_target(docs)

  find_package(Doxygen)
  find_package(LATEX COMPONENTS PDFLATEX BIBTEX)

  add_custom_target(sphinx
    COMMAND ${SPHINX_EXECUTABLE} -M html
    "${CMAKE_CURRENT_SOURCE_DIR}/sphinx" "${CMAKE_CURRENT_BINARY_DIR}/sphinx"
    -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
  add_dependencies(docs sphinx)

  if(LATEX_FOUND)
    add_custom_target(sphinx-pdf
      COMMAND ${SPHINX_EXECUTABLE} -M latexpdf
      "${CMAKE_CURRENT_SOURCE_DIR}/sphinx" "${CMAKE_CURRENT_BINARY_DIR}/sphinx"
      -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx ${SPHINX_OPTS})
    add_dependencies(docs sphinx-pdf)
  endif()

  if(DOXYGEN_FOUND)
    configure_file(Doxyfile.in Doxyfile @ONLY)
    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" )
    add_dependencies(docs doxygen)
    add_dependencies(doxygen sphinx)
  endif()
else()
    message(WARNING "Skipping Sphinx documentation due to missing dependencies.")
endif()
