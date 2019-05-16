find_package(Doxygen)
find_package(Sphinx)

if (DOXYGEN_FOUND OR SPHINX_FOUND)
  add_custom_target(docs)
endif()

if(DOXYGEN_FOUND)
    configure_file(Doxyfile.in Doxyfile @ONLY)
    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile" )
    add_dependencies(docs doxygen)
else()
    message(WARNING "Skipping Doxygen documentation due to missing dependencies.")
endif()

if(SPHINX_FOUND)
  find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
  if (ENABLE_SPHINX_API_DOCS)
    set(SPHINX_OPTS "-t" "use_breathe")
    configure_file(Doxyfile.breathe.in Doxyfile.breathe @ONLY)
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
      DEPENDS ${HEADER}
      COMMAND ${DOXYGEN_EXECUTABLE} "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.breathe")
    set_source_files_properties(
      ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
      PROPERTIES GENERATED TRUE)
    add_custom_target(sphinx
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
      COMMAND ${SPHINX_EXECUTABLE} -M html
      "${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc" "${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc"
      -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc ${SPHINX_OPTS})
    add_dependencies(docs sphinx)

    # PDF documentation generation using sphinx -> latex -> pdflatex
    if (LATEX_FOUND)
      add_custom_target(sphinx-pdf
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen/xml/index.xml
        COMMAND ${SPHINX_EXECUTABLE} -M latexpdf
        "${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc" "${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc"
        -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc ${SPHINX_OPTS})
    endif()
  else()
    set(SPHINX_OPTS "")
    add_custom_target(sphinx
      COMMAND ${SPHINX_EXECUTABLE} -M html
      "${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc" "${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc"
      -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc ${SPHINX_OPTS})
    add_dependencies(docs sphinx)

    # PDF documentation generation using sphinx -> latex -> pdflatex
    if (LATEX_FOUND)
      add_custom_target(sphinx-pdf
        COMMAND ${SPHINX_EXECUTABLE} -M latexpdf
        "${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc" "${CMAKE_CURRENT_BINARY_DIR}/sphinx_doc"
        -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx_doc ${SPHINX_OPTS})
    endif()
  endif()
else()
    message(WARNING "Skipping Sphinx documentation due to missing dependencies.")
endif()
