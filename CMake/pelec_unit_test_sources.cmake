#Aggregate PeleC unit test source files
function(get_pelec_unit_test_sources)
  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Testing/unit_tests")
  add_sources(GlobalUnitSourceList
     ${PELEC_SOURCE_DIR}/unit-tests-${PELEC_DIM}D.C
     ${PELEC_SOURCE_DIR}/unit-test-${PELEC_DIM}D-1.C
  )
endfunction(get_pelec_unit_test_sources)
