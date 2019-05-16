#Aggregate PeleC unit test source files
function(get_pelec_unit_test_sources)
  set(PELEC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Testing/unit_tests")
  add_sources(GlobalUnitSourceList
     ${PELEC_SOURCE_DIR}/unit-tests-${PELEC_DIM}d.C
     ${PELEC_SOURCE_DIR}/unit-test-${PELEC_DIM}d-1.C
  )
endfunction(get_pelec_unit_test_sources)
