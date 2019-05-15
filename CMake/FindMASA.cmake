# Find MASA library
#
# Set MASA_DIR to the base directory where the package is installed
#
# Sets 4 variables
#   - MASA_INCLUDE_DIRS
#   - MASA_MOD_DIRS
#   - MASA_LIBRARY
#   - MASA_FORTRAN_LIBRARY
#

find_path(MASA_INCLUDE_DIRS
  masa.h
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES include)

find_path(MASA_MOD_DIRS
  masa.mod
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

find_library(MASA_LIBRARY
  NAMES masa
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

find_library(MASA_FORTRAN_LIBRARY
  NAMES fmasa
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MASA DEFAULT_MSG MASA_INCLUDE_DIRS MASA_MOD_DIRS MASA_LIBRARY MASA_FORTRAN_LIBRARY)
mark_as_advanced(MASA_INCLUDE_DIRS MASA_MOD_DIRS MASA_LIBRARY MASA_FORTRAN_LIBRARY)
