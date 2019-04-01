# Find MASA library
#
# Set MASA_DIR to the base directory where the package is installed
#
# Sets three variables
#   - MASA_INCLUDE_DIRS
#   - MASA_MOD_DIRS
#   - MASA_LIBRARIES
#

find_path(MASA_INCLUDE_DIRS
  masa.h
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES include)

find_path(MASA_MOD_DIRS
  masa.mod
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

find_library(MASA_LIBRARIES
  NAMES masa fmasa
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MASA DEFAULT_MSG MASA_INCLUDE_DIRS MASA_MOD_DIRS MASA_LIBRARIES)
mark_as_advanced(MASA_INCLUDE_DIRS MASA_MOD_DIRS MASA_LIBRARIES)
