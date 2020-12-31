# Find MASA library
#
# Set MASA_DIR to the base directory where the package is installed
#
# Sets 3 variables
#   - MASA_INCLUDE_DIRS
#   - MASA_LIBRARY
#   - MASA_VERSION
# Sets target
#    - MASA::MASA

find_package(PkgConfig QUIET)
pkg_check_modules(PC_MASA QUIET MASA)

find_path(MASA_INCLUDE_DIRS
    NAMES masa.h
    PATHS ${PC_MASA_INCLUDE_DIRS}
    HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES include
)

find_library(MASA_LIBRARY
  NAMES masa
  HINTS ${MASA_DIR} ${CMAKE_INSTALL_PREFIX}
  PATH_SUFFIXES lib)

set(MASA_VERSION ${PC_MASA_VERSION})

mark_as_advanced(MASA_FOUND MASA_INCLUDE_DIRS MASA_LIBRARY MASA_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MASA
    REQUIRED_VARS MASA_INCLUDE_DIRS MASA_LIBRARY
    VERSION_VAR MASA_VERSION
)

if(MASA_FOUND AND NOT TARGET MASA::MASA)
    add_library(MASA::MASA INTERFACE IMPORTED)
    set_target_properties(MASA::MASA PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${MASA_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${MASA_LIBRARY}"
    )
endif()
