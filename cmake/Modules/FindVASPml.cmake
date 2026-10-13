# Find VASPml library.
#
# This script tries to find the precompiled library libvaspml and the
# corresponding interface header inside the VASPml root directory which needs
# to be provided via:
#
# * a shell variable, e.g., "export VASPML_DIR=/path/to/vaspml", or,
# * a CMake cache variable, e.g. by passing "-DVASPML_DIR=/path/to/vaspml".

include(FindPackageHandleStandardArgs)

if (NOT DEFINED CACHE{VASPML_DIR})
  if (DEFINED ENV{VASPML_DIR})
    set(VASPML_DIR "$ENV{VASPML_DIR}" CACHE PATH "Path to VASPml root directory")
  else()
    message(FATAL_ERROR "Cannot search for VASPml library, help CMake by setting VASPML_DIR (-DVASPML_DIR=/path/to/vaspml or via environment variable).")
  endif()
endif()

find_library(VASPML_LIBVASPML NAMES vaspml HINTS "${VASPML_DIR}/lib")
find_path(VASPML_INCLUDE_DIR InterfaceLAMMPS.hpp HINTS "${VASPML_DIR}/build/libvaspml/")

find_package_handle_standard_args(VASPml DEFAULT_MSG
  VASPML_DIR
  VASPML_LIBVASPML
  VASPML_INCLUDE_DIR
)

if(VASPml_FOUND)
  if (NOT TARGET VASPml::VASPml)
    add_library(VASPml::LIBVASPML UNKNOWN IMPORTED)
    set_target_properties(VASPml::LIBVASPML PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES ${VASPML_INCLUDE_DIR}
      IMPORTED_LOCATION "${VASPML_LIBVASPML}")

    add_library(VASPml::VASPml INTERFACE IMPORTED)
    set_property(TARGET VASPml::VASPml PROPERTY
      INTERFACE_LINK_LIBRARIES VASPml::LIBVASPML)
  endif()
endif()

mark_as_advanced(
  VASPML_DIR
  VASPML_LIBVASPML
  VASPML_INCLUDE_DIR
)
