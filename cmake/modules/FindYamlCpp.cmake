# Find yaml-cpp libraries
# Once done this will define
#
# YAML_CPP_FOUND - System has yaml-cpp library
# YAML_CPP_INCLUDE_DIR - The path to yaml-cpp headers 
# YAML_CPP_LIBRARY - The path to yaml-cpp libraries

# Give priority to environmental variable
set (YAML_CPP_PREFIX "$ENV{YAML_CPP_PREFIX}" CACHE PATH "Path to yaml-cpp installation")
set (YAML_CPP_INCLUDE_PREFIX "")
set (YAML_CPP_LIB_PREFIX "")
if (${YAML_CPP_PREFIX})
    set (YAML_CPP_INCLUDE_PREFIX "${YAML_CPP_PREFIX}/include")
    set (YAML_CPP_LIB_PREFIX "${YAML_CPP_PREFIX}/lib")
endif()
message(STATUS "YAML_CPP_PREFIX: ${YAML_CPP_PREFIX}")

# if (YAML_CPP_STATIC_LIBRARY)
set(YAML_CPP_STATIC libyaml-cpp.a)
# endif()

find_path (YAML_CPP_INCLUDE_DIR
        NAMES "yaml-cpp"
        PATH_SUFFIXES "include"
        HINTS
        "${YAML_CPP_PREFIX}"
        "${YAML_CPP_INCLUDE_PREFIX}"
        "${CMAKE_INCLUDE_PATH}"
        /usr/local/include
        /usr/include
        ~/Library/Frameworks/yaml-cpp/include
        /Library/Frameworks/yaml-cpp/include
        /sw/yaml-cpp/           # Fink
        /opt/local/yaml-cpp/    # DarwinPorts
        /opt/csw/yaml-cpp/      # Blastwave
        /opt/yaml-cpp/
        )

find_library (YAML_CPP_LIBRARY
        NAMES "${YAML_CPP_STATIC}" "yaml-cpp"
        PATH_SUFFIXES "lib64" "lib"
        HINTS
        "${YAML_CPP_PREFIX}"
        "${YAML_CPP_LIB_PREFIX}"
        "${CMAKE_LIBRARY_PATH}"
        /usr/local
        /usr
        ~/Library/Frameworks
        /Library/Frameworks
        /sw
        /opt/local
        /opt/csw
        /opt
        )

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(YAML_CPP DEFAULT_MSG YAML_CPP_INCLUDE_DIR YAML_CPP_LIBRARY)
mark_as_advanced(YAML_CPP_INCLUDE_DIR YAML_CPP_LIBRARY)
message(STATUS "YAML_CPP_FOUND: ${YAML_CPP_FOUND}")
message(STATUS "YAML_CPP_INCLUDE_DIR: ${YAML_CPP_INCLUDE_DIR}")
message(STATUS "YAML_CPP_LIBRARY: ${YAML_CPP_LIBRARY}")
