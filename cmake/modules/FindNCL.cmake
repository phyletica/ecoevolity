# Find NCL libraries
# Once done this will define
#
# NCL_FOUND - System has NCL library
# NCL_STATIC_FOUND - System has statically linked NCL library (libncl.a)
# NCL_CBLAS_FOUND - System has nclcblas library
# NCL_CBLAS_STATIC_FOUND - System has libnclcblas.a library
# USING_BUNDLED_NCL - Boolean of whether the bundled NCL is being used
# NCL_INCLUDE_DIRS - The path to NCL headers 
# NCL_LIBRARIES - The path to NCL libraries
# NCL_STATIC_LIBRARIES - The path to libncl.a 
# NCL_CBLAS_LIBRARIES - The path to nclcblas
# NCL_CBLAS_STATIC_LIBRARIES - The path to libnclcblas.a
# NCL_DEFINITIONS - Compiler switches required for using ncl

# Get NCL location info from user
#   Commandline flag gets top priority, then environmental variable, and then
#   defaults to empty string
set (NCL_PREFIX "$ENV{NCL_PREFIX}" CACHE PATH "Path to NCL installation")
set (NCL_INCLUDE_PREFIX "")
set (NCL_LIB_PREFIX "")
if (${NCL_PREFIX})
    set (NCL_INCLUDE_PREFIX "${NCL_PREFIX}/include")
    set (NCL_LIB_PREFIX "${NCL_PREFIX}/lib")
endif()
message(STATUS "NCL_PREFIX: ${NCL_PREFIX}")

include (FindPkgConfig)

if (PKG_CONFIG_FOUND)
    if (NCL_FIND_VERSION)
        set (_PACKAGE_ARGS "nclv2.1>=${GSL_FIND_VERSION}")
    else ()
        set (_PACKAGE_ARGS "nclv2.1")
    endif ()
    if (NCL_FIND_REQUIRED)
        set(_PACKAGE_ARGS "${_PACKAGE_ARGS}" REQUIRED)
    endif ()
    pkg_check_modules (PC_NCL "${_PACKAGE_ARGS}")
    set (NCL_DEFINITIONS "${PC_NCL_CFLAGS_OTHER}")
    set (NCL_VERSION_STRING "${PC_NCL_VERSION}")
endif ()

find_path (NCL_INCLUDE_DIR
        NAMES "ncl.h"
        HINTS
        "${NCL_PREFIX}"
        "${PC_NCL_INCUDEDIR}"
        "${PC_NCL_INCLUDE_DIRS}"
        "${CMAKE_INCLUDE_PATH}"
        PATH_SUFFIXES "ncl"
        )
find_library (NCL_LIBRARY
        NAMES "ncl"
        HINTS
        "${NCL_LIB_PREFIX}"
        "${PC_NCL_LIBDIR}"
        "${PC_NCL_LIBRARY_DIRS}"
        "${CMAKE_LIBRARY_PATH}"
        PATH_SUFFIXES "ncl"
        )
find_library (NCL_STATIC_LIBRARY
        NAMES "libncl.a"
        HINTS
        "${NCL_LIB_PREFIX}"
        "${PC_NCL_LIBDIR}"
        "${PC_NCL_LIBRARY_DIRS}"
        "${CMAKE_LIBRARY_PATH}"
        PATH_SUFFIXES "ncl"
        )
find_library (NCL_CBLAS_LIBRARY
        NAMES "nclcblas"
        HINTS
        "${NCL_LIB_PREFIX}"
        "${PC_NCL_LIBDIR}"
        "${PC_NCL_LIBRARY_DIRS}"
        "${CMAKE_LIBRARY_PATH}"
        PATH_SUFFIXES "ncl"
        )
find_library (NCL_CBLAS_STATIC_LIBRARY
        NAMES "libnclcblas.a"
        HINTS
        "${NCL_LIB_PREFIX}"
        "${PC_NCL_LIBDIR}"
        "${PC_NCL_LIBRARY_DIRS}"
        "${CMAKE_LIBRARY_PATH}"
        PATH_SUFFIXES "ncl"
        )

if (EXISTS "${NCL_INCLUDE_DIR}")
    set (NCL_FOUND TRUE)
else()
    set (NCL_FOUND FALSE)
endif()

if (EXISTS "${NCL_LIBRARY}")
    set (NCL_LIBRARIES "${NCL_LIBRARY}")
else()
    set (NCL_LIBRARIES "NOTFOUND")
endif()

if (EXISTS "${NCL_STATIC_LIBRARY}")
    set (NCL_STATIC_LIBRARIES "${NCL_STATIC_LIBRARY}")
    set (NCL_STATIC_FOUND TRUE)
else()
    set (NCL_STATIC_LIBRARIES "NOTFOUND")
    set (NCL_STATIC_FOUND FALSE)
endif()

if (EXISTS "${NCL_CBLAS_LIBRARY}")
    set (NCL_CBLAS_LIBRARIES "${NCL_CBLAS_LIBRARY}")
    set (NCL_CBLAS_FOUND TRUE)
else()
    set (NCL_CBLAS_LIBRARIES "NOTFOUND")
    set (NCL_CBLAS_FOUND FALSE)
endif()

if (EXISTS "${NCL_CBLAS_STATIC_LIBRARY}")
    set (NCL_CBLAS_STATIC_LIBRARIES "${NCL_CBLAS_STATIC_LIBRARY}")
    set (NCL_CBLAS_STATIC_FOUND TRUE)
else()
    set (NCL_CBLAS_STATIC_LIBRARIES "NOTFOUND")
    set (NCL_CBLAS_STATIC_FOUND FALSE)
endif()

set (USING_BUNDLED_NCL FALSE)
if (DEFINED NCL_SOURCE_DIR AND EXISTS "${NCL_SOURCE_DIR}" AND NOT ${NCL_FOUND})
    set(NCL_LIBRARY "ncl")
    add_library (${NCL_LIBRARY}
        "${NCL_SOURCE_DIR}/nxsassumptionsblock.cpp"
        "${NCL_SOURCE_DIR}/nxsblock.cpp"
        "${NCL_SOURCE_DIR}/nxscharactersblock.cpp"
        "${NCL_SOURCE_DIR}/nxscxxdiscretematrix.cpp"
        "${NCL_SOURCE_DIR}/nxsdatablock.cpp"
        "${NCL_SOURCE_DIR}/nxsdistancesblock.cpp"
        "${NCL_SOURCE_DIR}/nxsexception.cpp"
        "${NCL_SOURCE_DIR}/nxsmultiformat.cpp"
        "${NCL_SOURCE_DIR}/nxspublicblocks.cpp"
        "${NCL_SOURCE_DIR}/nxsreader.cpp"
        "${NCL_SOURCE_DIR}/nxssetreader.cpp"
        "${NCL_SOURCE_DIR}/nxsstring.cpp"
        "${NCL_SOURCE_DIR}/nxstaxaassociationblock.cpp"
        "${NCL_SOURCE_DIR}/nxstaxablock.cpp"
        "${NCL_SOURCE_DIR}/nxstoken.cpp"
        "${NCL_SOURCE_DIR}/nxstreesblock.cpp"
        "${NCL_SOURCE_DIR}/nxsunalignedblock.cpp"
        )
    set (NCL_INCLUDE_DIR "${NCL_SOURCE_DIR}")
    set (NCL_LIBRARIES "${NCL_LIBRARY}")
    set (USING_BUNDLED_NCL TRUE)
endif()

# Remove "ncl" from include path to support "#include "ncl/ncl.h" imports
if (EXISTS "${NCL_INCLUDE_DIR}")
    get_filename_component(NCL_INCLUDE_BASE "${NCL_INCLUDE_DIR}" DIRECTORY)
    set (NCL_INCLUDE_DIR "${NCL_INCLUDE_BASE}")
    set (NCL_INCLUDE_DIRS "${NCL_INCLUDE_DIR}")
else()
    set (NCL_INCLUDE_DIRS "NOTFOUND")
endif()
        
message(STATUS "NCL_FOUND: ${NCL_FOUND}")
message(STATUS "NCL_STATIC_FOUND: ${NCL_STATIC_FOUND}")
message(STATUS "NCL_CBLAS_FOUND: ${NCL_CBLAS_FOUND}")
message(STATUS "NCL_CBLAS_STATIC_FOUND: ${NCL_CBLAS_STATIC_FOUND}")
message(STATUS "USING_BUNDLED_NCL: ${USING_BUNDLED_NCL}")
message(STATUS "NCL_INCLUDE_DIRS: ${NCL_INCLUDE_DIRS}")
message(STATUS "NCL_LIBRARIES: ${NCL_LIBRARIES}")
message(STATUS "NCL_STATIC_LIBRARIES: ${NCL_STATIC_LIBRARIES}")
message(STATUS "NCL_CBLAS_LIBRARIES: ${NCL_CBLAS_LIBRARIES}")
message(STATUS "NCL_CBLAS_STATIC_LIBRARIES: ${NCL_CBLAS_STATIC_LIBRARIES}")
message(STATUS "NCL_DEFINITIONS: ${NCL_DEFINITIONS}")
