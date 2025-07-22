# Find PLL libraries
# Once done this will define
#
# PLL_FOUND - System has PLL library
# PLL_STATIC_FOUND - System has statically linked PLL library (libpll.a)
# PLL_INCLUDE_DIR - The path to PLL headers 
# PLL_LIBRARY - The path to PLL libraries
# PLL_STATIC_LIBRARY - The path to libpll.a 
# PLL_CFLAGS - Compiler -I switches for pll 
# PLL_LIBS - Compiler -L switches for ppll

# Get PLL location info from user
#   Commandline flag gets top priority, then environmental variable, and then
#   defaults to empty string
set (PLL_PREFIX "$ENV{PLL_PREFIX}" CACHE PATH "Path to PLL installation")
set (PLL_INCLUDE_PREFIX "")
set (PLL_LIB_PREFIX "")
set (PLL_INCLUDE_DIR "")
set (PLL_LIBRARY "")
set (PLL_STATIC_LIBRARY "")
if (NOT ${PLL_PREFIX} STREQUAL "")
    message(STATUS "PLL_PREFIX was set from command line or shell environment")
    set (PLL_INCLUDE_PREFIX "${PLL_PREFIX}/include")
    set (PLL_LIB_PREFIX "${PLL_PREFIX}/lib")
    if (EXISTS "${PLL_INCLUDE_PREFIX}/libpll")
        set (PLL_INCLUDE_DIR "${PLL_INCLUDE_PREFIX}/libpll")
    endif()
    if (EXISTS "${PLL_LIB_PREFIX}/libpll.a")
        set (PLL_LIBRARY "${PLL_LIB_PREFIX}/libpll.a")
        set (PLL_STATIC_LIBRARY "${PLL_LIB_PREFIX}/libpll.a")
    endif()
    if (EXISTS "${PLL_LIB_PREFIX}/libpll.so")
        set (PLL_LIBRARY "${PLL_LIB_PREFIX}/libpll.so")
    endif()
else()
    message(STATUS "Could not find PLL_PREFIX from command line or shell environment")
endif()
message(STATUS "PLL_PREFIX: ${PLL_PREFIX}")

# find_package(PkgConfig)

# if (PKG_CONFIG_FOUND)
#     pkg_check_modules(
#     message(STATUS "Pkg config found")
#     message(STATUS "PLL_INCUDEDIR: ${PC_PLL_INCUDEDIR}")
#     message(STATUS "PLL_INCUDEDIRS: ${PC_PLL_INCUDEDIRS}")
#     message(STATUS "PLL_LIBDIR: ${PC_PLL_LIBDIR}")
#     message(STATUS "PLL_LIBRARY_DIRS: ${PC_PLL_LIBRARY_DIRS}")
#     message(STATUS "PLL_CFLAGS: ${PC_PLL_CFLAGS}")
#     message(STATUS "PLL_CFLAGS_OTHER: ${PC_PLL_CFLAGS_OTHER}")
#     message(STATUS "PLL_LIBS: ${PC_PLL_LIBS}")
#     message(STATUS "PLL_LIBS_OTHER: ${PC_PLL_LIBS_OTHER}")
#     set (PLL_CFLAGS "${PC_PLL_CFLAGS_OTHER}")
#     set (PLL_LIBS "${PC_PLL_LIBS_OTHER}")
# endif ()

if ("${PLL_INCLUDE_DIR}" STREQUAL "")
    find_path (PLL_INCLUDE_DIR
            NAMES "pll.h"
            PATH_SUFFIXES "include/libpll" "libpll"
            HINTS
            "${PLL_PREFIX}"
            "${PLL_INCLUDE_PREFIX}"
            # "${PC_PLL_INCUDEDIR}"
            # "${PC_PLL_INCLUDE_DIRS}"
            "${CMAKE_INCLUDE_PATH}"
            /usr/local/include
            /usr/include
            ~/Library/Frameworks/include/libpll
            /Library/Frameworks/include/libpll
            /sw/include/libpll/         # Fink
            /opt/local/include/libpll/  # DarwinPorts
            /opt/csw/include/libpll     # Blastwave
            /opt/include/libpll
            )
endif()

if ("${PLL_LIBRARY}" STREQUAL "")
    find_library (PLL_LIBRARY
            NAMES "libpll"
            PATH_SUFFIXES "lib"
            HINTS
            "${PLL_LIB_PREFIX}"
            # "${PC_PLL_LIBDIR}"
            # "${PC_PLL_LIBRARY_DIRS}"
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
endif()

if ("${PLL_STATIC_LIBRARY}" STREQUAL "")
    find_library (PLL_STATIC_LIBRARY
            NAMES "libpll.a"
            PATH_SUFFIXES "lib"
            HINTS
            "${PLL_LIB_PREFIX}"
            # "${PC_PLL_LIBDIR}"
            # "${PC_PLL_LIBRARY_DIRS}"
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
endif()

if (EXISTS "${PLL_INCLUDE_DIR}")
    set (PLL_FOUND TRUE)
else()
    set (PLL_FOUND FALSE)
endif()

if (EXISTS "${PLL_STATIC_LIBRARY}")
    set (PLL_STATIC_FOUND TRUE)
else()
    set (PLL_STATIC_FOUND FALSE)
endif()

message(STATUS "PLL_FOUND: ${PLL_FOUND}")
message(STATUS "PLL_STATIC_FOUND: ${PLL_STATIC_FOUND}")
message(STATUS "PLL_INCLUDE_DIR: ${PLL_INCLUDE_DIR}")
message(STATUS "PLL_LIBRARY: ${PLL_LIBRARY}")
message(STATUS "PLL_STATIC_LIBRARY: ${PLL_STATIC_LIBRARY}")
# message(STATUS "PLL_CFLAGS: ${PLL_CFLAGS}")
# message(STATUS "PLL_LIBS: ${PLL_LIBS}")
