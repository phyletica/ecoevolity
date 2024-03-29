cmake_minimum_required (VERSION 3.1)

project("ecoevolity" C CXX)


#######################################################################
# Get local cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
include(PreventInSourceBuild)
include(FinalFile)


#######################################################################
# Course versioning
set(PROJECT_MAJOR_VERSION 1)
set(PROJECT_MINOR_VERSION 0)
set(PROJECT_PATCH_VERSION 0)
set(PROJECT_VERSION
    ${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION})


#######################################################################
# Testing
enable_testing(true)  # Enables unit-testing.


#######################################################################
# Project paths
get_filename_component(BASE_DIR "${PROJECT_SOURCE_DIR}" REALPATH)
set(SOURCE_DIR "${BASE_DIR}/src")
set(ECOEVOLITY_SOURCE_DIR "${SOURCE_DIR}/ecoevolity")
set(EXTERNAL_PROJECT_DIR "${SOURCE_DIR}/external")
set(NCL_SOURCE_DIR "${EXTERNAL_PROJECT_DIR}/ncl")
set(YAML_CPP_SOURCE_DIR "${EXTERNAL_PROJECT_DIR}/yaml-cpp-master-ce056ac")
set(TEST_DIR "${BASE_DIR}/test")
set(TEST_DATA_DIR "${TEST_DIR}/data")


#######################################################################
# Find and set up NCL variables
option (FORCE_BUNDLED_NCL "Use included NCL library" ON)
include(FindNCL)
include_directories("${NCL_INCLUDE_DIRS}")


#######################################################################
# Find and set up yaml-cpp variables
option (FORCE_BUNDLED_YAML_CPP "Use included yaml-cpp library" ON)
if (NOT FORCE_BUNDLED_YAML_CPP)
    include(FindYamlCpp)
else()
    set(YAML_CPP_FOUND OFF)
endif()

set(USING_BUNDLED_YAML_CPP OFF)
if (NOT YAML_CPP_FOUND)
    set(YAML_CPP_LIBRARY yaml)
    set(YAML_CPP_INCLUDE_DIR "${YAML_CPP_SOURCE_DIR}/include")
    set(USING_BUNDLED_YAML_CPP ON)
endif()
MESSAGE(STATUS "USING_BUNDLED_YAML_CPP: ${USING_BUNDLED_YAML_CPP}")


#######################################################################
# Find and set up threads variables
option (BUILD_WITH_THREADS "Build to use multithreading" OFF)
if (BUILD_WITH_THREADS)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    add_definitions(-DBUILD_WITH_THREADS)
endif()
MESSAGE(STATUS "BUILD_WITH_THREADS: ${BUILD_WITH_THREADS}")


#######################################################################
# Find and set up variables for parameterization of root size
option (BUILD_WITH_ABSOLUTE_ROOT_SIZE "Build so prior on root size is absolute" OFF)
if (BUILD_WITH_ABSOLUTE_ROOT_SIZE)
    add_definitions(-DBUILD_WITH_ABSOLUTE_ROOT_SIZE)
endif()
MESSAGE(STATUS "BUILD_WITH_ABSOLUTE_ROOT_SIZE: ${BUILD_WITH_ABSOLUTE_ROOT_SIZE}")


#######################################################################
# Compiler settings 
option (STATIC_LINKING "Build with static links" OFF)
if (STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()
MESSAGE(STATUS "Using static linking: ${STATIC_LINKING}")


if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
MESSAGE(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-null-sentinel -Wundef -Werror -Wno-unused -Wno-deprecated-copy -Wstrict-overflow=2")
    # NCL errors: -Wold-style-cast -Wsign-conversion 
    # yaml-cpp errors: -Wshadow -Wswitch-default
    if (${USING_BUNDLED_NCL} OR ${USING_BUNDLED_YAML_CPP})
        set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic -Wall -Wno-unused")
    endif()
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic -Wall -Wno-unused")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELEASE} -g")
endif()


#######################################################################
# Fine versioning

# configure a version file to pass CMake and Git settings to source code
include(TrackGitRevision)
get_git_revision_info(
    PROJECT_GIT_SHA1
    PROJECT_GIT_SHA1_SHORT
    PROJECT_GIT_BRANCH
    PROJECT_GIT_COMMIT_DATE
    PROJECT_GIT_COMMIT_TIMESTAMP)
if (PROJECT_GIT_SHA1)
    set(PROJECT_SOURCE_IDENTIFIER
        "${PROJECT_GIT_BRANCH} ${PROJECT_GIT_SHA1_SHORT}: ${PROJECT_GIT_COMMIT_DATE}")
    set(PROJECT_DETAILED_VERSION
        "${PROJECT_VERSION} (${PROJECT_SOURCE_IDENTIFIER})")
    set(PROJECT_SOURCE_IDENTIFIED 1)
else()
    set(PROJECT_SOURCE_IDENTIFIED 0)
    set(PROJECT_DETAILED_VERSION
        "${PROJECT_VERSION}")
endif()

# #cmakedefine creates a valid macro with the contents of the variable if
# it is not empty, and undefines and comments out the line if not
file(WRITE ${ECOEVOLITY_SOURCE_DIR}/version.hpp.in
"/**
 * This file is automatically generated during cmake configuration.
 * Do NOT modify the contents here; modify the root 'CMakelists.txt'
 * file instead.
 */

#ifndef ECOEVOLITY_VERSION_HPP
#define ECOEVOLITY_VERSION_HPP
#cmakedefine01 PROJECT_SOURCE_IDENTIFIED
#cmakedefine PROJECT_GIT_SHA1               \"${PROJECT_GIT_SHA1}\"
#cmakedefine PROJECT_GIT_SHA1_SHORT         \"${PROJECT_GIT_SHA1_SHORT}\"
#cmakedefine PROJECT_GIT_BRANCH             \"${PROJECT_GIT_BRANCH}\"
#cmakedefine PROJECT_GIT_COMMIT_DATE        \"${PROJECT_GIT_COMMIT_DATE}\"
#cmakedefine PROJECT_GIT_COMMIT_TIMESTAMP   \"${PROJECT_GIT_COMMIT_TIMESTAMP}\"
#cmakedefine PROJECT_SOURCE_IDENTIFIER      \"${PROJECT_SOURCE_IDENTIFIER}\"
#cmakedefine PROJECT_VERSION                \"${PROJECT_VERSION}\"
#cmakedefine PROJECT_DETAILED_VERSION       \"${PROJECT_DETAILED_VERSION}\"
#define PROJECT_NAME                   \"Ecoevolity\"
#endif
")
configure_file(
    "${ECOEVOLITY_SOURCE_DIR}/version.hpp.in"
    "${ECOEVOLITY_SOURCE_DIR}/version.hpp" @ONLY)

# add binary tree to search path so version.hpp can be found
include_directories("${PROJECT_BINARY_DIR}")


#######################################################################
# Add subdirectories 
add_subdirectory("${SOURCE_DIR}")
add_subdirectory("${ECOEVOLITY_SOURCE_DIR}")
add_subdirectory("${EXTERNAL_PROJECT_DIR}")
add_subdirectory("${TEST_DIR}")
add_subdirectory("${TEST_DATA_DIR}")


#######################################################################
# CPack settings 
set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set (CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_MAJOR_VERSION})
set (CPACK_PACKAGE_VERSION_MINOR ${PROJECT_MINOR_VERSION})
set (CPACK_PACKAGE_VERSION_PATCH ${PROJECT_PATCH_VERSION})
set (CPACK_GENERATOR "TGZ")
set (CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${PROJECT_MAJOR_VERSION}.${PROJECT_MINOR_VERSION}.${PROJECT_PATCH_VERSION}-${PROJECT_GIT_BRANCH}-${PROJECT_GIT_SHA1_SHORT}-${PROJECT_GIT_COMMIT_TIMESTAMP}"
    )
set (CPACK_SOURCE_IGNORE_FILES
    "${CMAKE_CURRENT_BINARY_DIR}/;/.git/;~$;.swp$;${CPACK_SOURCE_IGNORE_FILES}"
    )
include (CPack)

add_custom_target (dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
