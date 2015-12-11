# Get VCS info and make available to program(s)

#
# This function forces a reconfigure on each Git commit
# so that the values of the variables propogated to the build
# system will be current.
#
# Derived from:
#       https://github.com/jeetsukumaran/pstrudel/blob/master/cmake/Modules/TrackGitRevision.cmake
#       2013 Jeet Sukumaran https://github.com/jeetsukumaran
#       Distributed under the GNU GPL Version 2
#       (See accompanying license at 
#       https://github.com/jeetsukumaran/pstrudel/blob/master/LICENSE.txt)
#
# Which was derived from:
#       GetGitRevisionDescription:get_git_head_revision()
#       2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
#       http://academic.cleardefinition.com
#       Iowa State University HCI Graduate Program/VRAC
#       Copyright Iowa State University 2009-2010.
#       Distributed under the Boost Software License, Version 1.0.
#       (See accompanying file LICENSE_1_0.txt or copy at
#       http://www.boost.org/LICENSE_1_0.txt)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component(_gitdescmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)
function(track_git_revision _gitrepofound)
	set(GIT_PARENT_DIR "${CMAKE_SOURCE_DIR}")
	set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	while(NOT EXISTS "${GIT_DIR}")	# .git dir not found, search parent directories
		set(GIT_PREVIOUS_PARENT "${GIT_PARENT_DIR}")
		get_filename_component(GIT_PARENT_DIR ${GIT_PARENT_DIR} PATH)
		if(GIT_PARENT_DIR STREQUAL GIT_PREVIOUS_PARENT)
			# We have reached the root directory, we are not in git
            set(${_gitrepofound} "GITREPO-NOTFOUND" PARENT_SCOPE)
			return()
		endif()
		set(GIT_DIR "${GIT_PARENT_DIR}/.git")
	endwhile()
	set(GIT_DATA "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/git-data")
	if(NOT EXISTS "${GIT_DATA}")
		file(MAKE_DIRECTORY "${GIT_DATA}")
	endif()
	if(NOT EXISTS "${GIT_DIR}/HEAD")
		return()
	endif()
	set(HEAD_FILE "${GIT_DATA}/HEAD")
	configure_file("${GIT_DIR}/HEAD" "${HEAD_FILE}" COPYONLY)
    configure_file("${_gitdescmoddir}/TrackGitRevision.cmake.in"
		"${GIT_DATA}/grabRef.cmake"
		@ONLY)
	include("${GIT_DATA}/grabRef.cmake")
    set(${_gitrepofound} "${GIT_DIR}" PARENT_SCOPE)
endfunction()

function(set_git_not_found _sha1 _shorthashvar _branchvar _commitdatevar _committimestampvar)
    set(${_sha1} "GITREPO-NOTFOUND" PARENT_SCOPE)
    set(${_shorthashvar} "GITREPO-NOTFOUND" PARENT_SCOPE)
    set(${_branchvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
    set(${_commitdatevar} "GITREPO_NOTFOUND" PARENT_SCOPE)
    set(${_committimestampvar} "GITREPO_NOTFOUND" PARENT_SCOPE)
endfunction()

# Retrieves revision information
function(get_git_revision_info _sha1 _shorthashvar _branchvar _commitdatevar _committimestampvar)
    track_git_revision(GIT_FOUND)
    IF (NOT GIT_FOUND)
        set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar _committimestampvar)
        return()
    ENDIF()
    FIND_PACKAGE(Git)
    IF(GIT_FOUND)
        EXECUTE_PROCESS(
            COMMAND ${GIT_EXECUTABLE} status 2>/dev/null
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            RESULT_VARIABLE GIT_STATUS_RV
            OUTPUT_QUIET
            ERROR_QUIET
            )
        IF (${GIT_STATUS_RV})
            MESSAGE(STATUS "Not a Git Repository")
            set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar _committimestampvar)
            return()
        ELSE()
            EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                # chain COMMAND's in a pipe by listing them one after the other; stdout
                # of one will be passed to stdin of subsequent
                # COMMAND sed etc. etc.
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_BRANCH_NAME
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_SHA1
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_SHA1_SHORT
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%cI
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_DATE
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%ct
                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_COMMIT_UNIX_TIMESTAMP
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
            SET(GIT_SOURCE_DESC "${GIT_BRANCH_NAME}:${GIT_SHA1_SHORT}, ${GIT_COMMIT_DATE}")
            MESSAGE(STATUS "${GIT_SOURCE_DESC}")
        ENDIF()
    ELSE()
        MESSAGE(STATUS "Git not found")
        set_git_not_found(_sha1 _shorthashvar _branchvar _commitdatevar _committimestampvar)
        return()
    ENDIF()
    set(${_sha1} "${GIT_SHA1}" PARENT_SCOPE)
    set(${_shorthashvar} "${GIT_SHA1_SHORT}" PARENT_SCOPE)
    set(${_branchvar} "${GIT_BRANCH_NAME}" PARENT_SCOPE)
    set(${_commitdatevar} "${GIT_COMMIT_DATE}" PARENT_SCOPE)
    set(${_committimestampvar} "${GIT_COMMIT_UNIX_TIMESTAMP}" PARENT_SCOPE)
endfunction()
