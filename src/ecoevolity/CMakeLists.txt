include_directories("${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

set(ECOEVOLITY_SOURCES
    error.cpp
    matrix.cpp
    likelihood.cpp
    tree.cpp
    operator_schedule.cpp
    operator.cpp
    collection.cpp
    "${EXTERNAL_PROJECT_DIR}/cpp-optparse/OptionParser.cpp"
    ecoevolity.cpp
    simcoevolity.cpp
    dpprobs.cpp
    wdprobs.cpp
    sumcoevolity.cpp
    data.cpp)

set(ECOEVOLITY_LIBRARY libecoevolity)
set(ECOEVOLITY_INCLUDE_DIR "${ECOEVOLITY_SOURCE_DIR}")
CREATE_FINAL_FILE ("${CMAKE_CURRENT_BINARY_DIR}/libecoevolity.cpp" "${ECOEVOLITY_SOURCES}")
add_library(libecoevolity STATIC "${CMAKE_CURRENT_BINARY_DIR}/libecoevolity.cpp")
target_include_directories(libecoevolity PUBLIC "${ECOEVOLITY_SOURCE_DIR}")


if (BUILD_WITH_THREADS)
    set(LIBRARIES_TO_LINK ${ECOEVOLITY_LIBRARY} ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
else()
    set(LIBRARIES_TO_LINK ${ECOEVOLITY_LIBRARY} ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY})
endif()

add_executable(ecoevolity ecoevolity_main.cpp)
target_link_libraries(ecoevolity ${LIBRARIES_TO_LINK})
install(TARGETS ecoevolity
        DESTINATION bin)

add_executable(simcoevolity simcoevolity_main.cpp)
target_link_libraries(simcoevolity ${LIBRARIES_TO_LINK})
install(TARGETS simcoevolity
        DESTINATION bin)

add_executable(dpprobs dpprobs_main.cpp)
target_link_libraries(dpprobs ${LIBRARIES_TO_LINK})
install(TARGETS dpprobs
        DESTINATION bin)

add_executable(wdprobs wdprobs_main.cpp)
target_link_libraries(wdprobs ${LIBRARIES_TO_LINK})
install(TARGETS wdprobs
        DESTINATION bin)

add_executable(sumcoevolity sumcoevolity_main.cpp)
target_link_libraries(sumcoevolity ${LIBRARIES_TO_LINK})
install(TARGETS sumcoevolity
        DESTINATION bin)
