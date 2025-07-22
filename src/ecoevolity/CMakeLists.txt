include_directories("${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

set(ECOEVOLITY_SOURCES
    error.cpp
    matrix.cpp
    likelihood.cpp
    poptree.cpp
    operator_schedule.cpp
    operator.cpp
    collection.cpp
    "${EXTERNAL_PROJECT_DIR}/cpp-optparse/OptionParser.cpp"
    ecoevolity.cpp
    simcoevolity.cpp
    phycoeval.cpp
    simphycoeval.cpp
    sumphycoeval.cpp
    nex2yml.cpp
    dpprobs.cpp
    swprobs.cpp
    sumcoevolity.cpp
    data.cpp)

if (USING_PLL)
    set(ECOEVOLITY_SOURCES ${ECOEVOLITY_SOURCES} nephila.cpp)
endif()

set(ECOEVOLITY_LIBRARY libecoevolity)
set(ECOEVOLITY_INCLUDE_DIR "${ECOEVOLITY_SOURCE_DIR}")
CREATE_FINAL_FILE ("${CMAKE_CURRENT_BINARY_DIR}/libecoevolity.cpp" "${ECOEVOLITY_SOURCES}")
add_library(libecoevolity STATIC "${CMAKE_CURRENT_BINARY_DIR}/libecoevolity.cpp")
target_include_directories(libecoevolity PUBLIC "${ECOEVOLITY_SOURCE_DIR}")


set(LIBRARIES_TO_LINK ${ECOEVOLITY_LIBRARY} ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY})
if (BUILD_WITH_THREADS)
    set(LIBRARIES_TO_LINK ${LIBRARIES_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
endif()
if (USING_PLL)
    set(LIBRARIES_TO_LINK ${LIBRARIES_TO_LINK} ${PLL_LIBRARY})
endif()

set(ECOEVOLITY_EXE_NAME "ecoevolity")
set(SIMCOEVOLITY_EXE_NAME "simcoevolity")
if (BUILD_WITH_ABSOLUTE_ROOT_SIZE)
    set(ECOEVOLITY_EXE_NAME "${ECOEVOLITY_EXE_NAME}-ar")
    set(SIMCOEVOLITY_EXE_NAME "${SIMCOEVOLITY_EXE_NAME}-ar")
endif()

add_executable("${ECOEVOLITY_EXE_NAME}" ecoevolity_main.cpp)
target_link_libraries("${ECOEVOLITY_EXE_NAME}" ${LIBRARIES_TO_LINK})
install(TARGETS "${ECOEVOLITY_EXE_NAME}"
        DESTINATION bin)

add_executable("${SIMCOEVOLITY_EXE_NAME}" simcoevolity_main.cpp)
target_link_libraries("${SIMCOEVOLITY_EXE_NAME}" ${LIBRARIES_TO_LINK})
install(TARGETS "${SIMCOEVOLITY_EXE_NAME}"
        DESTINATION bin)

add_executable(dpprobs dpprobs_main.cpp)
target_link_libraries(dpprobs ${LIBRARIES_TO_LINK})
install(TARGETS dpprobs
        DESTINATION bin)

add_executable(swprobs swprobs_main.cpp)
target_link_libraries(swprobs ${LIBRARIES_TO_LINK})
install(TARGETS swprobs
        DESTINATION bin)

add_executable(sumcoevolity sumcoevolity_main.cpp)
target_link_libraries(sumcoevolity ${LIBRARIES_TO_LINK})
install(TARGETS sumcoevolity
        DESTINATION bin)

add_executable(nex2yml nex2yml_main.cpp)
target_link_libraries(nex2yml ${LIBRARIES_TO_LINK})
install(TARGETS nex2yml
        DESTINATION bin)

add_executable(phycoeval phycoeval_main.cpp)
target_link_libraries(phycoeval ${LIBRARIES_TO_LINK})
install(TARGETS phycoeval
        DESTINATION bin)

add_executable(simphycoeval simphycoeval_main.cpp)
target_link_libraries(simphycoeval ${LIBRARIES_TO_LINK})
install(TARGETS simphycoeval
        DESTINATION bin)

add_executable(sumphycoeval sumphycoeval_main.cpp)
target_link_libraries(sumphycoeval ${LIBRARIES_TO_LINK})
install(TARGETS sumphycoeval
        DESTINATION bin)

if (USING_PLL)
    add_executable(nephila nephila_main.cpp)
    target_link_libraries(nephila ${LIBRARIES_TO_LINK})
    install(TARGETS nephila
            DESTINATION bin)
endif()
