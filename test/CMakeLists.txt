# Make sure compiler optimization doesn't bork catch
set(CMAKE_BUILD_TYPE Debug)

include_directories("${SOURCE_DIR}" "${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

file(GLOB ECOEVOLITY_TEST_SOURCES "test_*.cpp")

set(ECOEVOLITY_SOURCES
    "${ECOEVOLITY_SOURCE_DIR}/error.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/matrix.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/likelihood.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/poptree.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/operator_schedule.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/operator.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/collection.cpp"
    "${EXTERNAL_PROJECT_DIR}/cpp-optparse/OptionParser.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/ecoevolity.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/simcoevolity.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/phycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/simphycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/sumphycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/data.cpp")

set(TEST_LIBRARIES_TO_LINK ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY})
if (BUILD_WITH_THREADS)
    set(TEST_LIBRARIES_TO_LINK ${TEST_LIBRARIES_TO_LINK} ${CMAKE_THREAD_LIBS_INIT})
endif()
if (USING_PLL)
    set(TEST_LIBRARIES_TO_LINK ${TEST_LIBRARIES_TO_LINK} ${PLL_LIBRARY})
endif()

#set(ECOEVOLITY_TEST_SOURCES "${ECOEVOLITY_TEST_SOURCES} ${TEST_DIR}/ecoevolity_testing.cpp")

# Using EXCLUDE_FROM_ALL and ADD_CUSTOM_TARGET so that tests compile and run
# ONLY when 'make check' target is called
add_executable(test_ecoevolity EXCLUDE_FROM_ALL ecoevolity_testing.cpp ${ECOEVOLITY_SOURCES} ${ECOEVOLITY_TEST_SOURCES})
target_link_libraries(test_ecoevolity ${TEST_LIBRARIES_TO_LINK})

# add_test(NAME run_tests COMMAND test_ecoevolity)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS test_ecoevolity)
