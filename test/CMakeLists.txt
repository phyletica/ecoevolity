# Make sure compiler optimization doesn't bork catch
set(CMAKE_BUILD_TYPE Debug)

include_directories("${SOURCE_DIR}" "${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

if (EXHAUSTIVE_TESTING)
    MESSAGE(STATUS "EXHAUSTIVE_TESTING: ${EXHAUSTIVE_TESTING}")
    MESSAGE(STATUS "  Compiling all test_*.cpp unit tests")
    file(GLOB ECOEVOLITY_TEST_SOURCES "test_*.cpp")
else()
    MESSAGE(STATUS "EXHAUSTIVE_TESTING: ${EXHAUSTIVE_TESTING}")
    MESSAGE(STATUS "  Compiling subset of fast unit tests")
    set(ECOEVOLITY_TEST_SOURCES
        "${TEST_DIR}/test_error.cpp"
        "${TEST_DIR}/test_data.cpp"
        "${TEST_DIR}/test_general_tree_settings.cpp"
        "${TEST_DIR}/test_math_util.cpp"
        "${TEST_DIR}/test_matrix.cpp"
        # "${TEST_DIR}/test_parameter.cpp"
        "${TEST_DIR}/test_path.cpp"
        # "${TEST_DIR}/test_probability.cpp"
        # "${TEST_DIR}/test_rng.cpp"
        "${TEST_DIR}/test_settings.cpp"
        "${TEST_DIR}/test_split.cpp"
        "${TEST_DIR}/test_spreadsheet.cpp"
        # "${TEST_DIR}/test_stats_util.cpp"
        "${TEST_DIR}/test_string_util.cpp"
        "${TEST_DIR}/test_treecomp.cpp"
        "${TEST_DIR}/test_tree_io.cpp"
        "${TEST_DIR}/test_vector2d.cpp"
        # "${TEST_DIR}/test_wattersons_theta.cpp"
        "${TEST_DIR}/test_phycoeval_cli.cpp"
        "${TEST_DIR}/test_sumphycoeval_cli.cpp"
        "${TEST_DIR}/test_simphycoeval_cli.cpp"
        "${TEST_DIR}/test_ecoevolity_cli.cpp"
        "${TEST_DIR}/test_simcoevolity_cli.cpp"
        "${TEST_DIR}/test_wtheta_cli.cpp"
        )
endif()

set(ECOEVOLITY_SOURCES
    "${ECOEVOLITY_SOURCE_DIR}/error.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/matrix.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/likelihood.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/tree.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/operator_schedule.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/operator.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/collection.cpp"
    "${EXTERNAL_PROJECT_DIR}/cpp-optparse/OptionParser.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/ecoevolity.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/simcoevolity.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/phycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/simphycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/sumphycoeval.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/wtheta.cpp"
    "${ECOEVOLITY_SOURCE_DIR}/data.cpp")

if (BUILD_WITH_THREADS)
    set(TEST_LIBRARIES_TO_LINK ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
else()
    set(TEST_LIBRARIES_TO_LINK ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY})
endif()

#set(ECOEVOLITY_TEST_SOURCES "${ECOEVOLITY_TEST_SOURCES} ${TEST_DIR}/ecoevolity_testing.cpp")

# Using EXCLUDE_FROM_ALL and ADD_CUSTOM_TARGET so that tests compile and run
# ONLY when 'make check' target is called
add_executable(test_ecoevolity EXCLUDE_FROM_ALL ecoevolity_testing.cpp ${ECOEVOLITY_SOURCES} ${ECOEVOLITY_TEST_SOURCES})
target_link_libraries(test_ecoevolity ${TEST_LIBRARIES_TO_LINK})

# add_test(NAME run_tests COMMAND test_ecoevolity)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS test_ecoevolity)
