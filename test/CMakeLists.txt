# Make sure compiler optimization doesn't bork catch
set(CMAKE_BUILD_TYPE Debug)

include_directories("${SOURCE_DIR}" "${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

file(GLOB ECOEVOLITY_TEST_SOURCES "test_*.cpp")

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
    "${ECOEVOLITY_SOURCE_DIR}/data.cpp")

#set(ECOEVOLITY_TEST_SOURCES "${ECOEVOLITY_TEST_SOURCES} ${TEST_DIR}/ecoevolity_testing.cpp")

# Using EXCLUDE_FROM_ALL and ADD_CUSTOM_TARGET so that tests compile and run
# ONLY when 'make check' target is called
add_executable(test_ecoevolity EXCLUDE_FROM_ALL ecoevolity_testing.cpp ${ECOEVOLITY_SOURCES} ${ECOEVOLITY_TEST_SOURCES})
target_link_libraries(test_ecoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})

# add_test(NAME run_tests COMMAND test_ecoevolity)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS test_ecoevolity)
