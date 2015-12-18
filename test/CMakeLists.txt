# Make sure compiler optimization doesn't bork catch
set(CMAKE_BUILD_TYPE Debug)

include_directories("${SOURCE_DIR}")

file(GLOB ECOEVOLITY_TEST_SOURCES "test_*.cpp")

#set(ECOEVOLITY_TEST_SOURCES "${ECOEVOLITY_TEST_SOURCES} ${TEST_DIR}/ecoevolity_testing.cpp")

# Using EXCLUDE_FROM_ALL and ADD_CUSTOM_TARGET so that tests compile and run
# ONLY when 'make check' target is called
add_executable(test_ecoevolity EXCLUDE_FROM_ALL ecoevolity_testing.cpp ${ECOEVOLITY_TEST_SOURCES})

add_test(NAME run_tests COMMAND test_ecoevolity)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS test_ecoevolity)
