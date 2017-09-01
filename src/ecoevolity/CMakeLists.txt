include_directories("${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

set(ECOEVOLITY_MAIN_SOURCES
    ecoevolity_main.cpp
    simcoevolity_main.cpp
    dpprobs_main.cpp
    sumcoevolity_main.cpp)

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
    sumcoevolity.cpp
    data.cpp)

add_executable(ecoevolity
        ecoevolity_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(ecoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS ecoevolity
        DESTINATION bin
        )

add_executable(simcoevolity
        simcoevolity_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(simcoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS simcoevolity
        DESTINATION bin
        )

add_executable(dpprobs
        dpprobs_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(dpprobs ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS dpprobs
        DESTINATION bin
        )

add_executable(sumcoevolity
        sumcoevolity_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(sumcoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
install(TARGETS sumcoevolity
        DESTINATION bin
        )
