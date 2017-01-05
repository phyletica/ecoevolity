include_directories("${YAML_CPP_INCLUDE_DIR}" "${EXTERNAL_PROJECT_DIR}")

set(ECOEVOLITY_MAIN_SOURCES
    ecoevolity_main.cpp
    simcoevolity_main.cpp)

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
    data.cpp)

add_executable(ecoevolity
        ecoevolity_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(ecoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} Threads::Threads)
install(TARGETS ecoevolity
        DESTINATION bin
        )

add_executable(simcoevolity
        simcoevolity_main.cpp
        ${ECOEVOLITY_SOURCES}
        )
target_link_libraries(simcoevolity ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY} Threads::Threads)
install(TARGETS simcoevolity
        DESTINATION bin
        )
