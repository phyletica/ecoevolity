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

foreach(main_source ${ECOEVOLITY_MAIN_SOURCES})
    get_filename_component(main_file_name ${main_source} NAME_WE)
    string(REGEX REPLACE "_main" "" stripped_name ${main_file_name})
    string(REGEX REPLACE "_" "-" exe_name ${stripped_name})
    add_executable(${exe_name}
            ${main_source}
            ${ECOEVOLITY_SOURCES}
            )
    TARGET_LINK_LIBRARIES(${exe_name} ${NCL_LIBRARIES} ${YAML_CPP_LIBRARY})
    INSTALL(TARGETS ${exe_name}
            DESTINATION bin
            )
endforeach()
