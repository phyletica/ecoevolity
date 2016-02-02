set(ECOEVOLITY_MAIN_SOURCES
    ecoevolity_main.cpp)

set(ECOEVOLITY_SOURCES
    util.cpp
    data.cpp)

foreach(main_source ${ECOEVOLITY_MAIN_SOURCES})
    get_filename_component(main_file_name ${main_source} NAME_WE)
    string(REGEX REPLACE "_main" "" stripped_name ${main_file_name})
    string(REGEX REPLACE "_" "-" exe_name ${stripped_name})
    add_executable(${exe_name}
            ${main_source}
            ${ECOEVOLITY_SOURCES}
            )
    TARGET_LINK_LIBRARIES(${exe_name} ${NCL_LIBRARIES})
    INSTALL(TARGETS ${exe_name}
            DESTINATION bin
            )
endforeach()