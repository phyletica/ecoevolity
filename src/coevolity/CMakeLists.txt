set(COEVOLITY_SOURCES
    util.cpp)

file(GLOB COEVOLITY_MAIN_SOURCES "*_main.cpp")

foreach(main_source ${COEVOLITY_MAIN_SOURCES})
    get_filename_component(main_file_name ${main_source} NAME_WE)
    string(REGEX REPLACE "_main" "" stripped_name ${main_file_name})
    string(REGEX REPLACE "_" "-" exe_name ${stripped_name})
    add_executable(${exe_name}
            ${main_source}
            ${COEVOLITY_SOURCES}
            )
    # TARGET_LINK_LIBRARIES(${exe_name} ${})
    INSTALL(TARGETS ${exe_name}
            DESTINATION bin
            )
endforeach()
