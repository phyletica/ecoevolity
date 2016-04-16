if (NOT YAML_CPP_FOUND)
    message (STATUS "Using bundled libyaml-cpp")

    aux_source_directory("${YAML_CPP_SOURCE_DIR}/src" yaml_src)
    CREATE_FINAL_FILE ("${CMAKE_CURRENT_BINARY_DIR}/libyamlcpp.cpp" "${yaml_src}")

    add_library(yaml STATIC "${CMAKE_CURRENT_BINARY_DIR}/libyamlcpp.cpp")
    target_include_directories(yaml PUBLIC "${YAML_CPP_SOURCE_DIR}/include")
    set_property(TARGET yaml PROPERTY FOLDER 3rdparty)

    if (CMAKE_COMPILER_IS_GNUCXX)
        target_compile_options(yaml PRIVATE -Wno-deprecated-declarations)
    endif()
endif()
