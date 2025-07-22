file (GLOB test_data_files RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.yml"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.nex"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.phy"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.tre")

foreach (test_data_file ${test_data_files})
    configure_file(${test_data_file} ${CMAKE_CURRENT_BINARY_DIR}/${test_data_file} COPYONLY)
endforeach (test_data_file)

