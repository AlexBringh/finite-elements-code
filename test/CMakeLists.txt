include_directories(${PROJECT_SOURCE_DIR}/include)

enable_testing()

# Find all .c test files in the test directory
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

# Create a test executable and add the test for each .c file
foreach(TEST_SRC ${TEST_SOURCES})
    # Get the name of the test source file without the extension
    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    
    # Create a test executable
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} femlibrary)

    # Add the test case, explicitly specifying the path to the executable
    add_test(NAME ${TEST_NAME} COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${TEST_NAME})
endforeach()