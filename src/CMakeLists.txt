include_directories(${PROJECT_SOURCE_DIR}/include)

set(SRC_FILES main.c)

add_executable(${PROJECT_NAME} ${SRC_FILES})

add_subdirectory(elements) # REMOVE: Remove this once restructuring is done
add_subdirectory(postprocessing)
add_subdirectory(mesh)
add_subdirectory(elasticity)
add_subdirectory(materials)
add_subdirectory(cfd)
add_subdirectory(solver)
add_subdirectory(math)
add_subdirectory(utils)
add_subdirectory(data)
add_subdirectory(plasticity)
add_subdirectory(codes)

# Make library for tests
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.c")
add_library(femlibrary STATIC ${LIB_SOURCES})
# Create a library from the collected source files

# Specify the include directories for the library target
target_include_directories(femlibrary PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Path to the sample folder and output location
set(SAMPLE_DIR "${CMAKE_SOURCE_DIR}/sample")
set(OUTPUT_SAMPLE_DIR "${CMAKE_BINARY_DIR}/sample")


