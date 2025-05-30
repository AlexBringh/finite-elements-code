include_directories(${PROJECT_SOURCE_DIR}/include)

set(SRC_FILES main.c)

add_executable(${PROJECT_NAME} ${SRC_FILES})

add_subdirectory(elements)
add_subdirectory(postprocessing)

if (BUILD_WITH_SOLVER)
    add_subdirectory(solver)
endif ()

if (BUILD_WITH_MATH)
    add_subdirectory(math)
endif ()

if (BUILD_WITH_UTILS)
    add_subdirectory(utils)
endif ()

if (BUILD_WITH_DATA)
    add_subdirectory(data)
endif ()

if (BUILD_WITH_PLASTICITY)
    add_subdirectory(plasticity)
endif()

if (BUILD_ALL)
    add_subdirectory(solver)
    add_subdirectory(math)
    add_subdirectory(utils)
    add_subdirectory(data)
    add_subdirectory(plasticity)
endif ()


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


