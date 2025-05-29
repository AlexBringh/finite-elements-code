cmake_minimum_required(VERSION 3.22)
project(FEM VERSION 1.0.0 LANGUAGES C)

# Set the C Standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_DIR}/lib)

option(BUILD_WITH_SOLVER "Build with solver module" OFF)
option(BUILD_WITH_MATH "Build with math module" OFF)
option(BUILD_WITH_UTILS "Build with utils module" OFF)
option(BUILD_WITH_DATA "Build with data module" OFF)
option(BUILD_WITH_ELASTICITY "Build with elasticity module" OFF)
option(BUILD_WITH_PLASTICITY "Build with plasticity module" OFF)
option(BUILD_ALL "Build all the modules" ON)

add_subdirectory(src)
add_subdirectory(test)

# Configuring the install destinations
install(TARGETS FEM EXPORT FEMlibraryTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

# Install the public headers
install(DIRECTORY include/ DESTINATION include)

# Create and install Config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file("FEMlibraryConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_file(FEMlibraryConfig.cmake.in "FEMlibraryConfig.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/FEMlibraryConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/FEMlibraryConfigVersion.cmake"
        DESTINATION lib/cmake/FEMlibrary)

# Export the package
export(EXPORT FEMlibraryTargets FILE "${CMAKE_CURRENT_BINARY_DIR}/FEMlibraryTargets.cmake" NAMESPACE FEMlibrary::)
install(EXPORT FEMlibraryTargets FILE FEMlibraryTargets.cmake NAMESPACE FEMlibrary:: DESTINATION lib/cmake/FEMlibrary)

# Install CMake configuration files
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FEMlibraryConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FEMlibraryConfigVersion.cmake"
    DESTINATION lib/cmake/FEMlibrary
)

# Export the CMake targets
install(EXPORT FEMlibraryTargets
    FILE FEMlibraryTargets.cmake
    NAMESPACE FEMlibrary::
    DESTINATION lib/cmake/FEMlibrary
)

set(SOURCE_SAMPLE_DIR "${CMAKE_SOURCE_DIR}/sample")
set(BUILD_SAMPLE_DIR "${CMAKE_BINARY_DIR}/sample")

add_custom_target(copy_sample_dir ALL
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${BUILD_SAMPLE_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${SOURCE_SAMPLE_DIR} ${BUILD_SAMPLE_DIR}
    COMMENT "Refreshing sample directory inside the build folder"
)
