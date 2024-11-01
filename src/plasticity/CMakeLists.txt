# Define where to look for header files
include_directories(${PROJECT_SOURCE_DIR}/include)

# Collect all source files in this directory
file(GLOB PLASTICITY_SOURCES "*.c")

# Add teh source files to the main executable
target_sources(${PROJECT_NAME} PRIVATE ${PLASTICITY_SOURCES})