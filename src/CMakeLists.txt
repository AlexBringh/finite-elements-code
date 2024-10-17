include_directories(${PROJECT_SOURCE_DIR}/include)

set(SRC_FILES main.c)

add_executable(${PROJECT_NAME} ${SRC_FILES})

if (BUILD_WITH_SOLVER)
    add_subdirectory(solver)
endif ()

if (BUILD_WITH_MATH)
    add_subdirectory(math)
endif ()

if (BUILD_ALL)
    add_subdirectory(solver)
    add_subdirectory(math)
endif ()
