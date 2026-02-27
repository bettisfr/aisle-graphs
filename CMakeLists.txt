cmake_minimum_required(VERSION 3.23)

project(aisle-graphs CXX)

set(CMAKE_CXX_STANDARD 17)

option(CXX "enable C++ compilation" ON)

if (CXX)
    enable_language(CXX)
endif ()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(sources
        src/core/algorithms.cpp
        src/core/algorithms.h
        src/io/input.cpp
        src/io/input.h
        src/io/output.cpp
        src/io/output.h
        src/core/deployment.cpp
        src/core/deployment.h
        src/util/util.cpp
        src/util/util.h
        src/main.cpp
)

add_executable(${CMAKE_PROJECT_NAME} ${sources})

find_package(GUROBI REQUIRED)

if (NOT GUROBI_CXX_LIBRARY)
    message(FATAL_ERROR "GUROBI C++ library not found (GUROBI_CXX_LIBRARY).")
endif ()

target_compile_definitions(${CMAKE_PROJECT_NAME} PRIVATE WITH_GUROBI=1)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE ${GUROBI_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif ()
