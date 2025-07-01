cmake_minimum_required(VERSION 3.10)

project(MDTools VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Werror -pthread -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE}")

configure_file(MDToolsConfig.h.in MDToolsConfig.h)

set(SOURCES rdf.cpp tools.cpp settings.cpp system.cpp pbc.cpp)
add_executable(MDTools ${SOURCES})

list(APPEND PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

target_include_directories(MDTools PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_INCLUDE_DIR}"
                           "${EXTRA_INCLUDES}"
                           )

