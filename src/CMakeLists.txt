cmake_minimum_required(VERSION 3.10)
project(MDTools VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Optional: Find OpenMP if you're using it
find_package(OpenMP)

set(CMAKE_CXX_FLAGS_DEBUG
    "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic -Werror -pthread -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE
    "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Configure header file
configure_file(MDToolsConfig.h.in MDToolsConfig.h)

# Source files
set(SOURCES rdf.cpp tools.cpp settings.cpp system.cpp pbc.cpp main.cpp)

# Create executable
add_executable(MDTools ${SOURCES})

# Include directories
list(APPEND PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
target_include_directories(MDTools PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           "${PROJECT_INCLUDE_DIR}"
                           "${EXTRA_INCLUDES}"
                           )

# Link libraries
target_link_libraries(MDTools Eigen3::Eigen)

# Optional: Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(MDTools OpenMP::OpenMP_CXX)
endif()

# Optional: Print configuration info
message(STATUS "MDTools Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  Eigen3 found: ${Eigen3_FOUND}")
if(OpenMP_CXX_FOUND)
    message(STATUS "  OpenMP found: ${OpenMP_CXX_FOUND}")
endif()
