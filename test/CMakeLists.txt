# Version number of the application
set (CTDET_LO_VERSION "1.5")

# Start configuration
cmake_minimum_required(VERSION 2.8)
project(ctdet_lo CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py REQUIRED)

# Default Install directory to TRIQS_ROOT
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
    message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
    set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# Macro defined in TRIQS which picks the hash of repo.
triqs_get_git_hash_of_source_dir(CTDET_GIT_HASH)

# Options for compilation that apply to everythings (C++ solver, Python interface, tests).
add_compile_options( -std=c++17 -DHAS_OPTIONAL_HEADER
        -DCTDET_GIT_HASH=${CTDET_GIT_HASH} -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
        $<$<CONFIG:Debug>:-DDEBUG_CTDET>
        )
# Common to all targets (ctint, python, tests).
include_directories(c++)

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(c++)
add_subdirectory(python/ctdet_lo)

# Tests
enable_testing()
add_subdirectory(test)

# Build the documentation
option(Build_Documentation "Build documentation" OFF)
if(${Build_Documentation})
    if (NOT ${TRIQS_WITH_DOCUMENTATION})
        message("Error: TRIQS library has not been compiled with its documentation")
    endif()
    add_subdirectory(doc)
endif()
