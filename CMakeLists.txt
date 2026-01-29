cmake_minimum_required(VERSION 3.16)
project(SparseSolv VERSION 2.0.0 LANGUAGES CXX)

# C++17 required for structured bindings, etc.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(SPARSESOLV_BUILD_TESTS "Build tests" ON)
option(SPARSESOLV_BUILD_PYTHON "Build Python bindings" ON)
option(SPARSESOLV_USE_OPENMP "Use OpenMP for parallelization" ON)
option(SPARSESOLV_USE_NGSOLVE "Build with NGSolve integration" ON)
option(SPARSESOLV_BUILD_LEGACY "Build legacy SparseSolv library (deprecated)" OFF)

# Find OpenMP
if(SPARSESOLV_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    else()
        message(WARNING "OpenMP not found, disabling parallel support")
        set(SPARSESOLV_USE_OPENMP OFF)
    endif()
endif()

# Auto-detect NGSolve if not explicitly disabled
# First try to find NGSolve, then enable the option if found
if(NOT SPARSESOLV_USE_NGSOLVE)
    # Try to auto-detect NGSolve
    if(DEFINED ENV{NGSOLVE_DIR})
        set(NGSOLVE_DIR $ENV{NGSOLVE_DIR})
        list(APPEND CMAKE_PREFIX_PATH ${NGSOLVE_DIR})
    endif()
    find_package(NGSolve CONFIG QUIET)
    if(NGSolve_FOUND)
        message(STATUS "NGSolve auto-detected, enabling integration")
        set(SPARSESOLV_USE_NGSOLVE ON)
    endif()
endif()

# Find NGSolve (if explicitly enabled or auto-detected)
if(SPARSESOLV_USE_NGSOLVE)
    # Try to find NGSolve using its config
    find_package(NGSolve CONFIG)
    if(NGSolve_FOUND)
        message(STATUS "NGSolve found: ${NGSOLVE_INSTALL_DIR}")
        message(STATUS "NGSolve version: ${NGSOLVE_VERSION}")
    else()
        # Try manual search using NGSOLVE_DIR
        if(DEFINED ENV{NGSOLVE_DIR})
            set(NGSOLVE_DIR $ENV{NGSOLVE_DIR})
        endif()
        if(NGSOLVE_DIR)
            list(APPEND CMAKE_PREFIX_PATH ${NGSOLVE_DIR})
            find_package(NGSolve CONFIG REQUIRED)
            message(STATUS "NGSolve found at: ${NGSOLVE_DIR}")
        else()
            message(FATAL_ERROR "NGSolve not found. Set NGSOLVE_DIR or install NGSolve.")
        endif()
    endif()
endif()

# Header-only library (new interface)
add_library(sparsesolv INTERFACE)
target_include_directories(sparsesolv INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(sparsesolv INTERFACE cxx_std_17)

if(SPARSESOLV_USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(sparsesolv INTERFACE OpenMP::OpenMP_CXX)
endif()

# NGSolve integration
if(SPARSESOLV_USE_NGSOLVE AND NGSolve_FOUND)
    target_compile_definitions(sparsesolv INTERFACE SPARSESOLV_USE_NGSOLVE)
    # Use NGSolve's include directory
    target_include_directories(sparsesolv INTERFACE
        ${NGSOLVE_INCLUDE_DIR}
    )
    # Link to ngsolve target (provides all necessary libraries on Windows)
    target_link_libraries(sparsesolv INTERFACE ngsolve)
endif()

# Alias for consistent naming
add_library(SparseSolv::sparsesolv ALIAS sparsesolv)

# Legacy library (for backwards compatibility - deprecated)
if(SPARSESOLV_BUILD_LEGACY)
    # Set variables expected by legacy CMakeLists.txt
    set(SparseSolvBase_BINARY_DIR ${CMAKE_BINARY_DIR})
    set(SparseSolvBase_SOURCE_DIR ${CMAKE_SOURCE_DIR}/SparseSolv)
    add_subdirectory(SparseSolv)
endif()

# Tests
if(SPARSESOLV_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Python bindings
if(SPARSESOLV_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS sparsesolv
    EXPORT SparseSolvTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/sparsesolv
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT SparseSolvTargets
    FILE SparseSolvTargets.cmake
    NAMESPACE SparseSolv::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SparseSolv
)

# Package config
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/SparseSolvConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SparseSolvConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SparseSolv
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SparseSolvConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SparseSolvConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SparseSolvConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SparseSolv
)
