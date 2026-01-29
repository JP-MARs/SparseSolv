# Python bindings CMakeLists.txt

find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(sparsesolv_py bindings.cpp)
target_link_libraries(sparsesolv_py PRIVATE sparsesolv)
target_include_directories(sparsesolv_py PRIVATE ${CMAKE_SOURCE_DIR}/include)

# NGSolve integration for Python bindings
if(SPARSESOLV_USE_NGSOLVE AND NGSolve_FOUND)
    target_compile_definitions(sparsesolv_py PRIVATE SPARSESOLV_USE_NGSOLVE)
    target_include_directories(sparsesolv_py PRIVATE ${NGSOLVE_INCLUDE_DIR})
    target_link_libraries(sparsesolv_py PRIVATE ngsolve)
    message(STATUS "Python bindings: NGSolve integration enabled")
endif()

# Set output name to sparsesolv (without _py suffix)
set_target_properties(sparsesolv_py PROPERTIES
    OUTPUT_NAME "sparsesolv"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)

# Install Python module
install(TARGETS sparsesolv_py
    LIBRARY DESTINATION ${Python3_SITELIB}
)
