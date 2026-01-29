# Tests CMakeLists.txt

# Basic tests (no external dependencies)
add_executable(test_basic cpp/test_basic.cpp)
target_link_libraries(test_basic PRIVATE sparsesolv)
target_include_directories(test_basic PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_test(NAME basic_tests COMMAND test_basic)

# If Google Test is available, build more comprehensive tests
find_package(GTest QUIET)
if(GTest_FOUND)
    message(STATUS "Google Test found, building comprehensive tests")

    # TODO: Add Google Test based tests
    # add_executable(test_comprehensive
    #     cpp/test_cg_solver.cpp
    #     cpp/test_ic_preconditioner.cpp
    #     cpp/test_sparse_matrix.cpp
    # )
    # target_link_libraries(test_comprehensive PRIVATE sparsesolv GTest::gtest_main)
    # add_test(NAME comprehensive_tests COMMAND test_comprehensive)
else()
    message(STATUS "Google Test not found, only building basic tests")
endif()

# NGSolve integration tests
if(SPARSESOLV_USE_NGSOLVE)
    add_subdirectory(ngsolve)
endif()
