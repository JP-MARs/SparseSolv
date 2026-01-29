# NGSolve integration tests

if(SPARSESOLV_USE_NGSOLVE AND NGSolve_FOUND)
    find_package(Python3 COMPONENTS Development REQUIRED)

    add_executable(test_ngsolve_preconditioner test_ngsolve_preconditioner.cpp)
    target_link_libraries(test_ngsolve_preconditioner PRIVATE sparsesolv Python3::Python)
    target_include_directories(test_ngsolve_preconditioner PRIVATE ${Python3_INCLUDE_DIRS})

    add_test(NAME test_ngsolve_preconditioner COMMAND test_ngsolve_preconditioner)
endif()
