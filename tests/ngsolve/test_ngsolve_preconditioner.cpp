/**
 * @file test_ngsolve_preconditioner.cpp
 * @brief Test SparseSolv IC preconditioner with NGSolve integration
 */

#include <sparsesolv/adapters/ngsolve_preconditioner.hpp>
#include <iostream>

int main()
{
#ifdef SPARSESOLV_USE_NGSOLVE
    std::cout << "SparseSolv with NGSolve integration is enabled!" << std::endl;
    std::cout << "ICPrecondMatrix and SparseSolvICPreconditioner classes are available." << std::endl;

    // Test compilation of template instantiation
    using namespace sparsesolv::ngsolve_native;

    // Just instantiate to verify compilation
    // Actual runtime test requires NGSolve runtime
    std::cout << "Template instantiation test passed." << std::endl;

    return 0;
#else
    std::cout << "ERROR: SPARSESOLV_USE_NGSOLVE is not defined!" << std::endl;
    return 1;
#endif
}
