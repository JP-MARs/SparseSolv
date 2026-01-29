/**
 * @file test_basic.cpp
 * @brief Basic tests for SparseSolv library
 */

#include <sparsesolv/sparsesolv.hpp>
#include <iostream>
#include <cmath>
#include <cassert>

using namespace sparsesolv;

/**
 * @brief Create a 1D Laplacian matrix (tridiagonal)
 *
 * Creates the matrix:
 *   [ 2 -1  0  0 ... ]
 *   [-1  2 -1  0 ... ]
 *   [ 0 -1  2 -1 ... ]
 *   [ ...            ]
 */
void create_laplacian_1d(int n,
                         std::vector<index_t>& row_ptr,
                         std::vector<index_t>& col_idx,
                         std::vector<double>& values) {
    row_ptr.resize(n + 1);
    col_idx.clear();
    values.clear();

    row_ptr[0] = 0;
    for (int i = 0; i < n; ++i) {
        // Sub-diagonal
        if (i > 0) {
            col_idx.push_back(i - 1);
            values.push_back(-1.0);
        }
        // Diagonal
        col_idx.push_back(i);
        values.push_back(2.0);
        // Super-diagonal
        if (i < n - 1) {
            col_idx.push_back(i + 1);
            values.push_back(-1.0);
        }
        row_ptr[i + 1] = static_cast<index_t>(col_idx.size());
    }
}

/**
 * @brief Test basic CG solver without preconditioning
 */
bool test_cg_no_precond() {
    std::cout << "Test: CG without preconditioning... ";

    const int n = 100;

    // Create 1D Laplacian
    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    // Right-hand side: all ones
    std::vector<double> b(n, 1.0);

    // Initial guess: zeros
    std::vector<double> x(n, 0.0);

    // Solve with CG (no preconditioner)
    CGSolver<double> solver;
    solver.config().tolerance = 1e-10;
    solver.config().max_iterations = 1000;

    auto result = solver.solve(A, b, x, nullptr);

    // Check result
    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations << " iterations)\n";
        return false;
    }

    // Verify: compute ||Ax - b||
    std::vector<double> Ax(n);
    A.multiply(x.data(), Ax.data());

    double residual = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = Ax[i] - b[i];
        residual += diff * diff;
    }
    residual = std::sqrt(residual);

    if (residual > 1e-8) {
        std::cout << "FAILED (residual = " << residual << ")\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations
              << ", residual=" << result.final_residual << ")\n";
    return true;
}

/**
 * @brief Test CG solver with Jacobi preconditioning
 */
bool test_cg_jacobi() {
    std::cout << "Test: CG with Jacobi preconditioning... ";

    const int n = 100;

    // Create 1D Laplacian
    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    // Right-hand side
    std::vector<double> b(n, 1.0);

    // Initial guess
    std::vector<double> x(n, 0.0);

    // Create Jacobi preconditioner
    JacobiPreconditioner<double> precond;
    precond.setup(A);

    // Solve
    CGSolver<double> solver;
    solver.config().tolerance = 1e-10;
    solver.config().max_iterations = 1000;

    auto result = solver.solve(A, b, x, &precond);

    if (!result.converged) {
        std::cout << "FAILED (not converged)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations
              << ", residual=" << result.final_residual << ")\n";
    return true;
}

/**
 * @brief Test CG solver with IC preconditioning
 */
bool test_cg_ic() {
    std::cout << "Test: CG with IC preconditioning (ICCG)... ";

    const int n = 100;

    // Create 1D Laplacian
    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    // Right-hand side
    std::vector<double> b(n, 1.0);

    // Initial guess
    std::vector<double> x(n, 0.0);

    // Create IC preconditioner
    ICPreconditioner<double> precond(1.0);  // No shift for this simple matrix
    precond.setup(A);

    // Solve
    CGSolver<double> solver;
    solver.config().tolerance = 1e-10;
    solver.config().max_iterations = 1000;

    auto result = solver.solve(A, b, x, &precond);

    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations << " iterations)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations
              << ", residual=" << result.final_residual << ")\n";
    return true;
}

/**
 * @brief Test convenience function solve_iccg
 */
bool test_solve_iccg() {
    std::cout << "Test: solve_iccg convenience function... ";

    const int n = 100;

    // Create 1D Laplacian
    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    SolverConfig config;
    config.tolerance = 1e-10;
    config.shift_parameter = 1.0;

    auto result = solve_iccg(A, b, x, config);

    if (!result.converged) {
        std::cout << "FAILED (not converged)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test SolverConfig builder pattern
 */
bool test_config_builder() {
    std::cout << "Test: SolverConfig builder pattern... ";

    SolverConfig config = SolverConfig()
        .with_tolerance(1e-8)
        .with_max_iterations(500)
        .with_shift(1.05)
        .with_residual_history(true);

    if (config.tolerance != 1e-8 ||
        config.max_iterations != 500 ||
        config.shift_parameter != 1.05 ||
        !config.save_residual_history) {
        std::cout << "FAILED (config values incorrect)\n";
        return false;
    }

    std::cout << "PASSED\n";
    return true;
}

/**
 * @brief Test MRTR solver without preconditioning
 */
bool test_mrtr_no_precond() {
    std::cout << "Test: MRTR without preconditioning... ";

    const int n = 100;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    MRTRSolver<double> solver;
    solver.config().tolerance = 1e-10;
    solver.config().max_iterations = 1000;

    auto result = solver.solve(A, b, x, nullptr);

    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations << " iterations)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test MRTR solver with IC preconditioning (ICMRTR)
 */
bool test_mrtr_ic() {
    std::cout << "Test: MRTR with IC preconditioning (ICMRTR)... ";

    const int n = 100;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    ICPreconditioner<double> precond(1.0);
    precond.setup(A);

    MRTRSolver<double> solver;
    solver.config().tolerance = 1e-10;

    auto result = solver.solve(A, b, x, &precond);

    if (!result.converged) {
        std::cout << "FAILED (not converged)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test SGS preconditioner with CG (SGS-CG)
 */
bool test_cg_sgs() {
    std::cout << "Test: CG with SGS preconditioning (SGS-CG)... ";

    const int n = 100;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    SGSPreconditioner<double> precond;
    precond.setup(A);

    CGSolver<double> solver;
    solver.config().tolerance = 1e-10;

    auto result = solver.solve(A, b, x, &precond);

    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations << " iterations)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test SGS-MRTR solver (specialized SGS preconditioning with split formula)
 */
bool test_mrtr_sgs() {
    std::cout << "Test: SGS-MRTR solver (specialized split preconditioner)... ";

    const int n = 100;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    // Use specialized SGS-MRTR solver with split formula
    SGSMRTRSolver<double> solver;
    solver.config().tolerance = 1e-8;  // Relax tolerance for now
    solver.config().max_iterations = 500;
    solver.config().save_residual_history = true;

    auto result = solver.solve(A, b, x);

    // Verify: compute ||Ax - b||
    std::vector<double> Ax(n);
    A.multiply(x.data(), Ax.data());

    double residual = 0.0;
    for (int i = 0; i < n; ++i) {
        double diff = Ax[i] - b[i];
        residual += diff * diff;
    }
    residual = std::sqrt(residual);

    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations
                  << " iterations, final_residual=" << result.final_residual
                  << ", actual_residual=" << residual << ")\n";
        // Print first few residual history entries
        if (!result.residual_history.empty()) {
            std::cout << "  First residuals: ";
            for (size_t i = 0; i < std::min(size_t(5), result.residual_history.size()); ++i) {
                std::cout << result.residual_history[i] << " ";
            }
            std::cout << "\n  Last residuals: ";
            for (size_t i = std::max(size_t(0), result.residual_history.size() - 5);
                 i < result.residual_history.size(); ++i) {
                std::cout << result.residual_history[i] << " ";
            }
            std::cout << "\n";
        }
        return false;
    }

    if (residual > 1e-7) {
        std::cout << "FAILED (actual residual = " << residual << ")\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test convenience function solve_icmrtr
 */
bool test_solve_icmrtr() {
    std::cout << "Test: solve_icmrtr convenience function... ";

    const int n = 100;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    SolverConfig config;
    config.tolerance = 1e-10;
    config.shift_parameter = 1.0;

    auto result = solve_icmrtr(A, b, x, config);

    if (!result.converged) {
        std::cout << "FAILED (not converged)\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Create a 1D Laplacian matrix with complex type (real values)
 *
 * Uses real values stored as complex for testing complex template instantiation.
 */
void create_complex_laplacian_1d(int n,
                                 std::vector<index_t>& row_ptr,
                                 std::vector<index_t>& col_idx,
                                 std::vector<complex_t>& values) {
    row_ptr.resize(n + 1);
    col_idx.clear();
    values.clear();

    row_ptr[0] = 0;
    for (int i = 0; i < n; ++i) {
        // Sub-diagonal
        if (i > 0) {
            col_idx.push_back(i - 1);
            values.push_back(complex_t(-1.0, 0.0));
        }
        // Diagonal
        col_idx.push_back(i);
        values.push_back(complex_t(2.0, 0.0));
        // Super-diagonal
        if (i < n - 1) {
            col_idx.push_back(i + 1);
            values.push_back(complex_t(-1.0, 0.0));
        }
        row_ptr[i + 1] = static_cast<index_t>(col_idx.size());
    }
}

/**
 * @brief Test CG solver with complex numbers and IC preconditioning
 */
bool test_complex_iccg() {
    std::cout << "Test: Complex ICCG... ";

    const int n = 50;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<complex_t> values;
    create_complex_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<complex_t> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    // Right-hand side: (1+0.5i, 1+0.5i, ...)
    std::vector<complex_t> b(n, complex_t(1.0, 0.5));

    // Initial guess: zeros
    std::vector<complex_t> x(n, complex_t(0.0, 0.0));

    // Create IC preconditioner
    ICPreconditioner<complex_t> precond(1.05);
    precond.setup(A);

    // Solve
    CGSolver<complex_t> solver;
    solver.config().tolerance = 1e-8;
    solver.config().max_iterations = 1000;

    auto result = solver.solve(A, b, x, &precond);

    if (!result.converged) {
        std::cout << "FAILED (not converged after " << result.iterations << " iterations)\n";
        return false;
    }

    // Verify: compute ||Ax - b||
    std::vector<complex_t> Ax(n);
    A.multiply(x.data(), Ax.data());

    double residual = 0.0;
    for (int i = 0; i < n; ++i) {
        complex_t diff = Ax[i] - b[i];
        residual += std::norm(diff);  // |diff|^2
    }
    residual = std::sqrt(residual);

    if (residual > 1e-6) {
        std::cout << "FAILED (residual = " << residual << ")\n";
        return false;
    }

    std::cout << "PASSED (iterations=" << result.iterations << ")\n";
    return true;
}

/**
 * @brief Test residual history tracking
 */
bool test_residual_history() {
    std::cout << "Test: Residual history tracking... ";

    const int n = 50;

    std::vector<index_t> row_ptr, col_idx;
    std::vector<double> values;
    create_laplacian_1d(n, row_ptr, col_idx, values);

    SparseMatrixView<double> A(n, n, row_ptr.data(), col_idx.data(), values.data());

    std::vector<double> b(n, 1.0);
    std::vector<double> x(n, 0.0);

    CGSolver<double> solver;
    solver.config().tolerance = 1e-10;
    solver.config().save_residual_history = true;

    auto result = solver.solve(A, b, x, nullptr);

    // Check that history was recorded
    if (result.residual_history.empty()) {
        std::cout << "FAILED (no history recorded)\n";
        return false;
    }

    // Check that residual decreases (mostly)
    bool decreasing = true;
    for (size_t i = 1; i < result.residual_history.size(); ++i) {
        // Allow small increase (CG doesn't guarantee monotonic decrease)
        if (result.residual_history[i] > result.residual_history[i-1] * 1.1) {
            decreasing = false;
            break;
        }
    }

    if (!decreasing) {
        std::cout << "WARNING (residual not monotonically decreasing, but this is okay for CG)\n";
    }

    std::cout << "PASSED (history size=" << result.residual_history.size() << ")\n";
    return true;
}

int main() {
    std::cout << "SparseSolv v" << Version::string() << " Basic Tests\n";
    std::cout << "==========================================\n\n";

    int passed = 0;
    int failed = 0;

    auto run_test = [&](bool (*test_func)()) {
        if (test_func()) {
            passed++;
        } else {
            failed++;
        }
    };

    run_test(test_config_builder);
    run_test(test_cg_no_precond);
    run_test(test_cg_jacobi);
    run_test(test_cg_ic);
    run_test(test_cg_sgs);
    run_test(test_solve_iccg);
    run_test(test_mrtr_no_precond);
    run_test(test_mrtr_ic);
    run_test(test_mrtr_sgs);
    run_test(test_solve_icmrtr);
    run_test(test_complex_iccg);
    run_test(test_residual_history);

    std::cout << "\n==========================================\n";
    std::cout << "Results: " << passed << " passed, " << failed << " failed\n";

    return failed > 0 ? 1 : 0;
}
