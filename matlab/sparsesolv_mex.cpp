/*
 * SPARSESOLV_MEX - MATLAB MEX Interface for SparseSolv Library
 *
 * Provides MATLAB interface to SparseSolv's iterative solvers:
 *   - ICCG (Incomplete Cholesky Conjugate Gradient) for real symmetric matrices
 *   - ICCOCG (IC + Conjugate Orthogonal CG) for complex symmetric matrices
 *   - ICMRTR (IC + MRTR) for real symmetric matrices
 *   - SGSMRTR (Symmetric Gauss-Seidel + MRTR) for real symmetric matrices
 *
 * SYNTAX:
 *   [x, flag, relres, iter, resvec] = sparsesolv_mex(solver_type, vals, col_ind, row_ptr, b, tol, max_iter, accel, x0, options)
 *
 * INPUTS:
 *   solver_type - String: 'iccg', 'iccocg', 'icmrtr', 'sgsmrtr'
 *   vals        - Non-zero values (lower triangular part only for symmetric)
 *   col_ind     - Column indices (0-based)
 *   row_ptr     - Row pointers (0-based)
 *   b           - Right-hand side vector
 *   tol         - Convergence tolerance
 *   max_iter    - Maximum iterations
 *   accel       - Acceleration (shift) parameter (negative for auto-determination)
 *   x0          - Initial guess (can be empty [])
 *   options     - Structure with optional parameters:
 *                 .diag_scaling - Enable diagonal scaling (default: false)
 *                 .save_best    - Save best result for divergence recovery (default: true)
 *                 .diverge_type - Divergence detection type (default: 1)
 *                 .diverge_val  - Divergence threshold multiplier (default: 10.0)
 *                 .diverge_count- Divergence count threshold (default: 10)
 *
 * OUTPUTS:
 *   x       - Solution vector
 *   flag    - Convergence flag (0=converged, 1=not converged)
 *   relres  - Final relative residual
 *   iter    - Number of iterations
 *   resvec  - Residual history
 *
 * For COMSOL LiveLink compatibility, use the wrapper function sparsesolv.m
 *
 * Author: SparseSolv MEX Interface
 * Based on: SRLfem SparseSolv Library
 */

#include "mex.h"
#include "matrix.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <complex>
#include <string>

// Include SparseSolv headers
#include "../SparseSolv/SparseMat.hpp"
#include "../SparseSolv/SparseMatC.hpp"
#include "../SparseSolv/MatSolvers.hpp"

using namespace SRLfem;

// Options structure
struct SolverOptions {
    bool diag_scaling;
    bool save_best;
    bool save_residual_log;
    int diverge_type;
    double diverge_val;
    int diverge_count;

    SolverOptions() :
        diag_scaling(false),
        save_best(true),
        save_residual_log(true),
        diverge_type(1),
        diverge_val(10.0),
        diverge_count(10) {}
};

// Parse options structure from MATLAB
void parseOptions(const mxArray* opts, SolverOptions& options) {
    if (opts == nullptr || mxIsEmpty(opts)) return;
    if (!mxIsStruct(opts)) {
        mexWarnMsgIdAndTxt("SparseSolv:invalidOptions", "options must be a struct");
        return;
    }

    mxArray* field;

    field = mxGetField(opts, 0, "diag_scaling");
    if (field != nullptr) {
        if (mxIsLogical(field)) {
            options.diag_scaling = mxIsLogicalScalarTrue(field);
        } else {
            options.diag_scaling = (mxGetScalar(field) != 0.0);
        }
    }

    field = mxGetField(opts, 0, "save_best");
    if (field != nullptr) {
        if (mxIsLogical(field)) {
            options.save_best = mxIsLogicalScalarTrue(field);
        } else {
            options.save_best = (mxGetScalar(field) != 0.0);
        }
    }

    field = mxGetField(opts, 0, "diverge_type");
    if (field != nullptr) {
        options.diverge_type = static_cast<int>(mxGetScalar(field));
    }

    field = mxGetField(opts, 0, "diverge_val");
    if (field != nullptr) {
        options.diverge_val = mxGetScalar(field);
    }

    field = mxGetField(opts, 0, "diverge_count");
    if (field != nullptr) {
        options.diverge_count = static_cast<int>(mxGetScalar(field));
    }
}

// Build SparseMat from CSR data (real)
SparseMat buildSparseMat(mwSize n, mwSize nnz,
                          const double* vals,
                          const mwIndex* col_idx,
                          const mwIndex* row_ptr) {
    // Convert to SparseSolv format
    std::vector<slv_int> rows, cols;
    std::vector<double> values;

    for (mwSize i = 0; i < n; i++) {
        for (mwIndex j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            // Lower triangular: row >= col
            slv_int row = static_cast<slv_int>(i);
            slv_int col = static_cast<slv_int>(col_idx[j]);

            rows.push_back(row);
            cols.push_back(col);
            values.push_back(vals[j]);

            // Add symmetric part (upper triangular) if not diagonal
            if (row != col) {
                rows.push_back(col);
                cols.push_back(row);
                values.push_back(vals[j]);
            }
        }
    }

    SparseMat mat(values.size(), rows, cols, values);
    mat.fix(true);  // Make square matrix
    return mat;
}

// Build SparseMatC from CSR data (complex)
SparseMatC buildSparseMatC(mwSize n, mwSize nnz,
                            const double* vals_real,
                            const double* vals_imag,
                            const mwIndex* col_idx,
                            const mwIndex* row_ptr) {
    std::vector<slv_int> rows, cols;
    std::vector<dcomplex> values;

    for (mwSize i = 0; i < n; i++) {
        for (mwIndex j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            slv_int row = static_cast<slv_int>(i);
            slv_int col = static_cast<slv_int>(col_idx[j]);

            dcomplex val(vals_real[j], vals_imag ? vals_imag[j] : 0.0);

            rows.push_back(row);
            cols.push_back(col);
            values.push_back(val);

            if (row != col) {
                rows.push_back(col);
                cols.push_back(row);
                values.push_back(val);  // Symmetric (not conjugate)
            }
        }
    }

    SparseMatC mat(values.size(), rows, cols, values);
    mat.fix(true);
    return mat;
}

// MEX gateway function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check input arguments
    if (nrhs < 8) {
        mexErrMsgIdAndTxt("SparseSolv:invalidNumInputs",
            "Usage: [x, flag, relres, iter, resvec] = sparsesolv_mex(solver_type, vals, col_ind, row_ptr, b, tol, max_iter, accel, [x0], [options])");
    }

    // Get solver type
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "solver_type must be a string");
    }
    char solver_type_str[64];
    mxGetString(prhs[0], solver_type_str, 64);
    std::string solver_type(solver_type_str);

    // Get matrix data in CSR format
    if (!mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "vals must be double");
    }
    const double* vals = mxGetPr(prhs[1]);
    const double* vals_imag = mxIsComplex(prhs[1]) ? mxGetPi(prhs[1]) : nullptr;
    mwSize nnz = mxGetNumberOfElements(prhs[1]);

    if (!mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "col_ind must be double");
    }
    double* col_ind_double = mxGetPr(prhs[2]);

    if (!mxIsDouble(prhs[3])) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "row_ptr must be double");
    }
    double* row_ptr_double = mxGetPr(prhs[3]);
    mwSize n = mxGetNumberOfElements(prhs[3]) - 1;

    // Convert indices to mwIndex
    std::vector<mwIndex> col_idx(nnz);
    std::vector<mwIndex> row_ptr(n + 1);

    // Detect if 1-based indexing
    bool is_one_based = (row_ptr_double[0] == 1.0);

    for (mwSize i = 0; i < nnz; i++) {
        col_idx[i] = static_cast<mwIndex>(col_ind_double[i] - (is_one_based ? 1 : 0));
    }
    for (mwSize i = 0; i <= n; i++) {
        row_ptr[i] = static_cast<mwIndex>(row_ptr_double[i] - (is_one_based ? 1 : 0));
    }

    // Get right-hand side vector b
    if (!mxIsDouble(prhs[4])) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "b must be double");
    }
    const double* b = mxGetPr(prhs[4]);
    const double* b_imag = mxIsComplex(prhs[4]) ? mxGetPi(prhs[4]) : nullptr;
    mwSize b_size = mxGetNumberOfElements(prhs[4]);

    if (b_size != n) {
        mexErrMsgIdAndTxt("SparseSolv:invalidInput", "b must have same size as matrix dimension");
    }

    // Get solver parameters
    double tol = mxGetScalar(prhs[5]);
    int max_iter = static_cast<int>(mxGetScalar(prhs[6]));
    double accel = mxGetScalar(prhs[7]);

    // Get initial guess (optional)
    double* x0 = nullptr;
    double* x0_imag = nullptr;
    if (nrhs > 8 && !mxIsEmpty(prhs[8])) {
        x0 = mxGetPr(prhs[8]);
        x0_imag = mxIsComplex(prhs[8]) ? mxGetPi(prhs[8]) : nullptr;
    }

    // Get options (optional)
    SolverOptions options;
    if (nrhs > 9 && !mxIsEmpty(prhs[9])) {
        parseOptions(prhs[9], options);
    }

    // Print info
    mexPrintf("SparseSolv: solver=%s, n=%d, nnz=%d, tol=%.2e, max_iter=%d, accel=%.3f\n",
              solver_type.c_str(), (int)n, (int)nnz, tol, max_iter, accel);
    mexPrintf("Options: diag_scaling=%d, save_best=%d, diverge_type=%d\n",
              options.diag_scaling, options.save_best, options.diverge_type);

    // Create solver
    MatSolvers solver;
    solver.setDiagScale(options.diag_scaling);
    solver.setSaveBest(options.save_best);
    solver.setSaveLog(true);
    solver.setDirvegeType(options.diverge_type);
    solver.setBadDivVal(options.diverge_val);
    solver.setBadDivCount(options.diverge_count);

    bool is_complex = (vals_imag != nullptr || b_imag != nullptr ||
                       solver_type == "iccocg");
    bool converged = false;

    if (!is_complex) {
        // Real solver
        SparseMat matA = buildSparseMat(n, nnz, vals, col_idx.data(), row_ptr.data());

        // Allocate output
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
        double* x = mxGetPr(plhs[0]);

        // Set initial guess
        if (x0) {
            std::memcpy(x, x0, n * sizeof(double));
        } else {
            std::memset(x, 0, n * sizeof(double));
        }

        // Solve
        if (solver_type == "iccg") {
            converged = solver.solveICCG(n, tol, max_iter, accel, matA, b, x, (x0 == nullptr));
        } else if (solver_type == "icmrtr") {
            converged = solver.solveICMRTR(n, tol, max_iter, accel, matA, b, x, (x0 == nullptr));
        } else if (solver_type == "sgsmrtr") {
            converged = solver.solveSGSMRTR(n, tol, max_iter, matA, b, x, (x0 == nullptr));
        } else {
            mexErrMsgIdAndTxt("SparseSolv:invalidSolver",
                "Unknown solver type: %s. Use 'iccg', 'iccocg', 'icmrtr', or 'sgsmrtr'",
                solver_type.c_str());
        }
    } else {
        // Complex solver
        SparseMatC matA = buildSparseMatC(n, nnz, vals, vals_imag, col_idx.data(), row_ptr.data());

        // Allocate complex output
        plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
        double* x_real = mxGetPr(plhs[0]);
        double* x_imag = mxGetPi(plhs[0]);

        // Build complex vectors
        std::vector<dcomplex> b_vec(n);
        std::vector<dcomplex> x_vec(n);

        for (mwSize i = 0; i < n; i++) {
            b_vec[i] = dcomplex(b[i], b_imag ? b_imag[i] : 0.0);
            if (x0) {
                x_vec[i] = dcomplex(x0[i], x0_imag ? x0_imag[i] : 0.0);
            } else {
                x_vec[i] = dcomplex(0.0, 0.0);
            }
        }

        // Solve
        if (solver_type == "iccocg" || solver_type == "iccg") {
            converged = solver.solveICCG(n, tol, max_iter, accel, matA, b_vec, x_vec, (x0 == nullptr));
        } else if (solver_type == "icmrtr") {
            converged = solver.solveICMRTR(n, tol, max_iter, accel, matA, b_vec, x_vec, (x0 == nullptr));
        } else if (solver_type == "sgsmrtr") {
            converged = solver.solveSGSMRTR(n, tol, max_iter, matA, b_vec, x_vec, (x0 == nullptr));
        } else {
            mexErrMsgIdAndTxt("SparseSolv:invalidSolver",
                "Unknown solver type for complex: %s", solver_type.c_str());
        }

        // Copy results
        for (mwSize i = 0; i < n; i++) {
            x_real[i] = x_vec[i].real();
            x_imag[i] = x_vec[i].imag();
        }
    }

    // Get residual log
    std::vector<double> resvec;
    solver.getResidualLog(resvec);
    int iter = resvec.size();

    // Calculate final relative residual
    double relres = (iter > 0) ? resvec.back() : 1.0;

    // Output: flag
    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleScalar(converged ? 0.0 : 1.0);
    }

    // Output: relres
    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleScalar(relres);
    }

    // Output: iter
    if (nlhs >= 4) {
        plhs[3] = mxCreateDoubleScalar(static_cast<double>(iter));
    }

    // Output: resvec
    if (nlhs >= 5) {
        plhs[4] = mxCreateDoubleMatrix(iter, 1, mxREAL);
        double* resvec_out = mxGetPr(plhs[4]);
        for (int i = 0; i < iter; i++) {
            resvec_out[i] = resvec[i];
        }
    }

    mexPrintf("SparseSolv: %s, iter=%d, relres=%.2e\n",
              converged ? "Converged" : "Not converged", iter, relres);
}
