/**
 * @file bindings.cpp
 * @brief pybind11 bindings for SparseSolv library
 *
 * This module provides Python bindings for the SparseSolv library,
 * enabling use of ICCG, MRTR, and other solvers from Python.
 *
 * Example usage:
 * @code
 * import sparsesolv
 * import numpy as np
 * from scipy.sparse import csr_matrix
 *
 * # Create sparse matrix
 * A = csr_matrix(...)
 * b = np.array(...)
 *
 * # Configure solver
 * config = sparsesolv.SolverConfig()
 * config.tolerance = 1e-10
 *
 * # Solve
 * x, result = sparsesolv.solve_iccg(A, b, config)
 * print(f"Converged: {result.converged}, iterations: {result.iterations}")
 * @endcode
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <sparsesolv/sparsesolv.hpp>

#ifdef SPARSESOLV_USE_NGSOLVE
#include <sparsesolv/adapters/ngsolve_adapter.hpp>
#endif

namespace py = pybind11;
using namespace sparsesolv;

/**
 * @brief Extract CSR data from scipy.sparse.csr_matrix
 */
struct CSRMatrix {
    py::array_t<index_t> indptr;
    py::array_t<index_t> indices;
    py::array_t<double> data;
    index_t rows;
    index_t cols;

    static CSRMatrix from_scipy(py::object scipy_mat) {
        CSRMatrix csr;

        // Get shape
        py::tuple shape = scipy_mat.attr("shape").cast<py::tuple>();
        csr.rows = shape[0].cast<index_t>();
        csr.cols = shape[1].cast<index_t>();

        // Convert to CSR if needed
        py::object csr_mat = scipy_mat;
        if (!py::hasattr(scipy_mat, "indptr")) {
            csr_mat = scipy_mat.attr("tocsr")();
        }

        // Extract arrays
        csr.indptr = csr_mat.attr("indptr").cast<py::array_t<index_t>>();
        csr.indices = csr_mat.attr("indices").cast<py::array_t<index_t>>();
        csr.data = csr_mat.attr("data").cast<py::array_t<double>>();

        return csr;
    }

    SparseMatrixView<double> view() const {
        return SparseMatrixView<double>(
            rows, cols,
            indptr.data(), indices.data(), data.data()
        );
    }
};

/**
 * @brief Solve Ax=b using scipy sparse matrix
 */
std::pair<py::array_t<double>, SolverResult> solve_scipy(
    py::object scipy_matrix,
    py::array_t<double> b,
    const std::string& method,
    const SolverConfig& config
) {
    // Extract CSR data
    CSRMatrix csr = CSRMatrix::from_scipy(scipy_matrix);
    auto A = csr.view();

    // Get b data
    py::buffer_info b_info = b.request();
    if (b_info.ndim != 1) {
        throw std::runtime_error("b must be 1-dimensional");
    }
    index_t n = static_cast<index_t>(b_info.size);
    const double* b_ptr = static_cast<const double*>(b_info.ptr);

    // Allocate solution
    py::array_t<double> x(n);
    py::buffer_info x_info = x.request();
    double* x_ptr = static_cast<double*>(x_info.ptr);
    std::fill(x_ptr, x_ptr + n, 0.0);

    // Solve based on method
    SolverResult result;
    if (method == "ICCG" || method == "iccg") {
        result = solve_iccg(A, b_ptr, x_ptr, n, config);
    } else if (method == "ICMRTR" || method == "icmrtr") {
        result = solve_icmrtr(A, b_ptr, x_ptr, n, config);
    } else if (method == "SGSMRTR" || method == "sgsmrtr") {
        result = solve_sgsmrtr(A, b_ptr, x_ptr, n, config);
    } else if (method == "CG" || method == "cg") {
        // CG without preconditioning
        CGSolver<double> solver;
        solver.set_config(config);
        result = solver.solve(A, b_ptr, x_ptr, n, nullptr);
    } else {
        throw std::runtime_error("Unknown method: " + method + ". Use ICCG, ICMRTR, SGSMRTR, or CG.");
    }

    return std::make_pair(x, result);
}

/**
 * @brief Solve Ax=b using raw CSR arrays
 */
SolverResult solve_csr(
    py::array_t<index_t> row_ptr,
    py::array_t<index_t> col_idx,
    py::array_t<double> values,
    py::array_t<double> b,
    py::array_t<double> x,
    const std::string& method,
    const SolverConfig& config
) {
    // Get array info
    py::buffer_info rp_info = row_ptr.request();
    py::buffer_info ci_info = col_idx.request();
    py::buffer_info v_info = values.request();
    py::buffer_info b_info = b.request();
    py::buffer_info x_info = x.request();

    index_t n = static_cast<index_t>(b_info.size);

    // Create matrix view
    SparseMatrixView<double> A(
        n, n,
        static_cast<const index_t*>(rp_info.ptr),
        static_cast<const index_t*>(ci_info.ptr),
        static_cast<const double*>(v_info.ptr)
    );

    const double* b_ptr = static_cast<const double*>(b_info.ptr);
    double* x_ptr = static_cast<double*>(x_info.ptr);

    // Solve
    if (method == "ICCG" || method == "iccg") {
        return solve_iccg(A, b_ptr, x_ptr, n, config);
    } else if (method == "ICMRTR" || method == "icmrtr") {
        return solve_icmrtr(A, b_ptr, x_ptr, n, config);
    } else if (method == "SGSMRTR" || method == "sgsmrtr") {
        return solve_sgsmrtr(A, b_ptr, x_ptr, n, config);
    } else {
        CGSolver<double> solver;
        solver.set_config(config);
        return solver.solve(A, b_ptr, x_ptr, n, nullptr);
    }
}

PYBIND11_MODULE(sparsesolv, m) {
    m.doc() = R"pbdoc(
        SparseSolv: Sparse matrix solvers for finite element analysis

        This module provides iterative solvers for sparse linear systems:
        - ICCG: Incomplete Cholesky Conjugate Gradient
        - ICMRTR: Incomplete Cholesky MRTR
        - SGSMRTR: Symmetric Gauss-Seidel MRTR

        Example:
            >>> import sparsesolv
            >>> from scipy.sparse import csr_matrix
            >>> import numpy as np
            >>>
            >>> # Create tridiagonal matrix
            >>> n = 100
            >>> A = csr_matrix(...)
            >>> b = np.ones(n)
            >>>
            >>> # Solve
            >>> x, result = sparsesolv.solve_scipy(A, b, "ICCG")
            >>> print(f"Converged: {result.converged}")
    )pbdoc";

    // Version info
    m.attr("__version__") = Version::string();

    // NormType enum
    py::enum_<NormType>(m, "NormType")
        .value("RHS", NormType::RHS, "Normalize by right-hand side norm")
        .value("InitialResidual", NormType::InitialResidual, "Normalize by initial residual")
        .value("Custom", NormType::Custom, "Use custom normalization value")
        .export_values();

    // DivergenceCheck enum
    py::enum_<DivergenceCheck>(m, "DivergenceCheck")
        .value("NONE", DivergenceCheck::None, "No divergence check")
        .value("StagnationCount", DivergenceCheck::StagnationCount, "Stop on stagnation")
        .export_values();

    // SolverConfig class
    py::class_<SolverConfig>(m, "SolverConfig", "Configuration for iterative solvers")
        .def(py::init<>())
        .def_readwrite("tolerance", &SolverConfig::tolerance,
            "Relative tolerance for convergence (default: 1e-10)")
        .def_readwrite("max_iterations", &SolverConfig::max_iterations,
            "Maximum number of iterations (default: 1000)")
        .def_readwrite("abs_tolerance", &SolverConfig::abs_tolerance,
            "Absolute tolerance threshold")
        .def_readwrite("norm_type", &SolverConfig::norm_type,
            "Normalization type for convergence check")
        .def_readwrite("custom_norm", &SolverConfig::custom_norm,
            "Custom normalization value")
        .def_readwrite("shift_parameter", &SolverConfig::shift_parameter,
            "Shift parameter for IC decomposition (default: 1.05)")
        .def_readwrite("diagonal_scaling", &SolverConfig::diagonal_scaling,
            "Enable diagonal scaling")
        .def_readwrite("divergence_check", &SolverConfig::divergence_check,
            "Divergence detection strategy")
        .def_readwrite("divergence_threshold", &SolverConfig::divergence_threshold,
            "Threshold for divergence detection")
        .def_readwrite("divergence_count", &SolverConfig::divergence_count,
            "Count for divergence detection")
        .def_readwrite("save_best_result", &SolverConfig::save_best_result,
            "Save best result for non-converging solves")
        .def_readwrite("save_residual_history", &SolverConfig::save_residual_history,
            "Save residual history")
        .def_readwrite("num_threads", &SolverConfig::num_threads,
            "Number of threads (0 = auto)")
        .def("__repr__", [](const SolverConfig& c) {
            return "<SolverConfig tol=" + std::to_string(c.tolerance) +
                   " max_iter=" + std::to_string(c.max_iterations) +
                   " shift=" + std::to_string(c.shift_parameter) + ">";
        });

    // SolverResult class
    py::class_<SolverResult>(m, "SolverResult", "Result of an iterative solve")
        .def(py::init<>())
        .def_readonly("converged", &SolverResult::converged,
            "Whether the solver converged")
        .def_readonly("iterations", &SolverResult::iterations,
            "Number of iterations performed")
        .def_readonly("final_residual", &SolverResult::final_residual,
            "Final relative residual")
        .def_readonly("residual_history", &SolverResult::residual_history,
            "Residual at each iteration (if saved)")
        .def("__bool__", [](const SolverResult& r) { return r.converged; })
        .def("__repr__", [](const SolverResult& r) {
            return "<SolverResult converged=" + std::string(r.converged ? "True" : "False") +
                   " iter=" + std::to_string(r.iterations) +
                   " residual=" + std::to_string(r.final_residual) + ">";
        });

    // Main solve functions
    m.def("solve_scipy", &solve_scipy,
        R"pbdoc(
            Solve Ax=b using a scipy sparse matrix.

            Parameters
            ----------
            A : scipy.sparse matrix
                System matrix (will be converted to CSR if needed)
            b : numpy.ndarray
                Right-hand side vector
            method : str
                Solver method: "ICCG", "ICMRTR", "SGSMRTR", or "CG"
            config : SolverConfig, optional
                Solver configuration

            Returns
            -------
            tuple
                (x, result) where x is the solution and result is SolverResult
        )pbdoc",
        py::arg("A"),
        py::arg("b"),
        py::arg("method") = "ICCG",
        py::arg("config") = SolverConfig());

    m.def("solve_csr", &solve_csr,
        R"pbdoc(
            Solve Ax=b using raw CSR arrays.

            Parameters
            ----------
            row_ptr : numpy.ndarray
                CSR row pointers
            col_idx : numpy.ndarray
                CSR column indices
            values : numpy.ndarray
                CSR values
            b : numpy.ndarray
                Right-hand side vector
            x : numpy.ndarray
                Solution vector (initial guess, modified in place)
            method : str
                Solver method
            config : SolverConfig
                Solver configuration

            Returns
            -------
            SolverResult
                Solver result
        )pbdoc",
        py::arg("row_ptr"),
        py::arg("col_idx"),
        py::arg("values"),
        py::arg("b"),
        py::arg("x"),
        py::arg("method") = "ICCG",
        py::arg("config") = SolverConfig());

    // Convenience functions
    m.def("solve_iccg", [](py::object A, py::array_t<double> b, const SolverConfig& config) {
        return solve_scipy(A, b, "ICCG", config);
    }, "Solve Ax=b using ICCG", py::arg("A"), py::arg("b"), py::arg("config") = SolverConfig());

    m.def("solve_icmrtr", [](py::object A, py::array_t<double> b, const SolverConfig& config) {
        return solve_scipy(A, b, "ICMRTR", config);
    }, "Solve Ax=b using ICMRTR", py::arg("A"), py::arg("b"), py::arg("config") = SolverConfig());

    m.def("solve_sgsmrtr", [](py::object A, py::array_t<double> b, const SolverConfig& config) {
        return solve_scipy(A, b, "SGSMRTR", config);
    }, "Solve Ax=b using SGSMRTR", py::arg("A"), py::arg("b"), py::arg("config") = SolverConfig());

#ifdef SPARSESOLV_USE_NGSOLVE
    // NGSolve integration - solve using NGSolve BaseMatrix and BaseVector
    m.def("solve_ngsolve", [](py::object mat, py::object rhs, py::object sol,
                               const std::string& method, const SolverConfig& config) -> SolverResult {
        // Cast to NGSolve types
        auto& ng_mat = py::cast<ngla::SparseMatrix<double>&>(mat);
        auto& ng_rhs = py::cast<ngla::BaseVector&>(rhs);
        auto& ng_sol = py::cast<ngla::BaseVector&>(sol);

        return ngsolve_adapter::solve_ngsolve(ng_mat, ng_rhs, ng_sol, method, config);
    },
    R"pbdoc(
        Solve Ax=b using NGSolve BaseMatrix and BaseVector.

        Parameters
        ----------
        mat : ngsolve.BaseMatrix
            System matrix (must be SparseMatrix)
        rhs : ngsolve.BaseVector
            Right-hand side vector
        sol : ngsolve.BaseVector
            Solution vector (modified in place)
        method : str
            Solver method: "ICCG", "ICMRTR", "SGSMRTR", or "CG"
        config : SolverConfig
            Solver configuration

        Returns
        -------
        SolverResult
            Solver result

        Example
        -------
        >>> from ngsolve import *
        >>> import sparsesolv
        >>>
        >>> # After assembling a.mat and f.vec
        >>> config = sparsesolv.SolverConfig()
        >>> config.tolerance = 1e-10
        >>> result = sparsesolv.solve_ngsolve(a.mat, f.vec, gfu.vec, "ICCG", config)
    )pbdoc",
    py::arg("mat"), py::arg("rhs"), py::arg("sol"),
    py::arg("method") = "ICCG", py::arg("config") = SolverConfig());

    // IC Preconditioner class for use with NGSolve's CGSolver
    py::class_<ngsolve_adapter::ICPreconditionerNGS<double>,
               std::shared_ptr<ngsolve_adapter::ICPreconditionerNGS<double>>,
               ngla::BaseMatrix>(m, "ICPreconditioner",
        R"pbdoc(
            Incomplete Cholesky preconditioner for NGSolve.

            This preconditioner can be used with NGSolve's built-in solvers
            like CGSolver, GMRESSolver, etc.

            Example
            -------
            >>> from ngsolve import *
            >>> from ngsolve.krylovspace import CGSolver
            >>> import sparsesolv
            >>>
            >>> # Create and assemble bilinear form
            >>> a.Assemble()
            >>>
            >>> # Create IC preconditioner
            >>> pre = sparsesolv.ICPreconditioner(a.mat, shift=1.05)
            >>> pre.Update()
            >>>
            >>> # Use with NGSolve's CGSolver
            >>> inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
            >>> gfu.vec.data = inv * f.vec
        )pbdoc")
        .def(py::init([](py::object mat, double shift) {
            auto sp_mat = py::cast<std::shared_ptr<ngla::SparseMatrix<double>>>(mat);
            return std::make_shared<ngsolve_adapter::ICPreconditionerNGS<double>>(sp_mat, shift);
        }), py::arg("mat"), py::arg("shift") = 1.05,
            "Create IC preconditioner from NGSolve SparseMatrix")
        .def("Update", &ngsolve_adapter::ICPreconditionerNGS<double>::Update,
            "Update preconditioner (recompute factorization after matrix change)")
        .def_property("shift",
            &ngsolve_adapter::ICPreconditionerNGS<double>::GetShift,
            &ngsolve_adapter::ICPreconditionerNGS<double>::SetShift,
            "Shift parameter for IC decomposition");

    // SGS Preconditioner class for use with NGSolve's CGSolver
    py::class_<ngsolve_adapter::SGSPreconditionerNGS<double>,
               std::shared_ptr<ngsolve_adapter::SGSPreconditionerNGS<double>>,
               ngla::BaseMatrix>(m, "SGSPreconditioner",
        R"pbdoc(
            Symmetric Gauss-Seidel preconditioner for NGSolve.

            This preconditioner can be used with NGSolve's built-in solvers
            like CGSolver, GMRESSolver, etc.

            Example
            -------
            >>> pre = sparsesolv.SGSPreconditioner(a.mat)
            >>> pre.Update()
            >>> inv = CGSolver(a.mat, pre, printrates=True)
        )pbdoc")
        .def(py::init([](py::object mat) {
            auto sp_mat = py::cast<std::shared_ptr<ngla::SparseMatrix<double>>>(mat);
            return std::make_shared<ngsolve_adapter::SGSPreconditionerNGS<double>>(sp_mat);
        }), py::arg("mat"),
            "Create SGS preconditioner from NGSolve SparseMatrix")
        .def("Update", &ngsolve_adapter::SGSPreconditionerNGS<double>::Update,
            "Update preconditioner (recompute after matrix change)");

    // Attribute to indicate NGSolve support is available
    m.attr("has_ngsolve") = true;
#else
    m.attr("has_ngsolve") = false;
#endif
}
