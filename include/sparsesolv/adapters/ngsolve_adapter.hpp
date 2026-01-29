/**
 * @file ngsolve_adapter.hpp
 * @brief Adapter for using SparseSolv with NGSolve
 *
 * This adapter allows SparseSolv preconditioners and solvers to be used
 * seamlessly with NGSolve's linear algebra infrastructure.
 *
 * Example usage in Python (with NGSolve):
 * @code
 * from ngsolve import *
 * import sparsesolv
 *
 * mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))
 * fes = H1(mesh, order=3)
 * u, v = fes.TnT()
 *
 * a = BilinearForm(fes)
 * a += grad(u)*grad(v)*dx
 * a.Assemble()
 *
 * f = LinearForm(fes)
 * f += 1*v*dx
 * f.Assemble()
 *
 * gfu = GridFunction(fes)
 *
 * # Use SparseSolv ICCG
 * config = sparsesolv.SolverConfig()
 * config.tolerance = 1e-10
 * sparsesolv.solve_ngsolve(a.mat, f.vec, gfu.vec, "ICCG", config)
 * @endcode
 */

#ifndef SPARSESOLV_ADAPTERS_NGSOLVE_ADAPTER_HPP
#define SPARSESOLV_ADAPTERS_NGSOLVE_ADAPTER_HPP

#include "../sparsesolv.hpp"

// NGSolve headers (only included if NGSolve is available)
#ifdef SPARSESOLV_USE_NGSOLVE

#include <comp.hpp>
#include <la.hpp>

namespace sparsesolv {
namespace ngsolve_adapter {

// ============================================================================
// Vector Data Access Helpers
// ============================================================================

/**
 * @brief Helper to extract raw data pointer from NGSolve vectors
 *
 * Supports:
 * - VVector<SCAL>: Direct access to contiguous data
 * - BlockVector: Flattened view (all blocks must be contiguous)
 */
template<typename SCAL>
const SCAL* get_vector_data(const ngla::BaseVector& vec) {
    // FV<SCAL>() provides a FlatVector view into the data
    // This works for VVector and BlockVector (flattened)
    try {
        return vec.FV<SCAL>().Data();
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("SparseSolv: Cannot access vector data. ") +
            "Supported types: VVector, BlockVector (contiguous). " +
            "Error: " + e.what());
    }
}

template<typename SCAL>
SCAL* get_vector_data(ngla::BaseVector& vec) {
    try {
        return vec.FV<SCAL>().Data();
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("SparseSolv: Cannot access vector data. ") +
            "Supported types: VVector, BlockVector (contiguous). " +
            "Error: " + e.what());
    }
}

/**
 * @brief Get the size of a vector
 */
template<typename SCAL>
index_t get_vector_size(const ngla::BaseVector& vec) {
    return static_cast<index_t>(vec.FV<SCAL>().Size());
}

/**
 * @brief Apply preconditioner with proper vector type handling
 *
 * This helper handles different vector types and applies the preconditioner.
 */
template<typename SCAL, typename Precond>
void apply_preconditioner(
    const Precond& precond,
    const ngla::BaseVector& x,
    ngla::BaseVector& y,
    index_t size
) {
    const SCAL* x_data = get_vector_data<SCAL>(x);
    SCAL* y_data = get_vector_data<SCAL>(y);

    // Verify sizes match
    index_t x_size = get_vector_size<SCAL>(x);
    index_t y_size = get_vector_size<SCAL>(y);

    if (x_size != size || y_size != size) {
        throw std::runtime_error(
            "SparseSolv: Vector size mismatch. Expected " + std::to_string(size) +
            ", got x=" + std::to_string(x_size) + ", y=" + std::to_string(y_size));
    }

    precond->apply(x_data, y_data, size);
}

/**
 * @brief Apply preconditioner with MultAdd semantics: y += s * M^{-1} * x
 */
template<typename SCAL, typename Precond>
void apply_preconditioner_multadd(
    const Precond& precond,
    SCAL s,
    const ngla::BaseVector& x,
    ngla::BaseVector& y,
    index_t size
) {
    const SCAL* x_data = get_vector_data<SCAL>(x);
    SCAL* y_data = get_vector_data<SCAL>(y);

    // Verify sizes
    index_t x_size = get_vector_size<SCAL>(x);
    index_t y_size = get_vector_size<SCAL>(y);

    if (x_size != size || y_size != size) {
        throw std::runtime_error(
            "SparseSolv: Vector size mismatch in MultAdd. Expected " + std::to_string(size) +
            ", got x=" + std::to_string(x_size) + ", y=" + std::to_string(y_size));
    }

    // temp = M^{-1} * x
    std::vector<SCAL> temp(size);
    precond->apply(x_data, temp.data(), size);

    // y = y + s * temp
    for (index_t i = 0; i < size; ++i) {
        y_data[i] += s * temp[i];
    }
}

// ============================================================================
// Matrix Wrapper
// ============================================================================

/**
 * @brief Create a SparseMatrixView from NGSolve's SparseMatrix
 *
 * This creates a zero-copy view into NGSolve's sparse matrix data.
 *
 * @tparam SCAL Scalar type
 * @param mat NGSolve sparse matrix
 * @return SparseMatrixView
 */
template<typename SCAL = double>
SparseMatrixView<SCAL> wrap_ngsolve_matrix(const ngla::SparseMatrix<SCAL>& mat) {
    return SparseMatrixView<SCAL>(
        mat.Height(),
        mat.Width(),
        mat.GetRowIndices().Data(),
        mat.GetColIndices().Data(),
        mat.GetData().Data()
    );
}

/**
 * @brief Wrapper to make SparseSolv preconditioner usable as NGSolve preconditioner
 *
 * This class wraps a SparseSolv Preconditioner to implement NGSolve's
 * BaseMatrix interface, allowing it to be used with NGSolve solvers.
 *
 * Supports:
 * - VVector<SCAL>: Standard dense vectors
 * - BlockVector: Block vectors (flattened access)
 */
template<typename SCAL = double>
class NGSolvePreconditionerWrapper : public ngla::BaseMatrix {
public:
    NGSolvePreconditionerWrapper(std::shared_ptr<Preconditioner<SCAL>> precond,
                                  const ngla::SparseMatrix<SCAL>& mat)
        : precond_(precond)
        , height_(mat.Height())
        , width_(mat.Width())
    {
        auto view = wrap_ngsolve_matrix(mat);
        precond_->setup(view);
    }

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner<SCAL>(precond_, x, y, height_);
    }

    void MultAdd(SCAL s, const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner_multadd<SCAL>(precond_, s, x, y, height_);
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    ngla::AutoVector CreateRowVector() const override {
        return ngla::AutoVector(std::make_shared<ngla::VVector<SCAL>>(width_));
    }

    ngla::AutoVector CreateColVector() const override {
        return ngla::AutoVector(std::make_shared<ngla::VVector<SCAL>>(height_));
    }

private:
    std::shared_ptr<Preconditioner<SCAL>> precond_;
    index_t height_;
    index_t width_;
};

/**
 * @brief Solve NGSolve linear system using SparseSolv
 *
 * Supports various NGSolve vector types:
 * - VVector<SCAL>: Standard dense vectors
 * - BlockVector: Block vectors (flattened access)
 *
 * @tparam SCAL Scalar type
 * @param mat NGSolve sparse matrix
 * @param rhs Right-hand side vector
 * @param sol Solution vector (modified in place)
 * @param method Solver method: "ICCG", "ICMRTR", "SGSMRTR", "CG"
 * @param config Solver configuration
 * @return SolverResult
 */
template<typename SCAL = double>
SolverResult solve_ngsolve(
    const ngla::SparseMatrix<SCAL>& mat,
    const ngla::BaseVector& rhs,
    ngla::BaseVector& sol,
    const std::string& method = "ICCG",
    const SolverConfig& config = SolverConfig()
) {
    auto A = wrap_ngsolve_matrix(mat);

    const SCAL* b = get_vector_data<SCAL>(rhs);
    SCAL* x = get_vector_data<SCAL>(sol);
    index_t n = static_cast<index_t>(mat.Height());

    // Verify vector sizes
    index_t rhs_size = get_vector_size<SCAL>(rhs);
    index_t sol_size = get_vector_size<SCAL>(sol);
    if (rhs_size != n || sol_size != n) {
        throw std::runtime_error(
            "SparseSolv: Vector size mismatch. Matrix size=" + std::to_string(n) +
            ", rhs size=" + std::to_string(rhs_size) +
            ", sol size=" + std::to_string(sol_size));
    }

    if (method == "ICCG" || method == "iccg") {
        return solve_iccg(A, b, x, n, config);
    } else if (method == "ICMRTR" || method == "icmrtr") {
        return solve_icmrtr(A, b, x, n, config);
    } else if (method == "SGSMRTR" || method == "sgsmrtr") {
        return solve_sgsmrtr(A, b, x, n, config);
    } else {
        CGSolver<SCAL> solver;
        solver.set_config(config);
        return solver.solve(A, b, x, n, nullptr);
    }
}

/**
 * @brief Create IC preconditioner for NGSolve
 */
template<typename SCAL = double>
std::shared_ptr<ngla::BaseMatrix> create_ic_preconditioner(
    const ngla::SparseMatrix<SCAL>& mat,
    double shift = 1.05
) {
    auto precond = std::make_shared<ICPreconditioner<SCAL>>(shift);
    return std::make_shared<NGSolvePreconditionerWrapper<SCAL>>(precond, mat);
}

/**
 * @brief Create SGS preconditioner for NGSolve
 */
template<typename SCAL = double>
std::shared_ptr<ngla::BaseMatrix> create_sgs_preconditioner(
    const ngla::SparseMatrix<SCAL>& mat
) {
    auto precond = std::make_shared<SGSPreconditioner<SCAL>>();
    return std::make_shared<NGSolvePreconditionerWrapper<SCAL>>(precond, mat);
}

/**
 * @brief IC Preconditioner following NGSolve's Preconditioner pattern
 *
 * This class implements NGSolve's Preconditioner interface, allowing it to be
 * used with NGSolve's built-in solvers (CGSolver, GMRESSolver, etc.)
 *
 * Supports various NGSolve vector types:
 * - VVector<SCAL>: Standard dense vectors
 * - BlockVector: Block vectors (flattened access)
 *
 * Usage in Python:
 * @code
 * from ngsolve import *
 * import sparsesolv
 *
 * # Create bilinear form and assemble
 * a = BilinearForm(fes)
 * a += grad(u)*grad(v)*dx
 * a.Assemble()
 *
 * # Create IC preconditioner
 * pre = sparsesolv.ICPreconditioner(a.mat, shift=1.05)
 * pre.Update()
 *
 * # Use with NGSolve's CGSolver
 * from ngsolve.krylovspace import CGSolver
 * inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
 * gfu.vec.data = inv * f.vec
 * @endcode
 */
template<typename SCAL = double>
class ICPreconditionerNGS : public ngla::BaseMatrix {
public:
    /**
     * @brief Construct IC preconditioner from NGSolve sparse matrix
     * @param mat The sparse matrix to precondition
     * @param shift Shift parameter for IC decomposition (default: 1.05)
     */
    ICPreconditionerNGS(std::shared_ptr<ngla::SparseMatrix<SCAL>> mat, double shift = 1.05)
        : mat_(mat)
        , shift_(shift)
        , height_(static_cast<index_t>(mat->Height()))
        , width_(static_cast<index_t>(mat->Width()))
        , precond_(std::make_shared<ICPreconditioner<SCAL>>(shift))
    {}

    /**
     * @brief Update the preconditioner (recompute factorization)
     *
     * Call this after the matrix has been assembled or modified.
     */
    void Update() {
        auto view = wrap_ngsolve_matrix(*mat_);
        precond_->setup(view);
    }

    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     *
     * Supports VVector and BlockVector.
     */
    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner<SCAL>(precond_, x, y, height_);
    }

    void MultAdd(SCAL s, const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner_multadd<SCAL>(precond_, s, x, y, height_);
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    ngla::AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    ngla::AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }

    /// Get shift parameter
    double GetShift() const { return shift_; }

    /// Set shift parameter (requires Update() to take effect)
    void SetShift(double shift) { shift_ = shift; precond_->set_shift_parameter(shift); }

private:
    std::shared_ptr<ngla::SparseMatrix<SCAL>> mat_;
    double shift_;
    index_t height_;
    index_t width_;
    std::shared_ptr<ICPreconditioner<SCAL>> precond_;
};

/**
 * @brief SGS Preconditioner following NGSolve's Preconditioner pattern
 *
 * Supports various NGSolve vector types:
 * - VVector<SCAL>: Standard dense vectors
 * - BlockVector: Block vectors (flattened access)
 *
 * Usage in Python:
 * @code
 * pre = sparsesolv.SGSPreconditioner(a.mat)
 * pre.Update()
 * inv = CGSolver(a.mat, pre, printrates=True)
 * @endcode
 */
template<typename SCAL = double>
class SGSPreconditionerNGS : public ngla::BaseMatrix {
public:
    SGSPreconditionerNGS(std::shared_ptr<ngla::SparseMatrix<SCAL>> mat)
        : mat_(mat)
        , height_(static_cast<index_t>(mat->Height()))
        , width_(static_cast<index_t>(mat->Width()))
        , precond_(std::make_shared<SGSPreconditioner<SCAL>>())
    {}

    void Update() {
        auto view = wrap_ngsolve_matrix(*mat_);
        precond_->setup(view);
    }

    /**
     * @brief Apply preconditioner: y = M^{-1} * x
     *
     * Supports VVector and BlockVector.
     */
    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner<SCAL>(precond_, x, y, height_);
    }

    void MultAdd(SCAL s, const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        apply_preconditioner_multadd<SCAL>(precond_, s, x, y, height_);
    }

    int VHeight() const override { return static_cast<int>(height_); }
    int VWidth() const override { return static_cast<int>(width_); }

    ngla::AutoVector CreateRowVector() const override {
        return mat_->CreateRowVector();
    }

    ngla::AutoVector CreateColVector() const override {
        return mat_->CreateColVector();
    }

private:
    std::shared_ptr<ngla::SparseMatrix<SCAL>> mat_;
    index_t height_;
    index_t width_;
    std::shared_ptr<SGSPreconditioner<SCAL>> precond_;
};

} // namespace ngsolve_adapter
} // namespace sparsesolv

#endif // SPARSESOLV_USE_NGSOLVE

#endif // SPARSESOLV_ADAPTERS_NGSOLVE_ADAPTER_HPP
