/**
 * @file ngsolve_preconditioner.hpp
 * @brief NGSolve-native IC preconditioner using TaskManager for parallelization
 *
 * This file provides SparseSolv preconditioners that integrate natively with
 * NGSolve's Preconditioner base class, supporting:
 * - Update() mechanism for automatic preconditioner updates
 * - TaskManager-based parallelization
 * - Registration as NGSolve preconditioner type
 *
 * Usage in Python:
 * @code
 * from ngsolve import *
 *
 * mesh = Mesh(unit_cube.GenerateMesh(maxh=0.1))
 * fes = H1(mesh, order=3)
 * u, v = fes.TnT()
 *
 * a = BilinearForm(fes)
 * a += grad(u)*grad(v)*dx
 *
 * # Use SparseSolv IC preconditioner
 * c = Preconditioner(a, "sparsesolv_ic", shift=1.05)
 * a.Assemble()
 *
 * f = LinearForm(fes)
 * f += 1*v*dx
 * f.Assemble()
 *
 * gfu = GridFunction(fes)
 * solvers.CG(mat=a.mat, rhs=f.vec, sol=gfu.vec, pre=c.mat,
 *            maxsteps=1000, tol=1e-10)
 * @endcode
 */

#ifndef SPARSESOLV_ADAPTERS_NGSOLVE_PRECONDITIONER_HPP
#define SPARSESOLV_ADAPTERS_NGSOLVE_PRECONDITIONER_HPP

#ifdef SPARSESOLV_USE_NGSOLVE

#include <comp.hpp>
#include <la.hpp>
#include <core/taskmanager.hpp>

namespace sparsesolv {
namespace ngsolve_native {

using namespace ngcomp;
using namespace ngla;
using namespace ngcore;

/**
 * @brief NGSolve-native IC preconditioner matrix
 *
 * This class stores the IC factorization result and implements
 * the BaseMatrix interface for use with NGSolve solvers.
 */
template<typename SCAL = double>
class ICPrecondMatrix : public BaseMatrix
{
public:
    ICPrecondMatrix(int aheight, double ashift = 1.05)
        : height_(aheight), shift_(ashift)
    {}

    /**
     * @brief Setup IC factorization from sparse matrix
     *
     * Uses TaskManager for parallel computation of IC factors.
     */
    void Setup(const SparseMatrix<SCAL>& mat)
    {
        const int n = mat.Height();
        height_ = n;

        // Get CSR data
        const auto& first_array = mat.GetFirstArray();
        const auto& col_indices = mat.GetColIndices();
        const auto& values = mat.GetValues();

        // Allocate storage for L factor (CSR format)
        // First pass: count lower triangular entries per row
        L_row_ptr_.SetSize(n + 1);

        ParallelFor(Range(n), [&](int i) {
            int count = 0;
            for (int k = first_array[i]; k < first_array[i + 1]; k++) {
                if (col_indices[k] <= i) {
                    count++;
                }
            }
            L_row_ptr_[i] = count;
        });

        // Cumulative sum for row pointers
        int total = 0;
        for (int i = 0; i <= n; i++) {
            int count = L_row_ptr_[i];
            L_row_ptr_[i] = total;
            total += count;
        }

        // Allocate L factor arrays
        L_col_idx_.SetSize(total);
        L_values_.SetSize(total);
        inv_diag_.SetSize(n);

        // Second pass: copy lower triangular part
        ParallelFor(Range(n), [&](int i) {
            int pos = L_row_ptr_[i];
            for (int k = first_array[i]; k < first_array[i + 1]; k++) {
                int j = col_indices[k];
                if (j <= i) {
                    L_col_idx_[pos] = j;
                    L_values_[pos] = values[k];
                    pos++;
                }
            }
        });

        // Compute IC factorization (sequential due to dependencies)
        ComputeICFactorization();

        // Compute L^T for backward substitution
        ComputeTranspose();

        // Compute level scheduling for parallel triangular solves
        ComputeLevelScheduling();

        is_setup_ = true;
    }

    void Mult(const BaseVector& x, BaseVector& y) const override
    {
        const FlatVector<SCAL> fx = x.FV<SCAL>();
        FlatVector<SCAL> fy = y.FV<SCAL>();
        Apply(fx.Data(), fy.Data());
    }

    void MultAdd(SCAL s, const BaseVector& x, BaseVector& y) const override
    {
        const FlatVector<SCAL> fx = x.FV<SCAL>();
        FlatVector<SCAL> fy = y.FV<SCAL>();

        // temp = M^{-1} * x
        Array<SCAL> temp(height_);
        Apply(fx.Data(), temp.Data());

        // y += s * temp (parallel)
        ParallelForRange(Range(height_), [&](IntRange r) {
            for (int i : r) {
                fy(i) += s * temp[i];
            }
        });
    }

    int VHeight() const override { return height_; }
    int VWidth() const override { return height_; }

    AutoVector CreateRowVector() const override
    {
        return make_unique<VVector<SCAL>>(height_);
    }

    AutoVector CreateColVector() const override
    {
        return make_unique<VVector<SCAL>>(height_);
    }

private:
    int height_;
    double shift_;
    bool is_setup_ = false;

    // L factor (CSR format)
    Array<int> L_row_ptr_;
    Array<int> L_col_idx_;
    Array<SCAL> L_values_;

    // L^T factor (CSR format)
    Array<int> Lt_row_ptr_;
    Array<int> Lt_col_idx_;
    Array<SCAL> Lt_values_;

    // D^{-1} diagonal
    Array<SCAL> inv_diag_;

    // Level scheduling data for parallel triangular solves
    Array<int> L_levels_;           // Level assignment for each row (forward)
    Array<int> L_level_ptr_;        // Start of each level in L_level_order_
    Array<int> L_level_order_;      // Row order grouped by level

    Array<int> Lt_levels_;          // Level assignment for each row (backward)
    Array<int> Lt_level_ptr_;       // Start of each level in Lt_level_order_
    Array<int> Lt_level_order_;     // Row order grouped by level

    int num_L_levels_ = 0;
    int num_Lt_levels_ = 0;

    // Temporary vector for apply
    mutable Array<SCAL> temp_;

    /**
     * @brief Compute IC factorization in-place
     *
     * Sequential algorithm due to data dependencies, but individual
     * row operations are optimized.
     */
    void ComputeICFactorization()
    {
        const int n = height_;

        for (int i = 0; i < n; i++) {
            const int row_start = L_row_ptr_[i];
            const int row_end = L_row_ptr_[i + 1];

            // Process off-diagonal elements L(i, j) where j < i
            for (int kk = row_start; kk < row_end - 1; kk++) {
                const int j = L_col_idx_[kk];
                if (j >= i) break;

                SCAL s = L_values_[kk];

                // s -= sum over k < j of: L(i,k) * L(j,k) * D(k)
                const int j_start = L_row_ptr_[j];
                const int j_end = L_row_ptr_[j + 1];

                for (int ii = row_start; ii < kk; ii++) {
                    const int k = L_col_idx_[ii];
                    if (k >= j) break;

                    // Find L(j, k) in row j
                    for (int jj = j_start; jj < j_end; jj++) {
                        if (L_col_idx_[jj] == k) {
                            // Multiply by D^{-1}[k] (stored in inv_diag_)
                            s -= L_values_[ii] * L_values_[jj] * inv_diag_[k];
                            break;
                        } else if (L_col_idx_[jj] > k) {
                            break;
                        }
                    }
                }

                L_values_[kk] = s;
            }

            // Process diagonal element
            const int diag_pos = row_end - 1;
            SCAL s = L_values_[diag_pos] * static_cast<SCAL>(shift_);

            for (int kk = row_start; kk < diag_pos; kk++) {
                const int k = L_col_idx_[kk];
                if (k >= i) break;
                // Multiply by D^{-1}[k] (stored in inv_diag_)
                s -= L_values_[kk] * L_values_[kk] * inv_diag_[k];
            }

            L_values_[diag_pos] = s;

            // Store D^{-1}(i)
            if (std::abs(s) > 1e-15) {
                inv_diag_[i] = SCAL(1) / s;
            } else {
                inv_diag_[i] = SCAL(1e15);
            }
        }
    }

    /**
     * @brief Compute L^T (transpose of L)
     *
     * Uses TaskManager for parallel transpose computation.
     */
    void ComputeTranspose()
    {
        const int n = height_;
        const int nnz = L_values_.Size();

        Lt_row_ptr_.SetSize(n + 1);
        Lt_col_idx_.SetSize(nnz);
        Lt_values_.SetSize(nnz);

        // Count entries per column of L (= per row of L^T)
        Lt_row_ptr_ = 0;
        for (int k = 0; k < nnz; k++) {
            Lt_row_ptr_[L_col_idx_[k] + 1]++;
        }

        // Cumulative sum
        for (int i = 0; i < n; i++) {
            Lt_row_ptr_[i + 1] += Lt_row_ptr_[i];
        }

        // Fill values (parallel per source row)
        Array<std::atomic<int>> counter(n);
        for (int i = 0; i < n; i++) {
            counter[i].store(0);
        }

        ParallelFor(Range(n), [&](int i) {
            for (int k = L_row_ptr_[i]; k < L_row_ptr_[i + 1]; k++) {
                int j = L_col_idx_[k];
                int pos = Lt_row_ptr_[j] + counter[j].fetch_add(1);
                Lt_col_idx_[pos] = i;
                Lt_values_[pos] = L_values_[k];
            }
        });
    }

    /**
     * @brief Compute level scheduling for parallel triangular solves
     *
     * For forward substitution (L): row i depends on rows j where L[i,j] != 0 and j < i
     * For backward substitution (L^T): row i depends on rows j where L^T[i,j] != 0 and j > i
     *
     * Rows at the same level have no dependencies among them and can be computed in parallel.
     */
    void ComputeLevelScheduling()
    {
        const int n = height_;

        // Compute levels for forward substitution (L)
        L_levels_.SetSize(n);
        L_levels_ = 0;

        for (int i = 0; i < n; i++) {
            int max_dep_level = -1;
            const int row_start = L_row_ptr_[i];
            const int row_end = L_row_ptr_[i + 1] - 1; // Exclude diagonal

            for (int k = row_start; k < row_end; k++) {
                int j = L_col_idx_[k];
                if (j < i) {
                    max_dep_level = std::max(max_dep_level, L_levels_[j]);
                }
            }
            L_levels_[i] = max_dep_level + 1;
        }

        // Find number of levels
        num_L_levels_ = 0;
        for (int i = 0; i < n; i++) {
            num_L_levels_ = std::max(num_L_levels_, L_levels_[i] + 1);
        }

        // Count rows per level
        Array<int> level_count(num_L_levels_);
        level_count = 0;
        for (int i = 0; i < n; i++) {
            level_count[L_levels_[i]]++;
        }

        // Build level pointers and order
        L_level_ptr_.SetSize(num_L_levels_ + 1);
        L_level_ptr_[0] = 0;
        for (int l = 0; l < num_L_levels_; l++) {
            L_level_ptr_[l + 1] = L_level_ptr_[l] + level_count[l];
        }

        L_level_order_.SetSize(n);
        Array<int> level_pos(num_L_levels_);
        for (int l = 0; l < num_L_levels_; l++) {
            level_pos[l] = L_level_ptr_[l];
        }
        for (int i = 0; i < n; i++) {
            int l = L_levels_[i];
            L_level_order_[level_pos[l]++] = i;
        }

        // Compute levels for backward substitution (L^T)
        // Process in reverse order since L^T is upper triangular
        Lt_levels_.SetSize(n);
        Lt_levels_ = 0;

        for (int i = n - 1; i >= 0; i--) {
            int max_dep_level = -1;
            const int row_start = Lt_row_ptr_[i] + 1; // Skip diagonal
            const int row_end = Lt_row_ptr_[i + 1];

            for (int k = row_start; k < row_end; k++) {
                int j = Lt_col_idx_[k];
                if (j > i) {
                    max_dep_level = std::max(max_dep_level, Lt_levels_[j]);
                }
            }
            Lt_levels_[i] = max_dep_level + 1;
        }

        // Find number of levels for L^T
        num_Lt_levels_ = 0;
        for (int i = 0; i < n; i++) {
            num_Lt_levels_ = std::max(num_Lt_levels_, Lt_levels_[i] + 1);
        }

        // Count rows per level
        Array<int> lt_level_count(num_Lt_levels_);
        lt_level_count = 0;
        for (int i = 0; i < n; i++) {
            lt_level_count[Lt_levels_[i]]++;
        }

        // Build level pointers and order
        Lt_level_ptr_.SetSize(num_Lt_levels_ + 1);
        Lt_level_ptr_[0] = 0;
        for (int l = 0; l < num_Lt_levels_; l++) {
            Lt_level_ptr_[l + 1] = Lt_level_ptr_[l] + lt_level_count[l];
        }

        Lt_level_order_.SetSize(n);
        Array<int> lt_level_pos(num_Lt_levels_);
        for (int l = 0; l < num_Lt_levels_; l++) {
            lt_level_pos[l] = Lt_level_ptr_[l];
        }
        for (int i = 0; i < n; i++) {
            int l = Lt_levels_[i];
            Lt_level_order_[lt_level_pos[l]++] = i;
        }
    }

    /**
     * @brief Apply IC preconditioner: y = (LDL^T)^{-1} * x
     *
     * Uses level scheduling for parallel triangular solves.
     */
    void Apply(const SCAL* x, SCAL* y) const
    {
        const int n = height_;

        // Ensure temp vector is allocated
        if (temp_.Size() != n) {
            temp_.SetSize(n);
        }

        // Forward substitution with level scheduling
        // Solve L * temp = x
        for (int l = 0; l < num_L_levels_; l++) {
            const int level_start = L_level_ptr_[l];
            const int level_end = L_level_ptr_[l + 1];

            // All rows in this level can be processed in parallel
            ParallelFor(IntRange(level_start, level_end), [&](int idx) {
                int i = L_level_order_[idx];
                SCAL s = x[i];
                const int row_start = L_row_ptr_[i];
                const int row_end = L_row_ptr_[i + 1] - 1; // Exclude diagonal

                for (int k = row_start; k < row_end; k++) {
                    s -= L_values_[k] * temp_[L_col_idx_[k]];
                }

                temp_[i] = s / L_values_[row_end];
            });
        }

        // Backward substitution with level scheduling
        // Solve L^T * y = D^{-1} * temp
        for (int l = 0; l < num_Lt_levels_; l++) {
            const int level_start = Lt_level_ptr_[l];
            const int level_end = Lt_level_ptr_[l + 1];

            // All rows in this level can be processed in parallel
            ParallelFor(IntRange(level_start, level_end), [&](int idx) {
                int i = Lt_level_order_[idx];
                SCAL s = SCAL(0);
                const int row_start = Lt_row_ptr_[i] + 1; // Skip diagonal
                const int row_end = Lt_row_ptr_[i + 1];

                for (int k = row_start; k < row_end; k++) {
                    s -= Lt_values_[k] * y[Lt_col_idx_[k]];
                }

                y[i] = s * inv_diag_[i] + temp_[i];
            });
        }
    }
};


/**
 * @brief NGSolve-native IC Preconditioner with Update() support
 *
 * This class integrates with NGSolve's preconditioner framework:
 * - Inherits from ngcomp::Preconditioner
 * - Implements Update() for automatic updates when BilinearForm changes
 * - Uses TaskManager for parallel operations
 * - Registerable as "sparsesolv_ic" preconditioner type
 */
class SparseSolvICPreconditioner : public Preconditioner
{
public:
    SparseSolvICPreconditioner(shared_ptr<BilinearForm> abfa,
                               const Flags& aflags,
                               const string& aname = "sparsesolv_ic")
        : Preconditioner(abfa, aflags, aname)
        , shift_(aflags.GetNumFlag("shift", 1.05))
    {
        // Get bilinear form
        bfa_ = abfa;
    }

    /**
     * @brief Update preconditioner when matrix changes
     *
     * Called automatically by NGSolve when BilinearForm is reassembled.
     * Uses timestamp comparison to avoid redundant updates.
     */
    virtual void Update() override
    {
        // Check if update is needed
        if (GetTimeStamp() == bfa_->GetTimeStamp()) {
            return; // Already up to date
        }

        // Update timestamp
        timestamp = bfa_->GetTimeStamp();

        // Get system matrix
        const BaseMatrix& basemat = bfa_->GetMatrix();

        // Try to cast to SparseMatrix<double>
        auto sparsemat = dynamic_cast<const SparseMatrix<double>*>(&basemat);
        if (!sparsemat) {
            throw Exception("SparseSolvICPreconditioner requires SparseMatrix<double>");
        }

        // Create or update preconditioner matrix
        if (!precond_matrix_ || precond_matrix_->VHeight() != sparsemat->Height()) {
            precond_matrix_ = make_shared<ICPrecondMatrix<double>>(
                sparsemat->Height(), shift_);
        }

        // Setup IC factorization
        precond_matrix_->Setup(*sparsemat);

        // Track memory
        GetMemoryTracer().Track(*precond_matrix_, "SparseSolvIC");

        // Optional testing
        if (test) Test();
        if (timing) Timing();
    }

    virtual const BaseMatrix& GetMatrix() const override
    {
        if (!precond_matrix_) {
            ThrowPreconditionerNotReady();
        }
        return *precond_matrix_;
    }

    virtual shared_ptr<BaseMatrix> GetMatrixPtr() const
    {
        return precond_matrix_;
    }

    virtual const BaseMatrix& GetAMatrix() const override
    {
        return bfa_->GetMatrix();
    }

    virtual const char* ClassName() const override
    {
        return "SparseSolvICPreconditioner";
    }

    /// Get shift parameter
    double GetShift() const { return shift_; }

    /// Set shift parameter (requires Update() after change)
    void SetShift(double shift) { shift_ = shift; }

private:
    shared_ptr<BilinearForm> bfa_;
    shared_ptr<ICPrecondMatrix<double>> precond_matrix_;
    double shift_;
};


/**
 * @brief NGSolve-native SGS Preconditioner with Update() support
 */
class SparseSolvSGSPreconditioner : public Preconditioner
{
public:
    SparseSolvSGSPreconditioner(shared_ptr<BilinearForm> abfa,
                                const Flags& aflags,
                                const string& aname = "sparsesolv_sgs")
        : Preconditioner(abfa, aflags, aname)
    {
        bfa_ = abfa;
    }

    virtual void Update() override;
    virtual const BaseMatrix& GetMatrix() const override;
    virtual const BaseMatrix& GetAMatrix() const override;
    virtual const char* ClassName() const override { return "SparseSolvSGSPreconditioner"; }

private:
    shared_ptr<BilinearForm> bfa_;
    shared_ptr<BaseMatrix> precond_matrix_;
};


// Registration helper (call in module initialization)
inline void RegisterSparseSolvPreconditioners()
{
    GetPreconditionerClasses().AddPreconditioner(
        "sparsesolv_ic",
        [](shared_ptr<BilinearForm> bfa, const Flags& flags, const string name)
            -> shared_ptr<Preconditioner>
        {
            return make_shared<SparseSolvICPreconditioner>(bfa, flags, name);
        }
    );
}

} // namespace ngsolve_native
} // namespace sparsesolv

#endif // SPARSESOLV_USE_NGSOLVE

#endif // SPARSESOLV_ADAPTERS_NGSOLVE_PRECONDITIONER_HPP
