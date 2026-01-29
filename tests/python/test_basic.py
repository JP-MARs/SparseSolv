"""
Basic tests for SparseSolv Python bindings
"""

import pytest
import numpy as np
from scipy.sparse import diags, csr_matrix

# Import will fail if module not built, skip tests in that case
pytest.importorskip("sparsesolv")
import sparsesolv


def create_laplacian_1d(n):
    """Create 1D Laplacian matrix (tridiagonal)"""
    diagonals = [
        np.ones(n) * 2,      # main diagonal
        np.ones(n-1) * -1,   # sub-diagonal
        np.ones(n-1) * -1,   # super-diagonal
    ]
    return diags(diagonals, [0, -1, 1], format='csr')


def create_laplacian_2d(n):
    """Create 2D Laplacian matrix (5-point stencil)"""
    N = n * n
    diagonals = [
        np.ones(N) * 4,       # main diagonal
        np.ones(N-1) * -1,    # ±1
        np.ones(N-n) * -1,    # ±n
    ]
    offsets = [0, 1, n]
    A = diags(diagonals, offsets, format='csr')
    A = A + A.T - diags([A.diagonal()], [0])
    return A


class TestSolverConfig:
    """Tests for SolverConfig"""

    def test_default_values(self):
        config = sparsesolv.SolverConfig()
        assert config.tolerance == pytest.approx(1e-10)
        assert config.max_iterations == 1000
        assert config.shift_parameter == pytest.approx(1.05)

    def test_set_values(self):
        config = sparsesolv.SolverConfig()
        config.tolerance = 1e-8
        config.max_iterations = 500
        config.shift_parameter = 1.1

        assert config.tolerance == pytest.approx(1e-8)
        assert config.max_iterations == 500
        assert config.shift_parameter == pytest.approx(1.1)

    def test_repr(self):
        config = sparsesolv.SolverConfig()
        repr_str = repr(config)
        assert "SolverConfig" in repr_str
        assert "tol=" in repr_str


class TestSolverResult:
    """Tests for SolverResult"""

    def test_result_properties(self):
        A = create_laplacian_1d(10)
        b = np.ones(10)
        x, result = sparsesolv.solve_iccg(A, b)

        assert hasattr(result, 'converged')
        assert hasattr(result, 'iterations')
        assert hasattr(result, 'final_residual')

    def test_result_bool(self):
        A = create_laplacian_1d(10)
        b = np.ones(10)
        x, result = sparsesolv.solve_iccg(A, b)

        # SolverResult should be truthy if converged
        assert bool(result) == result.converged


class TestICCG:
    """Tests for ICCG solver"""

    def test_solve_1d_laplacian(self):
        n = 100
        A = create_laplacian_1d(n)
        b = np.ones(n)

        x, result = sparsesolv.solve_iccg(A, b)

        assert result.converged
        assert result.iterations < 100

        # Verify solution
        residual = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
        assert residual < 1e-8

    def test_solve_2d_laplacian(self):
        n = 20
        A = create_laplacian_2d(n)
        b = np.ones(n * n)

        config = sparsesolv.SolverConfig()
        config.tolerance = 1e-8

        x, result = sparsesolv.solve_iccg(A, b, config)

        assert result.converged

        residual = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
        assert residual < 1e-7

    def test_shift_parameter(self):
        n = 50
        A = create_laplacian_1d(n)
        b = np.ones(n)

        # Test with different shift parameters
        for shift in [1.0, 1.05, 1.1]:
            config = sparsesolv.SolverConfig()
            config.shift_parameter = shift
            config.tolerance = 1e-10

            x, result = sparsesolv.solve_iccg(A, b, config)
            assert result.converged


class TestICMRTR:
    """Tests for ICMRTR solver"""

    def test_solve_1d_laplacian(self):
        n = 100
        A = create_laplacian_1d(n)
        b = np.ones(n)

        x, result = sparsesolv.solve_icmrtr(A, b)

        assert result.converged

        residual = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
        assert residual < 1e-8

    def test_solve_2d_laplacian(self):
        n = 20
        A = create_laplacian_2d(n)
        b = np.ones(n * n)

        config = sparsesolv.SolverConfig()
        config.tolerance = 1e-8

        x, result = sparsesolv.solve_icmrtr(A, b, config)

        assert result.converged


class TestSGSMRTR:
    """Tests for SGS-MRTR solver"""

    def test_solve_1d_laplacian(self):
        n = 100
        A = create_laplacian_1d(n)
        b = np.ones(n)

        x, result = sparsesolv.solve_sgsmrtr(A, b)

        assert result.converged

        residual = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
        assert residual < 1e-8


class TestSolveScipy:
    """Tests for solve_scipy function"""

    def test_method_selection(self):
        n = 50
        A = create_laplacian_1d(n)
        b = np.ones(n)

        for method in ["ICCG", "ICMRTR", "SGSMRTR", "CG"]:
            x, result = sparsesolv.solve_scipy(A, b, method)
            # All methods should converge for this simple problem
            assert result.converged, f"Method {method} failed to converge"

    def test_invalid_method(self):
        n = 10
        A = create_laplacian_1d(n)
        b = np.ones(n)

        with pytest.raises(RuntimeError, match="Unknown method"):
            sparsesolv.solve_scipy(A, b, "InvalidMethod")

    def test_csc_matrix_auto_convert(self):
        """CSC matrix should be auto-converted to CSR"""
        n = 50
        A_csr = create_laplacian_1d(n)
        A_csc = A_csr.tocsc()
        b = np.ones(n)

        x, result = sparsesolv.solve_scipy(A_csc, b, "ICCG")
        assert result.converged


class TestResidualHistory:
    """Tests for residual history tracking"""

    def test_history_saved(self):
        n = 50
        A = create_laplacian_1d(n)
        b = np.ones(n)

        config = sparsesolv.SolverConfig()
        config.save_residual_history = True

        x, result = sparsesolv.solve_iccg(A, b, config)

        assert len(result.residual_history) > 0
        assert len(result.residual_history) == result.iterations + 1  # includes initial

    def test_history_not_saved_by_default(self):
        n = 50
        A = create_laplacian_1d(n)
        b = np.ones(n)

        x, result = sparsesolv.solve_iccg(A, b)

        # History should be empty when not requested
        assert len(result.residual_history) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
