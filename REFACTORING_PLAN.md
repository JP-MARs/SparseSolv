# SparseSolv リファクタリング計画

## 目標

1. **NGSolve統合**: NGSolveのPreconditioner/Solverとして組み込み可能な設計
2. **pybind11最適化**: Pythonから使いやすいAPI設計
3. **ソルバー共通化**: 各ソルバーの共通部分を抽出し、差分のみを実装

---

## 現状の問題点

### 1. コード重複
- `SparseSolv/` と `SparseSolvPy/` で全ファイルが重複
- 各ソルバー（ICCG, MRTR, SGS-MRTR）で似たようなコードが繰り返し

### 2. NGSolve非互換
- 独自の疎行列フォーマット（Eigen::SparseMatrix）
- NGSolveのBaseMatrix/BaseVectorと互換性なし

### 3. API設計の問題
- 24個のオーバーロード関数
- 設定が分散（7つのsetter関数）
- 日本語変数名（gyo, retu）

---

## 新設計

### アーキテクチャ

```
sparsesolv/
├── core/                      # コアライブラリ（ヘッダオンリー可）
│   ├── sparse_matrix.hpp      # 疎行列インターフェース
│   ├── preconditioner.hpp     # 前処理基底クラス
│   ├── solver_base.hpp        # ソルバー基底クラス
│   ├── solver_config.hpp      # 設定構造体
│   ├── ic_decomposition.hpp   # IC分解
│   ├── cg_solver.hpp          # CG法実装
│   └── mrtr_solver.hpp        # MRTR法実装
├── preconditioners/           # 具体的な前処理
│   ├── ic_preconditioner.hpp  # 不完全コレスキー
│   ├── sgs_preconditioner.hpp # SGS前処理
│   └── abmc_preconditioner.hpp # ABMC順序付きIC
├── adapters/                  # 外部ライブラリアダプタ
│   ├── eigen_adapter.hpp      # Eigen対応
│   └── ngsolve_adapter.hpp    # NGSolve対応（重要）
├── python/                    # pybind11バインディング
│   └── bindings.cpp
└── CMakeLists.txt
```

### クラス設計

#### 1. SolverConfig（設定の一元管理）

```cpp
namespace sparsesolv {

struct SolverConfig {
    // 収束判定
    double tolerance = 1e-10;
    int max_iterations = 1000;

    // 前処理設定
    double shift_parameter = 1.05;  // 旧: accera（加速係数）
    bool diagonal_scaling = false;

    // 収束判定の正規化
    enum class NormType { RHS, InitialResidual, Custom };
    NormType norm_type = NormType::RHS;
    double custom_norm = 1.0;

    // 発散検出
    bool detect_divergence = true;
    double divergence_threshold = 1000.0;
    int divergence_count = 100;

    // 結果保存
    bool save_best_result = true;
    bool save_residual_history = false;
};

} // namespace sparsesolv
```

#### 2. Preconditioner基底クラス（NGSolve互換）

```cpp
namespace sparsesolv {

// NGSolveのBasePreconditionerと同様のインターフェース
template<typename Scalar = double>
class Preconditioner {
public:
    virtual ~Preconditioner() = default;

    // 前処理の適用: y = M^{-1} * x
    virtual void Apply(const Scalar* x, Scalar* y, int size) const = 0;

    // Eigen対応
    virtual void Apply(const Eigen::VectorX<Scalar>& x,
                       Eigen::VectorX<Scalar>& y) const {
        Apply(x.data(), y.data(), x.size());
    }

    // 行列からの構築
    virtual void Setup(const SparseMatrixView<Scalar>& matrix) = 0;

    // 前処理名
    virtual std::string Name() const = 0;
};

} // namespace sparsesolv
```

#### 3. IC Preconditioner（不完全コレスキー）

```cpp
namespace sparsesolv {

template<typename Scalar = double>
class ICPreconditioner : public Preconditioner<Scalar> {
public:
    explicit ICPreconditioner(double shift = 1.05)
        : shift_parameter_(shift) {}

    void Setup(const SparseMatrixView<Scalar>& A) override {
        // IC分解: A ≈ LDL^T
        ComputeICDecomposition(A, L_, D_, shift_parameter_);
    }

    void Apply(const Scalar* x, Scalar* y, int size) const override {
        // y = (LDL^T)^{-1} * x
        // 1. 前進代入: L * z = x
        // 2. 対角スケーリング: w = D^{-1} * z
        // 3. 後退代入: L^T * y = w
        ForwardSubstitution(L_, x, temp_.data(), size);
        DiagonalScale(D_.data(), temp_.data(), size);
        BackwardSubstitution(Lt_, temp_.data(), y, size);
    }

    std::string Name() const override { return "IC"; }

private:
    double shift_parameter_;
    SparseMatrix<Scalar> L_;   // 下三角
    SparseMatrix<Scalar> Lt_;  // 上三角（L^T）
    std::vector<Scalar> D_;    // 対角
    mutable std::vector<Scalar> temp_;
};

} // namespace sparsesolv
```

#### 4. IterativeSolver基底クラス

```cpp
namespace sparsesolv {

template<typename Scalar = double>
class IterativeSolver {
public:
    struct Result {
        bool converged;
        int iterations;
        double final_residual;
        std::vector<double> residual_history;  // optional
    };

    virtual ~IterativeSolver() = default;

    // メイン解法インターフェース
    Result Solve(const SparseMatrixView<Scalar>& A,
                 const Scalar* b,
                 Scalar* x,
                 int size,
                 const Preconditioner<Scalar>* precond = nullptr) {
        // 共通の前処理
        PrepareIteration(A, b, x, size, precond);

        // 反復ループ（派生クラスで実装）
        return DoIterate();
    }

    // 設定
    void SetConfig(const SolverConfig& config) { config_ = config; }
    const SolverConfig& GetConfig() const { return config_; }

protected:
    // 派生クラスで実装
    virtual Result DoIterate() = 0;
    virtual void PrepareIteration(const SparseMatrixView<Scalar>& A,
                                   const Scalar* b, Scalar* x, int size,
                                   const Preconditioner<Scalar>* precond);

    SolverConfig config_;
    // 作業用ベクトル（派生クラスで使用）
    std::vector<Scalar> r_, p_, Ap_, z_;
};

} // namespace sparsesolv
```

#### 5. CG法実装

```cpp
namespace sparsesolv {

template<typename Scalar = double>
class CGSolver : public IterativeSolver<Scalar> {
public:
    std::string Name() const { return "CG"; }

protected:
    Result DoIterate() override {
        Result result;
        result.converged = false;

        // r = b - A*x
        ComputeResidual();

        // z = M^{-1} * r
        ApplyPreconditioner();

        // p = z
        p_ = z_;

        double rz_old = DotProduct(r_, z_);
        double norm_b = ComputeNorm(b_);

        for (int iter = 0; iter < config_.max_iterations; ++iter) {
            // Ap = A * p
            MatVec(A_, p_, Ap_);

            // alpha = (r, z) / (p, Ap)
            double pAp = DotProduct(p_, Ap_);
            double alpha = rz_old / pAp;

            // x = x + alpha * p
            // r = r - alpha * Ap
            UpdateSolution(alpha);

            // 収束判定
            double norm_r = ComputeNorm(r_);
            if (config_.save_residual_history) {
                result.residual_history.push_back(norm_r / norm_b);
            }

            if (norm_r / norm_b < config_.tolerance) {
                result.converged = true;
                result.iterations = iter + 1;
                result.final_residual = norm_r / norm_b;
                return result;
            }

            // z = M^{-1} * r
            ApplyPreconditioner();

            // beta = (r_new, z_new) / (r_old, z_old)
            double rz_new = DotProduct(r_, z_);
            double beta = rz_new / rz_old;
            rz_old = rz_new;

            // p = z + beta * p
            UpdateSearchDirection(beta);
        }

        result.iterations = config_.max_iterations;
        result.final_residual = ComputeNorm(r_) / norm_b;
        return result;
    }
};

} // namespace sparsesolv
```

#### 6. MRTR法実装

```cpp
namespace sparsesolv {

template<typename Scalar = double>
class MRTRSolver : public IterativeSolver<Scalar> {
public:
    std::string Name() const { return "MRTR"; }

protected:
    Result DoIterate() override {
        // MRTR特有の実装
        // CG法との違い:
        // - 残差最小化アプローチ
        // - 異なるαの計算式
        // 共通部分はIterativeSolverから継承
    }
};

} // namespace sparsesolv
```

### NGSolve統合設計

#### NGSolveアダプタ

```cpp
// ngsolve_adapter.hpp
#ifdef USE_NGSOLVE

#include <ngsolve/bla.hpp>
#include <ngsolve/la.hpp>

namespace sparsesolv {
namespace ngsolve_adapter {

// NGSolveのBaseMatrixをSparseSolvのビューに変換
template<typename SCAL>
SparseMatrixView<SCAL> WrapNGSolveMatrix(const ngla::SparseMatrix<SCAL>& mat) {
    return SparseMatrixView<SCAL>(
        mat.Height(),
        mat.Width(),
        mat.GetFirstArray().Data(),  // row pointers
        mat.GetColIndices().Data(),  // column indices
        mat.GetData().Data()         // values
    );
}

// SparseSolvのPreconditionerをNGSolveのBasePreconditionerとして使用
template<typename SCAL>
class NGSolvePreconditionerWrapper : public ngla::BasePreconditioner {
public:
    NGSolvePreconditionerWrapper(std::shared_ptr<Preconditioner<SCAL>> pre)
        : precond_(pre) {}

    void Mult(const ngla::BaseVector& x, ngla::BaseVector& y) const override {
        auto& fx = dynamic_cast<const ngla::VVector<SCAL>&>(x);
        auto& fy = dynamic_cast<ngla::VVector<SCAL>&>(y);
        precond_->Apply(fx.FV().Data(), fy.FV().Data(), fx.Size());
    }

private:
    std::shared_ptr<Preconditioner<SCAL>> precond_;
};

// NGSolveから使えるソルバー関数
template<typename SCAL>
void SolveWithSparseSolv(
    const ngla::SparseMatrix<SCAL>& A,
    const ngla::BaseVector& b,
    ngla::BaseVector& x,
    const std::string& solver_type = "ICCG",
    const SolverConfig& config = SolverConfig()
) {
    auto mat_view = WrapNGSolveMatrix(A);

    // 前処理を作成
    auto precond = CreatePreconditioner<SCAL>(solver_type, config);
    precond->Setup(mat_view);

    // ソルバーを作成
    auto solver = CreateSolver<SCAL>("CG", config);

    // 解く
    auto& fb = dynamic_cast<const ngla::VVector<SCAL>&>(b);
    auto& fx = dynamic_cast<ngla::VVector<SCAL>&>(x);

    solver->Solve(mat_view, fb.FV().Data(), fx.FV().Data(),
                  fb.Size(), precond.get());
}

} // namespace ngsolve_adapter
} // namespace sparsesolv

#endif // USE_NGSOLVE
```

### pybind11バインディング

```cpp
// python/bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "sparsesolv/core/solver_config.hpp"
#include "sparsesolv/preconditioners/ic_preconditioner.hpp"
#include "sparsesolv/core/cg_solver.hpp"

namespace py = pybind11;

PYBIND11_MODULE(sparsesolv, m) {
    m.doc() = "SparseSolv: Sparse matrix solvers for FEM";

    // SolverConfig
    py::class_<sparsesolv::SolverConfig>(m, "SolverConfig")
        .def(py::init<>())
        .def_readwrite("tolerance", &sparsesolv::SolverConfig::tolerance)
        .def_readwrite("max_iterations", &sparsesolv::SolverConfig::max_iterations)
        .def_readwrite("shift_parameter", &sparsesolv::SolverConfig::shift_parameter)
        .def_readwrite("diagonal_scaling", &sparsesolv::SolverConfig::diagonal_scaling)
        .def_readwrite("save_residual_history", &sparsesolv::SolverConfig::save_residual_history);

    // Result
    py::class_<sparsesolv::IterativeSolver<double>::Result>(m, "SolverResult")
        .def_readonly("converged", &sparsesolv::IterativeSolver<double>::Result::converged)
        .def_readonly("iterations", &sparsesolv::IterativeSolver<double>::Result::iterations)
        .def_readonly("final_residual", &sparsesolv::IterativeSolver<double>::Result::final_residual)
        .def_readonly("residual_history", &sparsesolv::IterativeSolver<double>::Result::residual_history);

    // 便利関数: NumPy配列から直接解く
    m.def("solve_iccg", [](
        py::array_t<int> row_ptr,
        py::array_t<int> col_idx,
        py::array_t<double> values,
        py::array_t<double> b,
        py::array_t<double> x,
        sparsesolv::SolverConfig config
    ) {
        // CSR形式からSparseMatrixViewを作成
        auto mat = sparsesolv::SparseMatrixView<double>(
            b.size(), b.size(),
            row_ptr.data(), col_idx.data(), values.data()
        );

        // IC前処理
        sparsesolv::ICPreconditioner<double> precond(config.shift_parameter);
        precond.Setup(mat);

        // CG法で解く
        sparsesolv::CGSolver<double> solver;
        solver.SetConfig(config);

        return solver.Solve(mat, b.data(), x.mutable_data(), b.size(), &precond);
    }, "Solve Ax=b using ICCG method",
       py::arg("row_ptr"), py::arg("col_idx"), py::arg("values"),
       py::arg("b"), py::arg("x"), py::arg("config") = sparsesolv::SolverConfig());

    // scipy.sparse.csr_matrix互換
    m.def("solve_scipy", [](
        py::object scipy_matrix,
        py::array_t<double> b,
        std::string method,
        sparsesolv::SolverConfig config
    ) {
        // scipy.sparse.csr_matrixから抽出
        py::array_t<int> indptr = scipy_matrix.attr("indptr").cast<py::array_t<int>>();
        py::array_t<int> indices = scipy_matrix.attr("indices").cast<py::array_t<int>>();
        py::array_t<double> data = scipy_matrix.attr("data").cast<py::array_t<double>>();

        int n = b.size();
        std::vector<double> x(n, 0.0);

        auto mat = sparsesolv::SparseMatrixView<double>(
            n, n, indptr.data(), indices.data(), data.data()
        );

        // 解く
        auto result = sparsesolv::Solve(mat, b.data(), x.data(), n, method, config);

        return std::make_pair(py::array_t<double>(n, x.data()), result);
    }, "Solve using scipy.sparse.csr_matrix");
}
```

### 使用例

#### C++

```cpp
#include "sparsesolv/sparsesolv.hpp"

int main() {
    // 設定
    sparsesolv::SolverConfig config;
    config.tolerance = 1e-10;
    config.shift_parameter = 1.05;

    // IC前処理
    sparsesolv::ICPreconditioner<double> precond(config.shift_parameter);
    precond.Setup(matrix);

    // CG法で解く
    sparsesolv::CGSolver<double> solver;
    solver.SetConfig(config);

    auto result = solver.Solve(matrix, b.data(), x.data(), n, &precond);

    if (result.converged) {
        std::cout << "Converged in " << result.iterations << " iterations\n";
    }
}
```

#### Python

```python
import sparsesolv
import numpy as np
from scipy.sparse import csr_matrix

# scipy.sparse行列を使用
A = csr_matrix(...)
b = np.array(...)

# 設定
config = sparsesolv.SolverConfig()
config.tolerance = 1e-10
config.shift_parameter = 1.05

# 解く
x, result = sparsesolv.solve_scipy(A, b, "ICCG", config)

print(f"Converged: {result.converged}")
print(f"Iterations: {result.iterations}")
```

#### NGSolve統合

```python
from ngsolve import *
import sparsesolv

mesh = Mesh(...)
fes = H1(mesh, order=3)
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx
a.Assemble()

f = LinearForm(fes)
f += 1*v*dx
f.Assemble()

# SparseSolvのICCGで解く
gfu = GridFunction(fes)
config = sparsesolv.SolverConfig()
config.shift_parameter = 1.05

# NGSolve行列をSparseSolvで解く
sparsesolv.solve_ngsolve(a.mat, f.vec, gfu.vec, "ICCG", config)
```

---

## 実装順序

### Phase 1: 基盤整備
1. ディレクトリ構造の変更
2. SolverConfig実装
3. SparseMatrixView実装（コピー不要のビュー）

### Phase 2: コア機能
1. Preconditioner基底クラス
2. IC分解の実装（既存コードからリファクタ）
3. IterativeSolver基底クラス
4. CGSolver実装

### Phase 3: 拡張
1. MRTRSolver実装
2. SGSPreconditioner実装
3. ABMCPreconditioner実装（マルチカラー）

### Phase 4: 統合
1. pybind11バインディング
2. NGSolveアダプタ
3. テスト・ドキュメント

---

## 変更点まとめ

| 旧 | 新 | 理由 |
|----|----|----|
| `gyo`, `retu` | `row`, `col` | 英語化 |
| `accera` | `shift_parameter` | 明確な名前 |
| `is_fix` | `isAssembled()` | 動詞形 |
| 24個のオーバーロード | 1つの`Solve()`+ビュー | シンプル化 |
| 7つのsetter | `SolverConfig`構造体 | 一元管理 |
| 独自SparseMat | SparseMatrixView | ゼロコピー |
| コード重複 | 継承+テンプレート | DRY原則 |

---

## ドキュメント構造

```
docs/
├── index.md                    # ドキュメントトップ
├── getting_started.md          # クイックスタートガイド
├── user_guide/
│   ├── installation.md         # インストール方法
│   ├── basic_usage.md          # 基本的な使い方
│   ├── solver_selection.md     # ソルバー選択ガイド
│   ├── performance_tuning.md   # 性能チューニング
│   └── ngsolve_integration.md  # NGSolve統合ガイド
├── api/
│   ├── solver_config.md        # SolverConfig API
│   ├── preconditioners.md      # 前処理クラス API
│   ├── solvers.md              # ソルバークラス API
│   └── python_bindings.md      # Python API
├── theory/
│   ├── ic_decomposition.md     # 不完全コレスキー分解の理論
│   ├── cg_method.md            # CG法の説明
│   ├── mrtr_method.md          # MRTR法の説明
│   └── shift_parameter.md      # シフトパラメータの選び方
└── examples/
    ├── cpp/
    │   ├── basic_solve.cpp
    │   ├── custom_preconditioner.cpp
    │   └── ngsolve_adapter.cpp
    └── python/
        ├── basic_solve.py
        ├── scipy_integration.py
        └── ngsolve_solve.py
```

### ドキュメント生成

- **API ドキュメント**: Doxygen（C++）+ Sphinx（Python）
- **ユーザーガイド**: MkDocs または Sphinx
- **ホスティング**: GitHub Pages

```yaml
# mkdocs.yml
site_name: SparseSolv
theme:
  name: material
nav:
  - Home: index.md
  - Getting Started: getting_started.md
  - User Guide:
    - Installation: user_guide/installation.md
    - Basic Usage: user_guide/basic_usage.md
    - NGSolve Integration: user_guide/ngsolve_integration.md
  - API Reference:
    - SolverConfig: api/solver_config.md
    - Preconditioners: api/preconditioners.md
    - Solvers: api/solvers.md
  - Theory:
    - IC Decomposition: theory/ic_decomposition.md
    - CG Method: theory/cg_method.md
```

---

## テスト構造

```
tests/
├── CMakeLists.txt              # テストビルド設定
├── cpp/
│   ├── test_main.cpp           # Google Test main
│   ├── test_sparse_matrix.cpp  # 疎行列ビューのテスト
│   ├── test_ic_decomposition.cpp
│   ├── test_cg_solver.cpp
│   ├── test_mrtr_solver.cpp
│   ├── test_preconditioners.cpp
│   ├── test_ngsolve_adapter.cpp
│   └── fixtures/
│       ├── test_matrices.hpp   # テスト用行列生成
│       └── reference_data/     # 参照解データ
├── python/
│   ├── conftest.py             # pytest設定
│   ├── test_basic.py           # 基本機能テスト
│   ├── test_scipy_compat.py    # scipy互換性テスト
│   ├── test_ngsolve.py         # NGSolve統合テスト
│   └── benchmarks/
│       ├── bench_iccg.py
│       └── bench_comparison.py # 他ソルバーとの比較
└── data/
    ├── laplacian_2d_100.mtx    # Matrix Market形式
    ├── laplacian_3d_1000.mtx
    └── fem_poisson_10k.mtx
```

### C++テスト（Google Test）

```cpp
// tests/cpp/test_cg_solver.cpp
#include <gtest/gtest.h>
#include "sparsesolv/core/cg_solver.hpp"
#include "sparsesolv/preconditioners/ic_preconditioner.hpp"
#include "fixtures/test_matrices.hpp"

class CGSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // 2Dラプラシアン行列を生成
        matrix_ = TestMatrices::Laplacian2D(50);  // 50x50グリッド
        size_ = 50 * 50;

        // 右辺ベクトル
        b_.resize(size_, 1.0);
        x_.resize(size_, 0.0);
    }

    sparsesolv::SparseMatrix<double> matrix_;
    std::vector<double> b_, x_;
    int size_;
};

TEST_F(CGSolverTest, ConvergesWithICPreconditioner) {
    sparsesolv::SolverConfig config;
    config.tolerance = 1e-10;
    config.max_iterations = 1000;

    sparsesolv::ICPreconditioner<double> precond(1.05);
    precond.Setup(matrix_);

    sparsesolv::CGSolver<double> solver;
    solver.SetConfig(config);

    auto result = solver.Solve(matrix_, b_.data(), x_.data(), size_, &precond);

    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.final_residual, config.tolerance);
    EXPECT_LT(result.iterations, 100);
}

TEST_F(CGSolverTest, DivergenceDetection) {
    // 非対称行列でCGを使うと発散する
    auto nonsym_matrix = TestMatrices::NonSymmetric(100);

    sparsesolv::SolverConfig config;
    config.detect_divergence = true;
    config.max_iterations = 100;

    sparsesolv::CGSolver<double> solver;
    solver.SetConfig(config);

    auto result = solver.Solve(nonsym_matrix, b_.data(), x_.data(), 100, nullptr);

    EXPECT_FALSE(result.converged);
}

TEST_F(CGSolverTest, SavesResidualHistory) {
    sparsesolv::SolverConfig config;
    config.save_residual_history = true;

    sparsesolv::ICPreconditioner<double> precond(1.05);
    precond.Setup(matrix_);

    sparsesolv::CGSolver<double> solver;
    solver.SetConfig(config);

    auto result = solver.Solve(matrix_, b_.data(), x_.data(), size_, &precond);

    EXPECT_GT(result.residual_history.size(), 0);
    // 残差は単調減少すべき（CGは残差を最小化）
    for (size_t i = 1; i < result.residual_history.size(); ++i) {
        EXPECT_LE(result.residual_history[i], result.residual_history[i-1] * 1.01);
    }
}
```

### Pythonテスト（pytest）

```python
# tests/python/test_basic.py
import pytest
import numpy as np
from scipy.sparse import csr_matrix, diags
import sparsesolv

@pytest.fixture
def laplacian_2d():
    """2Dラプラシアン行列を生成"""
    n = 50
    N = n * n

    # 5点ステンシル
    diagonals = [
        np.ones(N) * 4,      # 主対角
        np.ones(N-1) * -1,   # ±1
        np.ones(N-n) * -1,   # ±n
    ]
    offsets = [0, 1, n]

    A = diags(diagonals, offsets, format='csr')
    A = A + A.T - diags([A.diagonal()], [0])

    return A

@pytest.fixture
def rhs_vector(laplacian_2d):
    """右辺ベクトル"""
    return np.ones(laplacian_2d.shape[0])

class TestBasicSolve:
    def test_iccg_converges(self, laplacian_2d, rhs_vector):
        config = sparsesolv.SolverConfig()
        config.tolerance = 1e-10

        x, result = sparsesolv.solve_scipy(laplacian_2d, rhs_vector, "ICCG", config)

        assert result.converged
        assert result.final_residual < 1e-10

        # 残差を確認
        residual = np.linalg.norm(laplacian_2d @ x - rhs_vector)
        assert residual / np.linalg.norm(rhs_vector) < 1e-9

    def test_shift_parameter_effect(self, laplacian_2d, rhs_vector):
        """シフトパラメータによる収束への影響"""
        iterations = []

        for shift in [1.0, 1.05, 1.1, 1.2]:
            config = sparsesolv.SolverConfig()
            config.shift_parameter = shift
            config.tolerance = 1e-8

            x, result = sparsesolv.solve_scipy(laplacian_2d, rhs_vector, "ICCG", config)
            iterations.append(result.iterations)

        # shift=1.0より1.05の方が速いはず（適度なシフトで安定化）
        assert iterations[1] <= iterations[0]

    def test_residual_history(self, laplacian_2d, rhs_vector):
        config = sparsesolv.SolverConfig()
        config.save_residual_history = True

        x, result = sparsesolv.solve_scipy(laplacian_2d, rhs_vector, "ICCG", config)

        assert len(result.residual_history) == result.iterations
        assert result.residual_history[-1] < result.residual_history[0]

class TestSciPyCompatibility:
    def test_csr_matrix_input(self, laplacian_2d, rhs_vector):
        """scipy.sparse.csr_matrix入力"""
        x, result = sparsesolv.solve_scipy(laplacian_2d, rhs_vector, "ICCG")
        assert result.converged

    def test_csc_matrix_input(self, laplacian_2d, rhs_vector):
        """scipy.sparse.csc_matrixは自動変換"""
        csc_matrix = laplacian_2d.tocsc()
        x, result = sparsesolv.solve_scipy(csc_matrix, rhs_vector, "ICCG")
        assert result.converged
```

### ベンチマーク

```python
# tests/python/benchmarks/bench_comparison.py
import pytest
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg, spilu, LinearOperator
import sparsesolv
import time

def create_3d_laplacian(n):
    """3Dラプラシアン行列（n^3 DOF）"""
    # 実装省略
    pass

@pytest.mark.benchmark
class TestPerformance:
    @pytest.fixture(params=[10, 20, 30, 40])
    def problem_size(self, request):
        n = request.param
        A = create_3d_laplacian(n)
        b = np.ones(n**3)
        return A, b, n**3

    def test_sparsesolv_iccg(self, problem_size, benchmark):
        A, b, ndof = problem_size

        def solve():
            config = sparsesolv.SolverConfig()
            config.tolerance = 1e-8
            return sparsesolv.solve_scipy(A, b, "ICCG", config)

        result = benchmark(solve)
        assert result[1].converged

    def test_scipy_cg_ilu(self, problem_size, benchmark):
        A, b, ndof = problem_size

        def solve():
            ilu = spilu(A.tocsc())
            M = LinearOperator(A.shape, ilu.solve)
            x, info = cg(A, b, M=M, tol=1e-8)
            return x, info

        result = benchmark(solve)
        assert result[1] == 0  # scipy: 0 = success
```

### CI/CD設定（GitHub Actions）

```yaml
# .github/workflows/test.yml
name: Test

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  test-cpp:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest]
        build_type: [Release, Debug]

    steps:
    - uses: actions/checkout@v4

    - name: Configure CMake
      run: cmake -B build -DCMAKE_BUILD_TYPE=${{ matrix.build_type }} -DBUILD_TESTING=ON

    - name: Build
      run: cmake --build build --config ${{ matrix.build_type }}

    - name: Test
      run: ctest --test-dir build -C ${{ matrix.build_type }} --output-on-failure

  test-python:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v4

    - name: Set up Python
      uses: actions/setup-python@v5
      with:
        python-version: '3.12'

    - name: Install dependencies
      run: |
        pip install numpy scipy pytest pytest-benchmark
        pip install ./  # install sparsesolv

    - name: Run tests
      run: pytest tests/python -v

    - name: Run benchmarks
      run: pytest tests/python/benchmarks --benchmark-only

  test-ngsolve:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v4

    - name: Install NGSolve
      run: pip install ngsolve

    - name: Install sparsesolv
      run: pip install ./

    - name: Test NGSolve integration
      run: pytest tests/python/test_ngsolve.py -v
```

---

## テスト方針

### 単体テスト
- 各クラス・関数に対するテスト
- 境界条件、エッジケースの確認
- カバレッジ目標: 80%以上

### 統合テスト
- ソルバー全体の動作確認
- NGSolveとの統合テスト
- scipy互換性テスト

### 回帰テスト
- 既知の問題に対する解の検証
- 参照解（直接法）との比較
- 収束履歴の検証

### 性能テスト
- 問題サイズに対するスケーリング
- 他ソルバー（scipy, Eigen）との比較
- メモリ使用量の計測
