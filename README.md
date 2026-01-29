# SparseSolv

SparseMatrix and Matrix Solvers including:
- shifted-ICCG
- multicolor ordering shifted ICCG
- shifted-IC+MRTR
- Eisenstat's Symmetric Gauss-Seidel-MRTR

## Origin

This library is based on [JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv), with additional features for NGSolve integration.

## Features

### Header-only Library (v2.0)

The new version provides a modern C++17 header-only library with:
- Zero-copy `SparseMatrixView` for CSR matrices
- Template-based preconditioners (IC, SGS)
- Iterative solvers (CG, MRTR)
- scipy.sparse compatibility via Python bindings

### NGSolve Integration

SparseSolv can be used as a preconditioner with NGSolve's Krylov solvers:

```python
from ngsolve import *
from ngsolve.krylovspace import CGSolver
import sparsesolv

# Create bilinear form and assemble
a = BilinearForm(fes)
a += grad(u)*grad(v)*dx
a.Assemble()

# Create IC preconditioner
pre = sparsesolv.ICPreconditioner(a.mat, shift=1.05)
pre.Update()

# Use with NGSolve's CGSolver
inv = CGSolver(a.mat, pre, printrates=True, tol=1e-10)
gfu.vec.data = inv * f.vec
```

Supported vector types:
- `VVector<double>`: Standard dense vectors
- `BlockVector`: Block vectors (flattened access)

### SetGeomInfo API Support

This fork is designed to work with [ksugahar/ngsolve](https://github.com/ksugahar/ngsolve) which includes the SetGeomInfo API (PR #232) for high-order curving of externally imported meshes.

## Building

### Requirements

- CMake 3.16+
- C++17 compiler
- Python 3.x (for Python bindings)
- pybind11 (for Python bindings)
- NGSolve (optional, for NGSolve integration)

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `SPARSESOLV_BUILD_TESTS` | ON | Build tests |
| `SPARSESOLV_BUILD_PYTHON` | ON | Build Python bindings |
| `SPARSESOLV_USE_OPENMP` | ON | Use OpenMP for parallelization |
| `SPARSESOLV_USE_NGSOLVE` | ON | Build with NGSolve integration |
| `SPARSESOLV_BUILD_LEGACY` | OFF | Build legacy SparseSolv library |

### Build Commands

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

With NGSolve:
```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DNGSOLVE_DIR=/path/to/ngsolve
cmake --build .
```

## Python API

### scipy.sparse Support

```python
import sparsesolv
import scipy.sparse as sp
import numpy as np

# Create sparse matrix
A = sp.csr_matrix(...)
b = np.array(...)

# Solve
config = sparsesolv.SolverConfig()
config.tolerance = 1e-10
config.max_iterations = 1000

result = sparsesolv.solve_iccg(A, b, config)
print(f"Converged: {result.converged}, iterations: {result.iterations}")
```

### Available Solvers

- `solve_iccg(A, b, config)` - Incomplete Cholesky Conjugate Gradient
- `solve_icmrtr(A, b, config)` - IC + MRTR method
- `solve_sgsmrtr(A, b, config)` - Symmetric Gauss-Seidel + MRTR method

### NGSolve Preconditioners

- `sparsesolv.ICPreconditioner(mat, shift=1.05)` - Incomplete Cholesky
- `sparsesolv.SGSPreconditioner(mat)` - Symmetric Gauss-Seidel

## Legacy Library

The original SparseSolv library is still available but deprecated. Enable with:
```bash
cmake .. -DSPARSESOLV_BUILD_LEGACY=ON
```

## Examples

### Pybind_example.py
A simple sample using Python with the legacy interface.

### Voxel FEM
A simple C++ sample based on voxel elements.
This method uses A-formula without coulomb gauge and regularization terms.
Original A-formula is solved by shifted-ICCG.

Voxel Mesh data (original format) can be obtained from the following links:
- https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_850000ele.zip
- https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_2000000ele.zip

Please unzip the zip-file, and put `*.vxldata` files at `examples/VoxelFEM/MeshData/`.
The `Defines_Mesh.h` at `examples/VoxelFEM/` is overwritten by the same file in the zip-file.

## Contributors

Original SparseSolv:
- Takahiro Sato (Muroran Institute of Technology, JAPAN)
- Shingo Hiruma (Kyoto University, JAPAN)
- Kengo Sugahara (Kindai University, JAPAN)
- Tomonori Tsuburaya (Fukuoka University, JAPAN)

NGSolve integration:
- Kengo Sugahara (Kindai University, JAPAN)

## 本ライブラリの説明 (Japanese)

本ライブラリは、日本の磁界系数値解析の研究者による疎行列ソルバのライブラリです。電磁界有限要素法のソルバとして広く使われている加速係数付きICCG法の線形ソルバです。また、MRTR法も実装されています。

本ライブラリはpythonから使えるように、pybind11で一部機能をpython用ライブラリとして公開しており、疎行列クラスとソルバをpythonから使うことができます。

### v2.0での変更点

- ヘッダオンリーライブラリとして再設計
- NGSolveとの統合（前処理としてNGSolveのCGSolverで使用可能）
- scipy.sparseとの互換性
- C++17対応
- メモリ管理の改善（RAII、Rule of Five）

## License

See LICENSE file for details.
