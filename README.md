# SparseSolv
SparseMatrix and Matrix Solvers including 
- shifted-ICCG
- multicolor ordering shifted ICCG
- shifted-IC+MRTR
- Eisenstat's Symmetric Gauss-Seidel-MRTR.

## SparseSolv
Provide sparse matrix and its linear solver.

## SparseSolvPy
Python binding of the SparseSolv Using Pybind11

## examples
example:
### Pybind_example.py
A simple sample using python
### Voxel FEM
A simple C++ sample based on voxel elements.
This method uses A-formula withough coulomb gauge and regularization terms.<br>
Original A-fomula is solved by shifted-ICCG.

Voxel Mesh data (original format) can be obtained from the following links:<br>
https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_850000ele.zip<br>
https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_2000000ele.zip<br>
Please unzip the zip-file, and put "*.vxldata" files at "examples/VoxelFEM/MeshData/".<br>
The  "Defines_Mesh.h" at "examples/VoxelFEM/" is overwrite by the same file in the zip-file.<br>

# Contributor
 - Takahiro Sato (Muroran institute of technology, JAPAN)
 - Shingo Hiruma (Kyoto University, JAPAN)
 - Kengo Sugahara (Kindai University, JAPAN)
 - Tomonori Tsuburaya (Fukuoka University, JAPAN)

# 本ライブラリの説明 (explain in Japanese)
 本ライブラリは、日本の磁界系数値解析の研究者による疎行列ソルバのライブラリです。電磁界有限要素法のソルバとして広く使われている。加速係数付きICCG法の線形ソルバです。また、MRTR法も実装されています。<br>
 本ライブラリはpythonから使えるように、pybind11で一部機能をpython用ライブラリとして公開しており、疎行列クラスとソルバをpythonから使うことができます。<br>
 


