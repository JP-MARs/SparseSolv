# SparseSolv
SparseMatrix and Matrix Solvers including the shifted-ICCG

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
This method uses A-formula withough coulomb gauge and regularization terms.
Original A-fomula is solved by shifted-ICCG.

Voxel Mesh data (original format) can be obtained from the following links:
https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_850000ele.zip
https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_2000000ele.zip
Please unzip the zip-file, and put "*.vxldata" files at "examples/VoxelFEM/MeshData/".
The  "Defines_Mesh.h" at "examples/VoxelFEM/" is overwrite by the same file in the zip-file.

