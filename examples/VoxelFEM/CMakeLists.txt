project(3DFEM CXX)
#出力先はベース
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Basedir}/examples/VoxelFEM)

#コンパイル設定
add_executable(
 voxel3DFEM.out 
   main.cpp
   Node.cpp
   Edge.cpp
   Element.cpp
   setBoundary.cpp
   setMaterial.cpp
   setT_Vector.cpp
   VoxelMatrix.cpp
   VoxelMesh.cpp
)
#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in Example !") 
 target_compile_options(voxel3DFEM.out PUBLIC -std=c++14 -O3 -I${Basedir} -I${Basedir}/SparseSolv -qopenmp -ipo)
 target_link_libraries(voxel3DFEM.out PUBLIC -qmkl -ipo SparseSolv)
#
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in Example !") 
 target_compile_options(voxel3DFEM.out PUBLIC -std=c++14 -O3 -I${Basedir} -I${Basedir}/SparseSolv -qopenmp)
 target_link_libraries(voxel3DFEM.out PUBLIC -qmkl -ipo SparseSolv)
#
else()
 find_package(OpenMP REQUIRED)
 target_compile_options(voxel3DFEM.out PUBLIC -std=c++14 -O3 -I${Basedir} -I${Basedir}/SparseSolv)
 target_link_libraries(voxel3DFEM.out PUBLIC OpenMP::OpenMP_CXX SparseSolv)
endif()
