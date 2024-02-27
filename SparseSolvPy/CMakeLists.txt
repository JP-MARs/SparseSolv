project(SparseSolvPy CXX)
# set base dir
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Basedir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Basedir})

message("Pybind export in SparseSolv !") 
find_package(pybind11 REQUIRED)
pybind11_add_module(SparseSolvPy 
  MODULE 
    SparseMatPy.cpp
    SparseMat.cpp
    SparseMatC.cpp
    SparseMatOperators.cpp
    MatSolvers_ICCG.cpp
    MatSolvers_ICCOCG.cpp
    MatSolvers_Eigen.cpp
    MatSolvers_MKL.cpp
    MatSolvers_ABMCICCG.cpp
)

# setting for old intel
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in SparseSolvPy !") 
 target_include_directories(SparseSolvPy PUBLIC ${Basedir})
 target_compile_options(SparseSolvPy PUBLIC -std=c++14 -O3 -qopenmp -ipo)
 target_link_libraries(SparseSolvPy PUBLIC -qmkl  SparseSolv)
#
# set for Intel OneAPI
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolvPy !") 
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 set(CMAKE_RANLIB "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ranlib")
 target_include_directories(SparseSolvPy PUBLIC ${Basedir})
 target_compile_options(SparseSolvPy PUBLIC -std=c++14 -O3 -qopenmp -fPIC)
 target_link_libraries(SparseSolvPy PUBLIC  SparseSolv -qmkl -ipo -shared)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolvPy PUBLIC ${Basedir})
  target_compile_options(SparseSolvPy PUBLIC -std=c++14 -O3 -fopenmp -fPIC)
  target_link_libraries(SparseSolvPy PUBLIC  SparseSolv OpenMP::OpenMP_CXX)
endif()


