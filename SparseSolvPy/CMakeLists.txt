project(SparseSolvPy CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message("Pybind export in SparseSolvPy !") 
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    set(PYBIND11_FINDPYTHON ON)
endif()
find_package(pybind11 REQUIRED)
pybind11_add_module(SparseSolvPy 
  MODULE 
    pybind_SparseMatPy.cpp
    pybind_MatSolversPy.cpp
    SparseMat.cpp
    SparseMatC.cpp
    SparseMatOperators.cpp
    MatSolvers_Base.cpp
    MatSolvers_ICCG.cpp
    MatSolvers_ICCOCG.cpp
    MatSolvers_ICMRTR.cpp
    MatSolvers_SGSMRTR.cpp
    MatSolvers_Eigen.cpp
    MatSolvers_MKL.cpp
    MatSolvers_ABMCICCG.cpp
)
# setting for old intel
if(CMAKE_CXX_COMPILER MATCHES "icpc")
 message("Intel Classical Compiler use in SparseSolvPy !") 
 target_include_directories(SparseSolvPy PUBLIC ${SparseSolvBase_SOURCE_DIR})
 target_compile_options(SparseSolvPy PUBLIC -O3 -qopenmp -ipo)
 target_link_libraries(SparseSolvPy PUBLIC -qmkl  SparseSolv)
#
# set for Intel OneAPI
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolvPy !") 
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 set(CMAKE_RANLIB "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ranlib")
 target_include_directories(SparseSolvPy PUBLIC ${SparseSolvBase_SOURCE_DIR})
 target_compile_options(SparseSolvPy PUBLIC -O3 -qopenmp -fPIC)
 target_link_libraries(SparseSolvPy PUBLIC  SparseSolv -qmkl -ipo -shared)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolvPy PUBLIC ${SparseSolvBase_SOURCE_DIR})
  target_compile_options(SparseSolvPy PUBLIC -O3 -fPIC)
  target_link_libraries(SparseSolvPy PUBLIC  SparseSolv OpenMP::OpenMP_CXX)
endif()
