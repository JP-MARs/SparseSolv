project(SparseSolv CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${Basedir})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${Basedir})

#コンパイル設定
add_library(SparseSolv
  STATIC
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
 message("Intel Classical Compiler use in SparseSolv !") 
 target_include_directories(SparseSolv PUBLIC ${Basedir})
 target_compile_options(SparseSolv PUBLIC -std=c++14 -O3 -qopenmp -ipo)
 target_link_libraries(SparseSolv PUBLIC -qmkl)
#
# set for Intel OneAPI
elseif(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolv !") 
 set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 target_include_directories(SparseSolv PUBLIC ${Basedir})
 target_compile_options(SparseSolv PUBLIC -std=c++14 -O3 -qopenmp)
 target_link_libraries(SparseSolv PUBLIC -qmkl -ipo)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolv PUBLIC ${Basedir})
  target_compile_options(SparseSolv PUBLIC -std=c++14 -O3 -fopenmp)
  target_link_libraries(SparseSolv PUBLIC OpenMP::OpenMP_CXX)
endif()
