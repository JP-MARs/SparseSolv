project(SparseSolvPy CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})

#message("BaseLIBOUTPY = ${SparseSolvBase_BINARY_DIR}" )
#message("DIRPY = ${PROJECT_SOURCE_DIR}" )

file(GLOB SRC_FILES *.cpp)

#コンパイル設定
message("Pybind export in SparseSolvPy !") 
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    set(PYBIND11_FINDPYTHON ON)
endif()
find_package(pybind11 REQUIRED)
pybind11_add_module(SparseSolvPy 
  MODULE  ${SRC_FILES}
)

# setting for old intel
# set for Intel OneAPI
if(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolvPy !") 
 #set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 #set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 target_include_directories(SparseSolvPy 
   PUBLIC ${PROJECT_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}
 )
 target_compile_options(SparseSolvPy PUBLIC -std=c++17 -O3 -qopenmp)
 target_link_libraries(SparseSolvPy PUBLIC -qmkl -ipo)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolvPy 
    PUBLIC ${PROJECT_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}
  )
  target_compile_options(SparseSolvPy PUBLIC -std=c++17 -O3 -fopenmp)
  target_link_libraries(SparseSolvPy PUBLIC OpenMP::OpenMP_CXX)
endif()
