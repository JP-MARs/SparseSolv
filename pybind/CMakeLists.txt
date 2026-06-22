project(SparseSolvPy CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})

message("BaseLIBOUTPY = ${SparseSolvBase_BINARY_DIR}" )
message("DIRPY = ${PROJECT_SOURCE_DIR}" )

file(GLOB SRC_FILES *.cpp)

#コンパイル設定
message("Pybind export in SparseSolvPy !") 
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    set(PYBIND11_FINDPYTHON ON)
endif()

# ------------------------------------------------------------
# 1. 仮想環境が activate されている場合、その Python を優先する
# ------------------------------------------------------------
if(DEFINED ENV{VIRTUAL_ENV})
    set(Python_EXECUTABLE "$ENV{VIRTUAL_ENV}/bin/python")
endif()
# ------------------------------------------------------------
# 2. Python を検出
#    pybind11 の Python 拡張モジュールを作るので Development.Module も要求
# ------------------------------------------------------------
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

message(STATUS "Python_EXECUTABLE = ${Python_EXECUTABLE}")

# ------------------------------------------------------------
# 3. 現在の Python 環境に入っている pybind11 の CMake 設定ディレクトリを取得
# ------------------------------------------------------------
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE pybind11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE pybind11_CMAKE_DIR_RESULT
)
message(STATUS "pybind11_CMAKE_DIR = ${pybind11_CMAKE_DIR}")

# ------------------------------------------------------------
# 4. CMake に pybind11Config.cmake の場所を教える
# ------------------------------------------------------------
set(pybind11_DIR "${pybind11_CMAKE_DIR}" CACHE PATH "pybind11 cmake directory" FORCE)

set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)
#find_package(pybind11 REQUIRED)



pybind11_add_module(SparseSolvPy 
  MODULE  ${SRC_FILES}
)

# set for Intel OneAPI
if(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolvPy !") 
 #set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 #set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 target_include_directories(SparseSolvPy 
   PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}/000_thirdparty  ${CMAKE_SOURCE_DIR}
 )
 target_compile_options(SparseSolvPy PUBLIC -std=c++17 -O3 -qopenmp)
 target_link_libraries(SparseSolvPy PUBLIC -qmkl -ipo)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolvPy 
     PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}/000_thirdparty  ${CMAKE_SOURCE_DIR}
  )
  target_compile_options(SparseSolvPy PUBLIC -std=c++17 -O3 -fopenmp)
  target_link_libraries(SparseSolvPy PUBLIC OpenMP::OpenMP_CXX)
endif()
