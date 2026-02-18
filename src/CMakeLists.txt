project(SparseSolv CXX)
#出力先はベース
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SparseSolvBase_BINARY_DIR})

#message("BaseLIBOUT = ${SparseSolvBase_BINARY_DIR}" )
#message("DIR = ${CMAKE_SOURCE_DIR}" )

file(GLOB SRC_FILES *.cpp)

#コンパイル設定
add_library(SparseSolv
  STATIC ${SRC_FILES}
)

# setting for old intel
# set for Intel OneAPI
if(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in SparseSolv !") 
 #set(CMAKE_AR "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-ar")
 #set(CMAKE_NM "/opt/intel/oneapi/compiler/latest/bin/compiler/llvm-nm")
 target_include_directories(SparseSolv 
   PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}
 )
 target_compile_options(SparseSolv PUBLIC -std=c++17 -O3 -qopenmp)
 target_link_libraries(SparseSolv PUBLIC -qmkl -ipo)
#
# gcc setting
else()
  find_package(OpenMP REQUIRED)
  target_include_directories(SparseSolv 
    PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}
  )
  target_compile_options(SparseSolv PUBLIC -std=c++17 -O3 -fopenmp)
  target_link_libraries(SparseSolv PUBLIC OpenMP::OpenMP_CXX)
endif()
