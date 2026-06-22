project(3DFEM CXX)
#出力先はベース
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/VoxelFEM)

file(GLOB SRC_FILES *.cpp)
#コンパイル設定
add_executable(
 voxel3DFEM.out ${SRC_FILES}
)

#オプションとライブラリリンクの設定
if(CMAKE_CXX_COMPILER MATCHES "icpx")
 message("Intel New Compiler use in Example !") 
 target_include_directories(SparseSolv 
   PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}/000_thirdparty  ${CMAKE_SOURCE_DIR}
 )
 target_compile_options(voxel3DFEM.out PUBLIC -std=c++14 -O3 -xhost -qopenmp)
 target_link_libraries(voxel3DFEM.out PUBLIC -qmkl -ipo SparseSolv)
#
else()
 target_include_directories(SparseSolv 
   PUBLIC ${CMAKE_SOURCE_DIR}/include  ${CMAKE_SOURCE_DIR}/000_thirdparty  ${CMAKE_SOURCE_DIR}
 )
 target_compile_options(voxel3DFEM.out PUBLIC -std=c++14 -O3 -fopenmp)
 target_link_libraries(voxel3DFEM.out PUBLIC -lgomp SparseSolv)
endif()
