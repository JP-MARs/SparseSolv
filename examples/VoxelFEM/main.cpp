
#include "VoxelMesh.h"
#include <fstream>
#include <SparseSolv/MatSolverICCG.hpp>
#include <SparseSolv/MatSolverMRTR.hpp>

#ifdef EIGEN_NO_DEBUG
#pragma message("EIGEN_NO_DEBUG=ON")
#else
#pragma message("EIGEN_NO_DEBUG=OFF")
#endif

#ifdef EIGEN_DONT_VECTORIZE
#pragma message("EIGEN_DONT_VECTORIZE=ON")
#else
#pragma message("EIGEN_DONT_VECTORIZE=OFF")
#endif

#ifdef NDEBUG
#pragma message("NDEBUG=ON")
#else
#pragma message("NDEBUG=OFF")
#endif

#include <chrono>
using clock_type = std::chrono::steady_clock;

using namespace std;

/*//=======================================================
  // ●　均質化ドライバ
  //=======================================================*/
void VoxelFEM_Driver(VoxelMesh *mesh){

	/* 材料設定 */
	cout << "set material" << endl;
	mesh->setMaterials();

	/* 境界条件設定 */
	cout << "set boundary" << endl;
	mesh->setBoundary();

	/* 電流設定T法 */
	cout << "set T vector" << endl;
	mesh->set_Tvec();

	/* 辺のサブIDセット */
	cout << "set edge id in solver" << endl;
	mesh->setSubIDs();

	/* 書き出し */
	//mesh->writeVTK(nullptr);

	/* サイズ計算 */
	const int edge_size = mesh->getEdgeNum();
	const int A0_size = mesh->getA0Num();
	const int all_size = edge_size - A0_size;
	cout << "total edge num : " << edge_size << endl;
	cout << "edge on boundary : " << A0_size << endl;
	cout << "unknwon num : " << all_size << endl;
	
	/* 行列系統作成 */	
	double *RighthandVec = new double[all_size];	
	/* 右辺作成 */
	cout << "make righthand vector " << endl;
	mesh->make_right_vector(RighthandVec);
	/* マトリックス作成 */
	cout << "make matrix" << endl;
	auto t0 = clock_type::now();
	JPMRspace::SparseSolv::SparseMat matA;
	mesh->make_matrix(matA);
	auto t1 = clock_type::now();
	double t_mat = std::chrono::duration<double, std::micro>(t1 - t0).count();
	cout << "mat time = " << t_mat / 1000.0 << endl;
	matA.printMat("mat.csv");
	
	/* ICCGスタート */
	double *result = new double[edge_size];
	for(int i = 0 ; i < edge_size ; i++){
		result[i] = 0;
	}
	JPMRspace::SparseSolv::MatSolverICCG solver;
	cout << "solve start " << endl;
	auto tA = clock_type::now();
	//omp_set_num_threads(1);
	solver.setSaveLog(true);
	solver.setAccelFactor(GA_IC);	
	solver.solve(EPS_SOLVER, MAX_ITR_ICCG, matA, RighthandVec, result);
	auto tB = clock_type::now();
	double t_slv = std::chrono::duration<double, std::micro>(tB - tA).count();
	cout << "mat t_slv = " << t_slv / 1000.0 << endl;
	cout << "solver end" << endl;
	vector<double> log;
	solver.getResidualLog(log);
	fstream fp("residual_log.csv", std::ios::out);
	for(int i = 0 ; i < log.size() ; i++){
		fp << i << ", " << log[i] << endl;
		//cout << i << ", " << log[i] << endl;
	}
	fp.close();


	/* 書き出し */
	cout << "output vtk file" << endl;
	mesh->writeVTK(result);
	
	/* インダクタンス計算 */
	/*double inductor = 0;
	for(int i = 0 ; i < all_size ; i++){
		inductor += RighthandVec[i] * result[i];
	}
	inductor *= 8.0;
	inductor /= (CURRENT*CURRENT);*/

	delete[] RighthandVec;
	delete[] result;
}
/*//=======================================================
  // ●　メイン関数
  //=======================================================*/
int main(int argc, char *argv[]){

	double length1[] = {X_L1, Y_L1, Z_L1};
	double length2[] = {X_L2, Y_L2, Z_L2};
	int change[] = {MESH_X1, MESH_Y1, MESH_Z1};
	//cout << X_L1 <<", " << X_L2 << ", " << MESH_X1 << endl;
	/* メッシュデータ生成 */
	cout << "read files " << endl;
	VoxelMesh mesh(length1, length2, change);

	VoxelFEM_Driver(&mesh);
	cout << "end!"<<endl;
	return 1;
}
