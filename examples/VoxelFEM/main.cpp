
#include "VoxelMesh.h"
#include <fstream>
#include <MatSolversICCG.hpp>

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
	SRLfem::SparseMat matA(all_size);
	mesh->make_matrix(matA);
	matA.fix();

	/* ICCGスタート */
	double *result = new double[edge_size];
	for(int i = 0 ; i < edge_size ; i++){
		result[i] = 0;
	}
	SRLfem::MatSolversICCG iccg;
	cout << "solve start " << endl;
	//omp_set_num_threads(1);
	iccg.setSaveLog(true);
	iccg.solveICCG(all_size, EPS_SOLVER, MAX_ITR_ICCG, GA_IC, matA, RighthandVec, result);
	cout << "solver end" << endl;
	vector<double> log;
	iccg.getResidualLog(log);
	fstream fp("residual_log.csv", std::ios::out);
	for(int i = 0 ; i < log.size() ; i++){
		fp << i << ", " << log[i] << endl;
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
