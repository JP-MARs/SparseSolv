#ifndef DEF_VOXEL_MAIN_VOXEL
#define DEF_VOXEL_MAIN_VOXEL

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

#include <vector>

#include "Defines_Mesh.h"
#include "Defines_Other.h"
#include "Node.h"
#include "Edge.h"
#include "Element.h"

#include <SparseMat.hpp>

/*
//=======================================================
// ■ VoxelMesh
//=======================================================
// FEMのボクセルメッシュ自動生成プログラム
//=======================================================*/
class VoxelMesh {
private:
	int x_num;												/* ｘ方向の要素数 */
	int y_num;												/* ｙ方向の要素数 */
	int z_num;												/* ｚ方向の要素数 */
	int node_num;											/* 節点の総数 */
	int edge_num;											/* 辺の総数 */
	int element_num;										/* 要素の総数 */
	double x_length1;										/* 各要素のｘ長さ */
	double y_length1;										/* 各要素のｙ長さ */
	double z_length1;										/* 各要素のｚ長さ */
	double x_length2;										/* 各要素のｘ長さ、その２ */
	double y_length2;										/* 各要素のｙ長さ、その２ */
	double z_length2;										/* 各要素のｚ長さ、その２ */
	int change_x;											/* 要素のサイズ変更になる位置Ｘ */
	int change_y;											/* 要素のサイズ変更になる位置Ｙ */
	int change_z;											/* 要素のサイズ変更になる位置Ｚ */
	Node *NodeData;											/* 節点データ保持配列 */
	Edge *EdgeData;											/* 辺データ保持配列 */
	Element *elements;										/* 要素データ保持配列 */
	int A0_size;											/* ゼロ既知になる辺のサイズ */
	int T_size;												/* 電流ポテンシャルを持つ辺のサイズ */
	/* サブルーチン群 */
	/* ------------------------------------------------------------------------------------- */
	void setNodeData();										/* Nodeデータ読み取り  */
	void setEdgeData();										/* Edgeデータ読み取り */
	void setElementData();									/* Elementデータ読み取り */
	void setMatCoefficients(int k, double& sk, double& rk, double& tk);				/* 行列作成に使う係数設定ルーチン */
	void calc_v_tensor(double *v_tensor, Element *Ele_ptr);							/* 磁気抵抗率ゲット */
	void setTValue(Edge *tempE, double xg, double yg, double zg, double dz);

public:
	VoxelMesh(double *xyz_L1, double *xyz_L2, int *change_xyz);				/* コンストラクタ */
	~VoxelMesh();
	Edge *getEdges(){return EdgeData;};
	Element *getElement(){return elements;};
	int getNodeNum(){return node_num;};
	int getEdgeNum(){return edge_num;};
	int getElementNum(){return element_num;};
	int getA0Num(){return A0_size;};
	void setSubIDs();
	void make_right_vector(double *vecB);									/* 右辺ベクトル作成 */
	void make_matrix(SRLfem::SparseMat& matA);								/* 全体係数マトリックス作成 */
	void setMaterials();													/* 各要素の材料情報作成 */
	void setBoundary();
	void set_Tvec();
	/**/
	void writeVTK(double* resultA);
};

#endif
