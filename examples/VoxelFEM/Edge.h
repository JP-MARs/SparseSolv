#ifndef DEF_EDGE_VOXEL
#define DEF_EDGE_VOXEL

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;
/*
//=======================================================
// ■ Edge
//=======================================================
// FEMの辺定義クラス
//=======================================================*/
class Edge {
private:
	int ID;												/* この辺のID */
	int subID;											/* 行列計算用のサブID */
	int stat_node_id;									/* 辺の始点NodeのID */
	int end_node_id;									/* 辺の終点NodeのID */
	bool zero_edge;										/* 既知ゼロの辺か */
	bool T_edge;										/* 電流ベクトルポテンシャルの辺か */
	double T_potential;									/* 電流ベクトルポテンシャルの値 */
public:
	Edge();												/* コンストラクタ */
	Edge(int sID, int eID);								/* コンストラクタ2 */
	void set(int sID, int eID);							/* 座標セット */
	void setSubID(int i){subID = i;};
	int getSubID(){return subID;};
	int getSID(){return stat_node_id;};
	int getEID(){return end_node_id;};
	int getID(){return ID;};
	void setID(int id){ID=id;};
	void setZeroEdge(bool bl){zero_edge=bl;};
	bool isZeroEdge(){return zero_edge;};
	void setT_Edge(bool bl){T_edge=bl;};
	bool isT_Edge(){return T_edge;};
	double get_Tvalue(){return T_potential;};
	void set_Tvalue(double x){T_potential = x;};
};

#endif

