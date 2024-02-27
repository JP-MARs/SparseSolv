#ifndef DEF_ELE_VOXEL
#define DEF_ELE_VOXEL

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;

#include "Node.h"
#include "Edge.h"


/*
//=======================================================
// ■ Element
//=======================================================
// FEMの要素定義クラス
//=======================================================*/
class Element {
private:
	int ID;																	/* この要素のID */
	int x;																	/* 原点から何番目の要素か(ｘ方向) */
	int y;																	/* 原点から何番目の要素か(ｙ方向) */
	int z;																	/* 原点から何番目の要素か(ｚ方向) */
	int *nodes;																/* この要素の節点たち */
	int *edges;																/* この要素の辺たち */
	int mat_type;															/* 要素の材質情報 */
	double permeability;													/* 透磁率 */
	static Node *NodeData;													/* 節点データ保持配列 */
	static Edge *EdgeData;													/* 辺データ保持配列 */
	void calcRot(int k, double r, double s, double t, double *rot);
public:
	Element();																/* コンストラクタ */
	~Element();
	void set(int xx, int yy, int zz, int *n, int *e, int id);				/* セット */
	int getID(){return ID;};
	int *getNodes(){return nodes;};
	int *getEdges(){return edges;};
	int getNodeID(int i){return nodes[i];};
	int getEdgeID(int i){return edges[i];};
	int getX(){return x;};
	int getY(){return y;};
	int getZ(){return z;};
	void getGrav(double *x);												/* 要素の中心位置返し */
	double getMinX(){return NodeData[nodes[0]].getX();};					/* 要素の最小X位置返し */
	double getMinY(){return NodeData[nodes[0]].getY();};					/* 要素の最小Y位置返し */
	double getMinZ(){return NodeData[nodes[0]].getZ();};					/* 要素の最小Z位置返し */
	double getMaxX(){return NodeData[nodes[6]].getX();};					/* 要素の最大X位置返し */
	double getMaxY(){return NodeData[nodes[6]].getY();};					/* 要素の最大Y位置返し */
	double getMaxZ(){return NodeData[nodes[6]].getZ();};					/* 要素の最大Z位置返し */	
	static void setNodeData(Node *n){NodeData = n;};
	static void setEdgeData(Edge *n){EdgeData = n;};
	void setMaterial(int x, double perm){mat_type = x;permeability=perm;};	/* 要素の材質指定 */
	int getMatType(){return mat_type;};										/* 要素の材質情報取得 */
	double getMatPermeability(){return permeability;};						/* 要素の透磁率返し */

	void calcB(double *B, double *A);										/* 要素中心の磁束密度計算 */
	void calcRot(int k, double *pos_rst, double *rst, double *rot);
	void calcB(double *B, double *A, double *pos);
	void calcJ(double *J, double *pos);
	bool isInclude(double *pos);
};

#endif

