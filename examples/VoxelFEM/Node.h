#ifndef DEF_NODE_VOXEL
#define DEF_NODE_VOXEL

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;
/*
//=======================================================
// ■ Node
//=======================================================
// FEMの節点定義クラス
//=======================================================*/
class Node {
private:
	int ID;													/* このノードのID */
	double x;												/* X座標 */
	double y;												/* Y座標 */
	double z;												/* Z座標 */
public:
	Node(){;};												/* コンストラクタ */
	Node(int I, double xx, double yy, double zz);			/* コンストラクタ2 */
	void set(int I, double xx, double yy, double zz);		/* 座標セット */
	double getX(){return x;};
	double getY(){return y;};
	double getZ(){return z;};
	int getID(){return ID;};
	void setID(int id){ID=id;};
};



#endif