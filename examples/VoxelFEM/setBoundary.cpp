#include "VoxelMesh.h"

/*//=======================================================
  // ● 各要素の材料情報作成
  //=======================================================*/
void VoxelMesh::setBoundary(){
	A0_size = 0;
	for(int i = 0 ; i < edge_num ; i++){
		Edge *tempE = EdgeData+i;
		/* 始点と終点ゲット */
		int sID,eID;
		sID = tempE->getSID();
		eID = tempE->getEID();
		double x1,x2,y1,y2,z1,z2;
		Node *tempN1 = NodeData+sID;
		Node *tempN2 = NodeData+eID;
		/* 二つの節点の座標ゲット */
		x1 = tempN1->getX();
		y1 = tempN1->getY();
		z1 = tempN1->getZ();
		x2 = tempN2->getX();
		y2 = tempN2->getY();
		z2 = tempN2->getZ();
		double xg = (x1+x2) / 2.0;
		double yg = (y1+y2) / 2.0;
		double zg = (z1+z2) / 2.0;
		/*  */
		bool zero = (xg < 1.0e-6 || yg < 1.0e-6);
		bool h_pos =  (xg > X_WIDTH - 1.0e-6 || yg > Y_WIDTH - 1.0e-6 || zg > Z_WIDTH - 1.0e-6);
		if(zero || h_pos){
			tempE->setZeroEdge(true);
			A0_size++;
		}else{
			tempE->setZeroEdge(false);
		}
	}
}
