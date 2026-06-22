#include "VoxelMesh.h"

/*//=======================================================
  // ● 各要素の材料情報作成
  //=======================================================*/
void VoxelMesh::set_Tvec(){
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

		bool incoil = ( (xg < COIL_MAX) && (yg < COIL_MAX) && (zg < COIL_Z) );
		/*  */
		if(incoil){
			double dz = z2 - z1;
			/* Z向きでなければTは無し */
			if(dz < 1.0e-10){
				tempE->setT_Edge(false);
			/* Z向きならTあり */
			}else{
				tempE->setT_Edge(true);
				setTValue(tempE, xg, yg, zg, dz);
			}
		}
	}
}

/*//=======================================================
  // ● 電流ポテンシャルの値計算
  //=======================================================*/
void VoxelMesh::setTValue(Edge *tempE, double xg, double yg, double zg, double dz){
	const double tMax = TOTAL_CURRENT / COIL_Z;
	double t_val;
	if((xg < COIL_MIN) && (yg < COIL_MIN)){
		t_val = dz*tMax;
	}else if((xg < COIL_MAX) && (yg < COIL_MIN)){
		t_val = ( dz*tMax*(1.0 -(xg - COIL_MIN)/(COIL_MAX -COIL_MIN)) );
	}else if((yg < COIL_MAX) && (xg < COIL_MIN)){
		t_val = ( dz*tMax*(1.0 -(yg - COIL_MIN)/(COIL_MAX -COIL_MIN)) );
	}else{
		double x = dz*tMax*(1.0 -(xg - COIL_MIN)/(COIL_MAX -COIL_MIN)) * (1.0 -(yg - COIL_MIN)/(COIL_MAX -COIL_MIN));
		t_val = x;
	}
	tempE->set_Tvalue(t_val);
}

