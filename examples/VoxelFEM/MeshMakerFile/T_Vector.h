/*
//=======================================================
// ■ 電流ベクトルポテンシャル設定専門staticクラス
//=======================================================
// 
//=======================================================*/
class T_Vector{
public:
	static bool is_T_Edge(Edge *tempE, Node_Mesh *nodes);
	static double setTValue(Edge *tempE, Node_Mesh *nodes);
};
/*//=======================================================
  // ● ゼロの既知ポテンシャルを持つかどうか
  //=======================================================*/
bool T_Vector::is_T_Edge(Edge *tempE, Node_Mesh *nodes){
	/* 始点と終点ゲット */
	int sID,eID;
	sID = tempE->getSID();
	eID = tempE->getEID();
	double x1,x2,y1,y2,z1,z2;
	Node_Mesh *tempN1 = nodes+sID;
	Node_Mesh *tempN2 = nodes+eID;
	/* 二つの節点の座標ゲット */
	x1 = tempN1->getPosX();
	y1 = tempN1->getPosY();
	z1 = tempN1->getPosZ();
	x2 = tempN2->getPosX();
	y2 = tempN2->getPosY();
	z2 = tempN2->getPosZ();
	double xg = (x1+x2) / 2.0;
	double yg = (y1+y2) / 2.0;
	double zg = (z1+z2) / 2.0;
	/* コイルの四角形ないにある、かつ… */
	if((xg < COIL_MAX) && (yg < COIL_MAX) && (zg < COIL_Z)){
		double dz = z2 - z1;
		/* Z向きでなければTは無し */
		if(dz < 1.0e-10){
			return false;
		/* Z向きならTあり */
		}else{
			return true;
		}
	}
	return false;
}
/*//=======================================================
  // ● 電流ポテンシャルの値計算
  //=======================================================*/
double T_Vector::setTValue(Edge *tempE, Node_Mesh *nodes){
	/* 始点と終点ゲット */
	int sID,eID;
	sID = tempE->getSID();
	eID = tempE->getEID();
	double x1,x2,y1,y2,z1,z2;
	Node_Mesh *tempN1 = nodes+sID;
	Node_Mesh *tempN2 = nodes+eID;
	/* 二つの節点の座標ゲット */
	x1 = tempN1->getPosX();
	y1 = tempN1->getPosY();
	z1 = tempN1->getPosZ();
	x2 = tempN2->getPosX();
	y2 = tempN2->getPosY();
	z2 = tempN2->getPosZ();
	const double xg = (x1+x2) / 2.0;
	const double yg = (y1+y2) / 2.0;
	const double zg = (z1+z2) / 2.0;
	const double tMax = TOTAL_CURRENT / (double)(COIL_Z);
	const double dz = z2 - z1;

	if((xg < COIL_MIN) && (yg < COIL_MIN)){
		return dz*tMax;
	}else if((xg < COIL_MAX) && (yg < COIL_MIN)){
		return( dz*tMax*(1.0 -(xg - COIL_MIN)/(COIL_MAX -COIL_MIN)) );
	}else if((yg < COIL_MAX) && (xg < COIL_MIN)){
		return( dz*tMax*(1.0 -(yg - COIL_MIN)/(COIL_MAX -COIL_MIN)) );
	}else{
		double x = dz*tMax*(1.0 -(xg - COIL_MIN)/(COIL_MAX -COIL_MIN)) * (1.0 -(yg - COIL_MIN)/(COIL_MAX -COIL_MIN));
		return x;
	}
	return 0;
}

