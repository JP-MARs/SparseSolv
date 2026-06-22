
#include "Element.h"

Node *Element::NodeData=NULL;
Edge *Element::EdgeData=NULL;
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Element::Element(){
	nodes = new int[8];
	edges = new int[12];
}
Element::~Element(){
	delete[] nodes;
	delete[] edges;
}
/*//=======================================================
  // ● データのセッタ
  //=======================================================*/
void Element::set(int xx, int yy, int zz, int *n, int *e, int id){
	ID = id;
	x = xx;
	y = yy; 
	z = zz;
	for(int i = 0 ; i < 8 ; i++){
		nodes[i] = n[i];
	}
	for(int i = 0 ; i < 12 ; i++){
		edges[i] = e[i];
	}
}
/*//=======================================================
  // ● 要素の中心位置返し
  //=======================================================*/
void Element::getGrav(double *ggg){
	double g_x = 0;
	double g_y = 0;
	double g_z = 0;
	for(int i = 0 ; i < 8 ;i++){
		g_x += NodeData[nodes[i]].getX();
		g_y += NodeData[nodes[i]].getY();
		g_z += NodeData[nodes[i]].getZ();
	}
	g_x /= 8.0; g_y /= 8.0; g_z /= 8.0;
	ggg[0] = g_x; ggg[1] = g_y; ggg[2] = g_z;
}
/*//=======================================================
  // ● 要素中心の磁束密度計算
  //=======================================================*/
void Element::calcB(double *B, double *A){
	double Data_rk[] = {0,  0,  0,  0, 1,  1, -1, -1, 1, -1, -1,  1};
	double Data_sk[] = {1, -1, -1,  1, 0,  0,  0,  0, 1,  1, -1, -1};
	double Data_tk[] = {1,  1, -1, -1, 1, -1, -1,  1, 0,  0,  0,  0};
	double ggg[3];
	getGrav(ggg);
	double vec_B[3]={0,0,0};
	double rot[3];
	double r,s,t;
	for(int i = 0 ; i < 12 ; i++){
		r = Data_rk[i];
		s = Data_sk[i];
		t = Data_tk[i];
		/* まず中心の回転計算 */
		calcRot(i, r, s, t, rot);
		/* 辺番号から、磁気ベクトルポテンシャルゲット */
		int subID = EdgeData[edges[i]].getSubID();
		double edge_A;
		if(subID < 0){
			edge_A = 0;
		}else{
			/* 代入 */
			edge_A = A[subID];
			vec_B[0] += edge_A * rot[0];
			vec_B[1] += edge_A * rot[1];
			vec_B[2] += edge_A * rot[2];		
		}
	}
	/* 要素中心の座標と磁束密度代入して終わり */
	B[0] = ggg[0];
	B[1] = ggg[1];
	B[2] = ggg[2];
	B[3] = vec_B[0];
	B[4] = vec_B[1];
	B[5] = vec_B[2];
}
/*//=======================================================
  // ● 要素中心の回転計算
  //=======================================================*/
void Element::calcRot(int k, double r, double s, double t, double *rot){
	double L_x = getMaxX() - getMinX();
	double L_y = getMaxY() - getMinY();
	double L_z = getMaxZ() - getMinZ();
	L_x /= 2.0; L_y /= 2.0; L_z /= 2.0;
	/* 中心の回転 */
	if(k < 4){
		rot[0] = 0;
		rot[1] = t / 8.0 / L_x / L_z;
		rot[2] = -s / 8.0 / L_x / L_y;
	}else if(k > 3 && k < 8){
		rot[0] = -t / 8.0 / L_y / L_z;
		rot[1] = 0;
		rot[2] = r / 8.0 / L_x / L_y;
	}else{
		rot[0] = s / 8.0 / L_y / L_z;
		rot[1] = -r / 8.0 / L_x / L_z;
		rot[2] = 0;
	}
}




/*//=======================================================
  // ● 要素内の任意の位置のrot算出
  //=======================================================*/
void Element::calcRot(int k, double *pos_rst, double *rst, double *rot){
	/* pos_rst：計算位置   rts：係数ri, si, ti */
	double L_x = getMaxX() - getMinX();
	double L_y = getMaxY() - getMinY();
	double L_z = getMaxZ() - getMinZ();
	L_x /= 2.0; L_y /= 2.0; L_z /= 2.0;
	/* 中心の回転 */
	if(k < 4){
		rot[0] = 0;
		rot[1] = rst[2] * (1.0 + pos_rst[1]*rst[1]) / 8.0 / L_x / L_z;
		rot[2] = -rst[1] * (1.0 + pos_rst[2]*rst[2]) / 8.0 / L_x / L_y;
	}else if(k > 3 && k < 8){
		rot[0] = -rst[2] * (1.0 + pos_rst[0]*rst[0]) / 8.0 / L_y / L_z;
		rot[1] = 0;
		rot[2] = rst[0] * (1.0 + pos_rst[2]*rst[2]) / 8.0 / L_x / L_y;
	}else{
		rot[0] = rst[1] * (1.0 + pos_rst[0]*rst[0]) / 8.0 / L_y / L_z;
		rot[1] = -rst[0] * (1.0 + pos_rst[1]*rst[1]) / 8.0 / L_x / L_z;
		rot[2] = 0;
	}	
}

/*//=======================================================
  // ● 要素の磁束密度計算
  //=======================================================*/
void Element::calcB(double *B, double *A, double *pos){
	double Data_rk[] = {0,  0,  0,  0, 1,  1, -1, -1, 1, -1, -1,  1};
	double Data_sk[] = {1, -1, -1,  1, 0,  0,  0,  0, 1,  1, -1, -1};
	double Data_tk[] = {1,  1, -1, -1, 1, -1, -1,  1, 0,  0,  0,  0};
	double vec_B[3]={0,0,0};
	double rot[3];
	double rst[3];
	for(int i = 0 ; i < 12 ; i++){
		rst[0] = Data_rk[i];
		rst[1] = Data_sk[i];
		rst[2] = Data_tk[i];
		/* まず中心の回転計算 */
		calcRot(i, pos, rst, rot);
		/* 辺番号から、磁気ベクトルポテンシャルゲット */
		int subID = EdgeData[edges[i]].getSubID();
		double edge_A;
		if(subID < 0){
			edge_A = 0;
		}else{
			/* 代入 */
			edge_A = A[subID];
			vec_B[0] += edge_A * rot[0];
			vec_B[1] += edge_A * rot[1];
			vec_B[2] += edge_A * rot[2];		
		}
	}
	B[0] = vec_B[0]; B[1] = vec_B[1]; B[2] = vec_B[2];
}

/*//=======================================================
  // ● 要素のJ
  //=======================================================*/
void Element::calcJ(double *J, double *pos){
	double Data_rk[] = {0,  0,  0,  0, 1,  1, -1, -1, 1, -1, -1,  1};
	double Data_sk[] = {1, -1, -1,  1, 0,  0,  0,  0, 1,  1, -1, -1};
	double Data_tk[] = {1,  1, -1, -1, 1, -1, -1,  1, 0,  0,  0,  0};
	double vec_J[3]={0,0,0};
	double rot[3];
	double rst[3];
	for(int i = 0 ; i < 12 ; i++){
		rst[0] = Data_rk[i];
		rst[1] = Data_sk[i];
		rst[2] = Data_tk[i];
		/* まず中心の回転計算 */
		calcRot(i, pos, rst, rot);
		bool isT = EdgeData[edges[i]].isT_Edge();
		if(isT){
			double t = EdgeData[edges[i]].get_Tvalue();
			vec_J[0] += t * rot[0];
			vec_J[1] += t * rot[1];
			vec_J[2] += t * rot[2];		
		}
	}
	J[0] = vec_J[0]; J[1] = vec_J[1]; J[2] = vec_J[2];
}

bool Element::isInclude(double *pos){
	double minX = getMinX(); double maxX = getMaxX();
	double minY = getMinY(); double maxY = getMaxY();
	double minZ = getMinZ(); double maxZ = getMaxZ();
	if(minX <= pos[0] && maxX >= pos[0]){
		if(minY <= pos[1] && maxY >= pos[1]){
			if(minZ <= pos[2] && maxZ >= pos[2]){
				return true;
			}
		}
	}
	return false;
}
