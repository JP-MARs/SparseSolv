/*
☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
　　拡大係数行列計算系統メソッド定義群
☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★
*/
#include "VoxelMesh.h"


/*//=======================================================
  // ● 行列作成に使う係数設定ルーチン
  //=======================================================*/
void VoxelMesh::setMatCoefficients(int k, double& sk, double& rk, double& tk){
	double Data_rk[] = {0,  0,  0,  0, 1,  1, -1, -1, 1, -1, -1,  1};
	double Data_sk[] = {1, -1, -1,  1, 0,  0,  0,  0, 1,  1, -1, -1};
	double Data_tk[] = {1,  1, -1, -1, 1, -1, -1,  1, 0,  0,  0,  0};
	rk = Data_rk[k];
	sk = Data_sk[k];
	tk = Data_tk[k];
}

/*//=======================================================
  // ● 右辺ベクター作成
  //=======================================================*/
void VoxelMesh::make_right_vector(double *vecB){
	const int size = edge_num - A0_size;
	int k_id, m_id;
	int *Edge_ids;
	double rm,sm,tm;
	for(int i = 0 ; i < size ; i++){
		vecB[i] = 0;
	}
	/* ループ開始 */
	Element *Ele_ptr = elements;
	double temp=0;
	for(int i = 0 ; i < element_num ; i++){
		/* ｉ番目の要素のどの辺もＴが無いなら、さようなら */
		int mat_coil = Ele_ptr->getMatType();
		if( mat_coil != MAT_COIL ){
			Ele_ptr++;
			continue;
		}
		/* ｉ番目要素の辺ゲット */
		Edge_ids = Ele_ptr->getEdges();
		/* ｉ番目要素の辺の数だけループ */
		for(int k = 0 ; k < 12 ; k++){
			/* 辺ｋでの情報セット */
			k_id = Edge_ids[k];
			int k_sub = EdgeData[k_id].getSubID();
			bool isZero = EdgeData[k_id].isZeroEdge();
			if(isZero){
				continue;
			}
			temp=0;
			/* 辺ｋの位置で分岐 */
			if(k < 4){
				for(int m = 4 ; m < 8 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += -tm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
				for(int m = 8 ; m < 12 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += sm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
			/* 分岐2 */
			}else if(k > 3 && k < 8){
				for(int m = 0 ; m < 4 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += tm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
				for(int m = 8 ; m < 12 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += -rm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
			/* 分岐3 */
			}else{
				for(int m = 0 ; m < 4 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += -sm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
				for(int m = 4 ; m < 8 ; m++){
					m_id = Edge_ids[m];
					bool is_T_e = EdgeData[m_id].isT_Edge();
					if( is_T_e ){
						setMatCoefficients(m, sm, rm, tm);
						temp += rm * EdgeData[m_id].get_Tvalue() / 8.0;
					}
				}
			}			
			/* 右辺ベクトルに代入 */
			vecB[k_sub] += temp;
		}
		Ele_ptr++;
	}
}
/*//=======================================================
  // ● 全体係数マトリックス作成
  //=======================================================*/
void VoxelMesh::make_matrix(SRLfem::SparseMat& matA){
	const int size = edge_num - A0_size;
	int k_id, m_id, k_sub,m_sub;
	double rk,sk,tk, rm,sm,tm;
	int *Edge_ids;
	double temp, L_x, L_y, L_z;
	/* 磁気抵抗率テンソル（v[0]:Vxx, v[1]:Vxy, v[2]:Vxz, v[3]:Vyy, v[4]:Vyz, v[5]:Vzz） */
	double v_tensor[9];

	Element *Ele_ptr = elements;
	for(int i = 0 ; i < element_num ; i++){
		/* ｉ番目要素の辺ゲット */
		Edge_ids = Ele_ptr->getEdges();
		/* 要素の長さゲット */
		L_x = ( Ele_ptr->getMaxX() > CHANGE_X ? x_length2/2.0 : x_length1/2.0 );
		L_y = ( Ele_ptr->getMaxY() > CHANGE_Y ? y_length2/2.0 : y_length1/2.0 );
		L_z = ( Ele_ptr->getMaxZ() > CHANGE_Z ? z_length2/2.0 : z_length1/2.0 );
		/* 磁気抵抗率ゲット */
		calc_v_tensor(v_tensor, Ele_ptr);
		/* ｋ、ｍ　ループ */
		for(int k = 0 ; k < 12 ; k++){
			k_id = Edge_ids[k];
			/* ゼロの既知辺なら、さようなら */
			if( EdgeData[k_id].isZeroEdge() ) continue;
			for(int m = k ; m < 12 ; m++){
				m_id = Edge_ids[m];
				/* ゼロの既知辺なら、さようなら */
				if( EdgeData[m_id].isZeroEdge() ) continue;
				/* 係数セット */
				setMatCoefficients(k, sk, rk, tk);
				setMatCoefficients(m, sm, rm, tm);
				temp = 0;
				/* (rotNm)・(rotNk) 計算 */
				if(k < 4){
					if(m < 4){
						temp =  v_tensor[4] * L_y*tk*tm * (2.0 + 2.0*sk*sm/3.0) / 16.0 / L_x / L_z;
						temp -= v_tensor[5] * sm * tk / 8.0 / L_x;
						temp -= v_tensor[7] * tm * sk / 8.0 / L_x;
						temp += v_tensor[8] * L_z*sk*sm * (2.0 + 2.0*tk*tm/3.0) / 16.0 / L_x / L_y;
					}else if( m > 3 && m < 8){
						temp = -v_tensor[3] * tk*tm / 8.0 / L_z;
						temp += v_tensor[5] * tk*rm / 8.0 / L_x;
						temp += v_tensor[6] * tm*sk / 8.0 / L_y;
						temp -= v_tensor[8] * L_z*sk*rm * (2.0 + 2.0*tk*tm/3.0) / 16.0 / L_x / L_y;
					}else{
						temp =  v_tensor[3] * tk*sm / 8.0 / L_z;
						temp -= v_tensor[6] * sk*sm / 8.0 / L_y;
						temp += v_tensor[7] * rm*sk / 8.0 / L_x;
						temp -= v_tensor[4] * L_y*tk*rm * (2.0 + 2.0*sk*sm/3.0) / 16.0 / L_x / L_z;
					}
				}else if(k > 3 && k < 8){
					if(m > 3 && m < 8){
						temp =  v_tensor[0] * L_x*tk*tm * (2.0 + 2.0*rk*rm/3.0) / 16.0 / L_y / L_z;
						temp -= v_tensor[2] * rm * tk / 8.0 / L_y;
						temp -= v_tensor[6] * tm * rk / 8.0 / L_y;
						temp += v_tensor[8] * L_z*rk*rm * (2.0 + 2.0*tk*tm/3.0) / 16.0 / L_x / L_y;
					}else{
						temp =  v_tensor[1] * tk*rm / 8.0 / L_z;
						temp += v_tensor[6] * rk*sm / 8.0 / L_y;
						temp -= v_tensor[7] * rm*rk / 8.0 / L_x;
						temp -= v_tensor[0] * L_x*tk*sm * (2.0 + 2.0*rk*rm/3.0) / 16.0 / L_y / L_z;
					}
				}else{
					temp =  v_tensor[0] * L_x*sk*sm * (2.0 + 2.0*rk*rm/3.0) / 16.0 / L_y / L_z;
					temp -= v_tensor[1] * rm * sk / 8.0 / L_z;
					temp -= v_tensor[3] * sm * rk / 8.0 / L_z;
					temp += v_tensor[4] * L_y*rk*rm * (2.0 + 2.0*sk*sm/3.0) / 16.0 / L_x / L_z;
				}
				/**/
				k_sub = EdgeData[k_id].getSubID();
				m_sub = EdgeData[m_id].getSubID();
				/* 代入～！ */
				matA.add(k_sub, m_sub, temp);
				if(k != m){
					matA.add(m_sub, k_sub, temp);
				}
			}
		}
		Ele_ptr++;
	}
}

