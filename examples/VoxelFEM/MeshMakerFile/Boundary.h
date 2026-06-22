/*
//=======================================================
// ■ ゼロ境界かどうかを判定する専門staticクラス
//=======================================================
// 
//=======================================================*/
class Boundary{
public:
	static bool isZeroEdge(Edge *tempE, Node_Mesh *nodes);
};
/*//=======================================================
  // ● ゼロの既知ポテンシャルを持つかどうか
  //=======================================================*/
bool Boundary::isZeroEdge(Edge *tempE, Node_Mesh *nodes){
	/* 始点と終点ゲット */
	int sID,eID;
	sID = tempE->getSID();
	eID = tempE->getEID();
	bool bl;


	//return false;

	/* X座標が０か端の時、ゼロ */
	bl = ( nodes[sID].getX() == 0 );
	bl &= ( nodes[eID].getX() == 0 );
	if(bl) return true;
	bl = ( nodes[sID].getX() == TOTAL_MESH_X );
	bl &= ( nodes[eID].getX() == TOTAL_MESH_X );
	if(bl) return true;
	/* Y座標が０か端の時、ゼロ */
	bl = ( nodes[sID].getY() == 0 );
	bl &= ( nodes[eID].getY() == 0 );
	if(bl) return true;
	bl = ( nodes[sID].getY() == TOTAL_MESH_Y );
	bl &= ( nodes[eID].getY() == TOTAL_MESH_Y );
	if(bl) return true;
	/* Z座標が上端の時、ゼロ（０は自然境界なのでいらん） */
	bl = ( nodes[sID].getZ() == TOTAL_MESH_Z );
	bl &= ( nodes[eID].getZ() == TOTAL_MESH_Z );
	if(bl) return true;
	return false;
}








#ifdef FFFFFFFFFFFFFFF/*
//=======================================================
// ■ ICCG
//=======================================================
// ICCG法ソルバ (実体なしstatic使用がメイン・スパース非対応)
//=======================================================*/
class ICCG {
private:
	static int size;																										/* 行列サイズ */
	static double conv;																										/* 収束判定 */
	static int max_iteration;																								/* 最大反復回数 */
	static void Cholesky(int *clmn_size, int **column, double **matrixA, int *clmnL_size, int **columnL, double **matrixL, int *clmL_tr_size, int **columnL_tr, double **matrixL_tr, double *diagD);	/* 不完全コレスキー分解！！ */
	static void preProcess(int *clmnL_size, int **columnL, double **matrixL, int *clmnL_tr_size, int **columnL_tr, double **matrixL_tr, double *diagD, double *vecR, double *vec);					/* 前処理 */
	static double inner_products(const double *a, const double *b);																		/* 内積計算 */
public:
	static void setSize(int s){size = s;};																								/* サイズセット */
	static void ICCGsolve(int s, int *clmn_size, int **column, double **matrixA, double *vecB, double *results);						/* ICCGで解く */
	static void ICCGsolve(int *clmn_size, int **column, double **matrixA, double *vecB, double *results);								/* ICCGで解くver2 */
};
int ICCG::size;
double ICCG::conv = EPS_SOLVER;
int ICCG::max_iteration = MAX_ITR_ICCG;
/*//=======================================================
  // ● ICCGで解く
  //=======================================================*/
void ICCG::ICCGsolve(int s, int *clmn_size, int **column, double **matrixA, double *vecB, double *results){
	size = s;
	ICCGsolve(clmn_size, column, matrixA, vecB, results);
}
/*//=======================================================
  // ● ICCGで解く
  //=======================================================*/
void ICCG::ICCGsolve(int *clmn_size, int **column, double **matrixA, double *vecB, double *results){
	int *clmnL_size = new int[size];
	int **columnL = new int*[size];
	double **matrixL = new double*[size];
	int *clmnL_tr_size = new int[size];
	int **columnL_tr = new int*[size];
	double **matrixL_tr = new double*[size];
	/* 行列Lの各行のサイズを計算＝Aの下三角分 */
	for(int i = 0 ; i < size ; i++){
		int inv_size = 0;
		const int C_size =clmn_size[i];
		int *C_ptr = column[i];
		for(int l = 0 ; l < C_size ; l++){
			if( *C_ptr == i){
				clmnL_size[i] = l+1;
			}
			if(*C_ptr < i){
				inv_size++;
			}
			C_ptr++;
		}	
		clmnL_tr_size[i] = C_size-inv_size;
	}
	/* 行列Lを確保 */
#ifdef OMP_USING
	#pragma omp parallel for
#endif
	for(int i = 0 ; i < size ; i++){
		columnL[i] = new int[clmnL_size[i]];
		matrixL[i] = new double[clmnL_size[i]];
		columnL_tr[i] = new int[clmnL_tr_size[i]];
		matrixL_tr[i] = new double[clmnL_tr_size[i]];
	}
	/* Lの列数を代入 */
#ifdef OMP_USING
	#pragma omp parallel for
#endif
	for(int i = 0 ; i < size ; i++){
		const int C_size = clmnL_size[i];
		for(int l = 0 ; l < C_size ; l++){
			columnL[i][l] = column[i][l];
		}
	}
	/* 要素確保 */
	double *diagD = new double[size];
	double *vecP = new double[size];
	double *vecR = new double[size];
	double *vecLDV = new double[size];
	double alpha;
	double beta;
	double *tempAP = new double[size];
	int j;
	
	/* 初期設定 */
	for(int i = 0 ; i < size ; i++){
		results[i] = 0;
	}
	/* コレスキー分解 */
	Cholesky(clmn_size, column, matrixA, clmnL_size, columnL, matrixL, clmnL_tr_size, columnL_tr, matrixL_tr, diagD);
	/* 前処理 */
	preProcess(clmnL_size, columnL, matrixL, clmnL_tr_size, columnL_tr, matrixL_tr, diagD, vecB, results);
	double normB=0;
	double *temp_ptr1 = tempAP;
	double *temp_ptr2 = vecB;
	for(int ii = 0 ; ii < size ; ii++){
		*temp_ptr1 = 0;//tempAP[ii] = 0;
		normB += (*temp_ptr2) * (*temp_ptr2);//vecB[ii]*vecB[ii];
		const int c_size = clmn_size[ii];
		double *A_pr = matrixA[ii];
		int *C_pr = column[ii];
		for(j = 0 ; j < c_size ; j++){
			//*(tempAP+ii) += (*A_pr) * (*(results+the_c));
			*temp_ptr1 += (*A_pr) * results[(*C_pr)];
			A_pr++;
			C_pr++;
		}
		/* 初期残差等計算*/
		vecR[ii] = vecB[ii] - tempAP[ii];
		if(isnan(vecR[ii])){
			cout << ii << ", " << vecR[ii] << endl;
			exit(1);
		}
		temp_ptr1++;
		temp_ptr2++;
	}
	normB = sqrt(normB);
	if(isnan(normB)){
		cout << "Norm B is " << normB << endl;
		exit(1);
	}
	/* 前処理 */
	preProcess(clmnL_size, columnL, matrixL, clmnL_tr_size, columnL_tr, matrixL_tr, diagD, vecR, vecP);
#ifdef OMP_USING		
	#pragma omp parallel for
#endif
	for(int i = 0 ; i < size ; i++){
		vecLDV[i] = vecP[i];
	}
	double normR=0;
	double temp, temp2;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_iteration ; It++){
		/* AP計算 */
		double *temp_ptr1 = tempAP;
#ifdef OMP_USING		
		#pragma omp parallel for private(j)
#endif
		for(int ii = 0 ; ii < size ; ii++){
			*temp_ptr1 = 0;//*(tempAP+ii) = 0;
			const int c_size = clmn_size[ii];
			double *A_pr = matrixA[ii];
			int *C_pr = column[ii];
			for(j = 0 ; j < c_size ; j++){
				//*(tempAP+ii) += (*A_pr) * (*(vecP+the_c));
				*temp_ptr1 += (*A_pr) * vecP[(*C_pr)];
				A_pr++;
				C_pr++;
			}
			temp_ptr1++;
		}
		/* α計算 */
		temp = 0;
		temp2 = 0;
		/* 内積計算 */
#ifdef OMP_USING		
		#pragma omp parallel for reduction(+: temp) reduction(+: temp2)
#endif
		for(int i = 0 ; i < size ; i++){
			temp += vecR[i] * vecLDV[i];
			temp2 += vecP[i] * tempAP[i];
			if(isnan(vecP[i])){
				cout << "vecP isnan " << i << endl;
				exit(1);
			}
		}
		alpha = temp / temp2;
		/* 解ベクトルと残差計算 */
		normR=0;
#ifdef OMP_USING		
		#pragma omp parallel for reduction (+:normR)
#endif
		for(int i = 0 ; i < size ; i++){
			results[i] += alpha * vecP[i];
			vecR[i] -= alpha * tempAP[i];
			normR += vecR[i] * vecR[i];
		}
		normR = sqrt(normR);

		if(isnan(normR)){
			cout << "Norm R is " << normR << endl;
			exit(1);
		}

		normR /= normB;
		if(normR < conv){
			cout << "Solved!!! -- " << normR  << " " << It << " \r\n";
			break;
		}else{
			cout << "NORMR is -- " << normR  << " " << It << " \r\n";
		}
		/* v=(LDLT)-1rk　を計算 */
		preProcess(clmnL_size, columnL, matrixL, clmnL_tr_size, columnL_tr, matrixL_tr, diagD, vecR, vecLDV);
		/* β計算 */
		beta = inner_products(vecR, vecLDV);
		beta /= temp;
		/* P計算 */
#ifdef OMP_USING		
		#pragma omp parallel for
#endif
		for(int i = 0 ; i < size ; i++){
			vecP[i] *= beta;
			vecP[i] += vecLDV[i];
		}
	}
	if(It >= max_iteration-1){
		cout << "not Convergence!!! -- " << normR << "\r\n";
	}
	delete[] matrixL;
	delete[] columnL;
	delete[] columnL_tr;
	delete[] matrixL_tr;
	delete[] diagD;
	delete[] vecP;
	delete[] vecR;
	delete[] vecLDV;
	delete[] tempAP;
}
/*//=======================================================
  // ● 内積計算
  //=======================================================*/
double ICCG::inner_products(const double *a, const double *b){
	double x = 0;
#ifdef OMP_USING		
	#pragma omp parallel for reduction (+:x)
#endif
	for(int i = 0 ; i < size ; i++){
		x += a[i] * b[i];
	}
	return x;
}
/*//=======================================================
  // ● 不完全コレスキー分解！！
  //=======================================================*/
void ICCG::Cholesky(int *clmn_size, int **column, double **matrixA, int *clmnL_size, int **columnL, double **matrixL, int *clmL_tr_size, int **columnL_tr, double **matrixL_tr, double *diagD){
	double s;
	int *J_count = new int[size];
	int *temp_ptr = J_count;
	for (int i = 0; i < size; i++) {
		*temp_ptr = 0;
		temp_ptr++;//	J_count[i] = 0;
	}
	/* LとDを求める */
	for (int i = 0; i < size; i++) {
		/* i行目にある非ゼロの数をゲット */
		const int c_size = clmn_size[i];
		double *A_ptr = matrixA[i];
		int *C_ptr = column[i];
		/* i行目の非ゼロの回数だけ列ループをまわす */
		for (int jj = 0; jj < c_size; jj++) {
			/* 現在の列番号jを取得 */
			const int j = *C_ptr;			
			//const int inv_j_size = clmL_tr_size[j];
			if(j > i) break;
			//if( fabs( *A_ptr ) > 1.0e-50 ){
				s = *A_ptr;
				/* i行目の列の非ゼロ列をサーチ */
				const int L_ksize = clmnL_size[i];
				double *L_ptr = matrixL[i];
				int *K_ptr = columnL[i];
				for(int kk = 0 ; kk < L_ksize ; kk++){
					/* i行目の列の非ゼロ列の列番号ｋを取得 */
					const int k = *K_ptr;
					/* kがj-1以上ならループエンド */
					if( k >= j-1 ) break;
					/* ｊ行目で、列番号kになる位置をサーチ */
					const int LJ_ksize = clmnL_size[j];
					double *LJ_ptr = matrixL[j];
					int *KJ_ptr = columnL[j];
					for(int l = 0 ; l < LJ_ksize ; l++){
						if( *KJ_ptr >= j-1 ) break;
						if( *KJ_ptr == k ){
							s -= (*L_ptr) * (*LJ_ptr) * diagD[k];
						}
						LJ_ptr++;
						KJ_ptr++;
					}
					L_ptr++;
					K_ptr++;
				}
				columnL[i][jj] = j;
				matrixL[i][jj] = s;
				columnL_tr[j][J_count[j]] = i;
				matrixL_tr[j][J_count[j]] = s;
				J_count[j]++;
//				cout << "Mat L is " << i << ", " << jj << ", " << matrixL[i][jj] << ", " << matrixA[i][jj] <<endl;
//				if(isnan( matrixL[i][jj] ) || isinf(matrixL[i][jj])) exit(1);
			//}else{
			//	cout << "!?!?!?!?!?!?! ==> "<< i << ", " << columnL[i][jj] << " & " << matrixL[i][jj] << endl;
			//}
			A_ptr++;
			C_ptr++;
		}
		const int last = clmnL_size[i] - 1;
		double x = matrixL[i][last];
		if(fabs(x) < 1.0e-20){
			cout << "TOOOO SMALLLL!!!!!!!!!!!! " << x << ", " <<  DBL_EPSILON<< endl;
			diagD[i] = 1.0 / (DBL_EPSILON);
		}else{
			diagD[i] = 1.0 / x;
		}
	}
	delete[] J_count; 
/*
	for (int j = 0; j < size; j++) {
		int j_size = clmL_tr_size[j];
		for (int i = 0; i < j_size; i++) {
			int the_i = columnL_tr[j][i];
			double temp1 = matrixL_tr[j][i];
			int c_size = clmnL_size[the_i];
			for (int jj = 0; jj < c_size; jj++) {
				if(columnL[the_i][jj] == j){
					double temp2 = matrixL[the_i][jj];

					cout << "LLL is " << j << ", " << the_i << " & " << the_i << ", " << columnL[the_i][jj] << " --> " << temp1 << ", " << temp2 << endl;
					if(fabs(temp1-temp2) > 1.0e-10){
						exit(1);
					}
				}
			}
		}
	}
*/
}
/*//=======================================================
  // ● 前処理
  //=======================================================*/
void ICCG::preProcess(int *clmnL_size, int **columnL, double **matrixL, int *clmnL_tr_size, int **columnL_tr, double **matrixL_tr, double *diagD, double *vecR, double *vec){
	double s;
	double *vec_ptr = vec;
	/* 第一方程式計算 */
	for (int i = 0; i < size; i++) {
		s = vecR[i];//*vecR;//
		const int c_size = clmnL_size[i] - 1;//(*clmnL_size) - 1;//clmnL_size[i] - 1;
		double *L_ptr = matrixL[i];;//*matrixL;//matrixL[i];
		int *C_ptr = columnL[i];//*columnL;//columnL[i];
		int j;
		for (int jj = 0; jj < c_size; jj++) {
			j = *C_ptr ;
			s -= (*L_ptr) * vec[j];
			L_ptr++;
			C_ptr++;
		}
		//s *= diagD[i];
		vec[i] = s * diagD[i];// / (*L_ptr);
		if( isnan(vec[i]) || isinf(vec[i])){
			cout << "NANANA" << vec[i] << endl;
			exit(1);
		}
		//*vec_ptr = s / (*L_ptr);
		//vec[i] = s / ( matrixL[i][c_size] );
		//vecR++;
		//clmnL_size++;
		//matrixL++;
		//columnL++;
		//vec_ptr++;
	}
	vecR -= size;
	clmnL_tr_size+=size-2; columnL_tr+=size-2; matrixL_tr+=size-2; diagD+=size-2;
	vec_ptr = vec+size-2;
	/* 第二方程式計算 */
	for (int i = size-2; i > -1; i--) {	
		s = 0;
		const int c_size = *clmnL_tr_size;//clmnL_tr_size[i]
		double *L_ptr = (*matrixL_tr)+1; //matrixL_tr[i];		
		int *C_ptr = (*columnL_tr)+1; //columnL_tr[i];
		int the_row;
		for (int j = 1; j < c_size; j++) {
			the_row = *C_ptr ;
			s += (*L_ptr) * vec[the_row];
			L_ptr++;
			C_ptr++;
		}
		s *= *diagD; //diagD[i]; 
		if( isnan(vec[i]) || isinf(vec[i])){
			cout << "NANANA222" << vec[i] << endl;
			exit(1);
		}
		*vec_ptr -= s; //vec[i] -= s;// * diagD[i];
		clmnL_tr_size--;
		matrixL_tr--;
		columnL_tr--;
		vec_ptr--;
		diagD--;
	}
}

#endif
