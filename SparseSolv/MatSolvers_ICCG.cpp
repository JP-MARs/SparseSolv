
#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
ICCGソルバ
//=======================================================
//=======================================================
*/


/*--------------------------------------------------------------------------*/



/*//=======================================================
// ● ICCGで解く(右辺ノルム内部計算パターン)
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init){
	const slv_int size = size0;
	double norm = 0;
	for(slv_int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	return solveICCG(size, conv_cri, max_ite, accera, norm, matA, vecB, results, init);
}

/*//=======================================================
// ● ICCGで解く(入力右辺がVector)
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init){
	double* vecBa = new double[size0];
	double* results_a = new double[size0];
	double norm=0;
	for(slv_int i = 0 ; i < size0 ; i++){
		const double temp = vecB[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = results[i];
	}
	norm = sqrt(norm);
	bool bl = solveICCG(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results_a, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}

/*//=======================================================
// ● ICCGで解く(入力右辺がEigen)
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init){
	double* vecBa = new double[size0];
	double norm=0;
	for(slv_int i = 0 ; i < size0 ; i++){
		const double temp = vecB(i);
		vecBa[i] = temp;
		norm += temp*temp;
	}
	norm = sqrt(norm);
	bool bl = solveICCG(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results, init);
	delete[] vecBa;
	return bl;
}



/*//=======================================================
// ● ICCGで解く・外部実行本体
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init){
	/* 対角スケーリングあり */
	if(is_diag_scale){
		return solveICCG_diag(size0, conv_cri, max_ite, accera, normB, matA, vecB, results, init);
	}
	/* 通常 */
	/* コレスキー用スパース行列作成 */
	double* diagD = new double[size0];
	double accela_val = accera;
	SparseMat matL;
	/* 加速係数が負なら自動決定モードへ */
	if(accera < 0){
		auto_accel_determine(size0, accera, matA, diagD, matL);
	}else{
		matL = matA.IC_decomp(diagD, accera);
	}
	/* 確定 */
	SparseMat matL_tr = matL.trans();

	bool bl = solveICCG(size0, conv_cri, max_ite, normB, diagD, matA.matrix, matL.matrix, matL_tr.matrix, vecB, results, init);
	delete[] diagD;
	return bl;
}

/*  */

/*//=======================================================
// ● ICCG専用内部処理(対角スケーリングあり)
//=======================================================*/
bool MatSolvers::solveICCG_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double* vecB, double* results, bool init){
	/* 対角スケーリング */
	double* vecB2 = new double[size0];
	SparseMat matD = matA.diagScaling(vecB2, vecB);
	SparseMat matDAD = matD*matA*matD;

	double normB2 = 0;
	for(slv_int i = 0; i < size0; i++){
		normB2 += vecB2[i]*vecB2[i];
	}
	normB2 = sqrt(normB2);

	/* コレスキー用スパース行列作成 */
	double* diagD = new double[size0];
	double accela_val = accera;
	SparseMat matL;
	/* 加速係数が負なら自動決定モードへ */
	if(accera < 0){
		auto_accel_determine(size0, accera, matA, diagD, matL);
	}else{
		matL = matA.IC_decomp(diagD, accera);
	}
	/* 確定 */
	SparseMat matL_tr = matL.trans();

	/* 解く */
	bool bl = solveICCG(size0, conv_cri, max_ite, normB2, diagD, matDAD.matrix, matL.matrix, matL_tr.matrix, vecB2, results, init);

	/* 元に戻す */
	double* result_true = matD*results;
	for(slv_int i = 0; i < size0; i++){
		results[i] = result_true[i];
	}
	delete[] diagD;
	delete[] vecB2;		
	delete[] result_true;
	return bl;
}

/*========================================*/
/*========================================*/
/*//=======================================================
// ● ICCGで解く（本体）
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
	const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results, bool init){

	/* 要素確保 */
	double alpha;
	double beta;
	Eigen::VectorXd EvecP(size);
	Eigen::VectorXd EvecR(size);
	Eigen::VectorXd EvecLDV(size);
	Eigen::VectorXd EtempAP(size);


	/* 初期設定 */
	if(init){
		for(slv_int i = 0 ; i < size ; i++){
			results[i] = 0;
		}
	}

		
	slv_int* start_posA = new slv_int[size];
	slv_int* end_posA = new slv_int[size];
	matA.getCols(start_posA, end_posA);
	auto col_ptrA = matA.getColPtr();
	auto val_ptrA = matA.getValuePtr();
#ifdef OMP_USING_ICCG		
	#pragma omp parallel for
#endif
	for(slv_int ii = 0 ; ii < size ; ii++){
		const slv_int c_size = end_posA[ii];
		double ap_temp=0;
		for(slv_int j = start_posA[ii] ; j < c_size ; j++){
			ap_temp += val_ptrA[j] * results[col_ptrA[j]];
		}
		EtempAP(ii) = ap_temp;
		/* 初期残差等計算*/
		EvecR(ii) = vecB[ii] - ap_temp;
	}

	/* 前処理 */
	IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecP);
	EvecLDV = EvecP;

	/* 絶対収束判定値をセット */
	const double abs_conv_cri = (normB*conv_cri*0.9 < small_abs_conv_val ? small_abs_conv_val : normB*conv_cri*0.9);

	/* 最初から答えだったら何もしない */
	const double first_normR = EvecR.norm() / normB;
	if(first_normR < conv_cri*0.1 || first_normR*normB < abs_conv_cri*0.1){
		delete[] start_posA;
		delete[] end_posA;
		return true;
	}
	/* 残差正規化方法をセット */
	double normalizer = normB;
	if(conv_normalize_type == 1){
		normalizer = first_normR;
	}else if(conv_normalize_type == 2){
		normalizer = conv_normalize_const;
	}

	/* 初期値 */
	double rur0 = EvecR.dot(EvecLDV);

	/* 最良結果の保存用（フラグがonなら） */
	double* best_results=nullptr;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = new double[size];
	}

	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* AP計算 */

//#ifdef OMP_USING_ICCG		
//		#pragma omp parallel for
//#endif
/*		for(slv_int ii = 0 ; ii < size ; ii++){
			const slv_int c_size = end_posA[ii];
			double tempVal = 0;
			for(slv_int j = start_posA[ii] ; j < c_size ; j++){
				tempVal	+= val_ptrA[j] * EvecP(col_ptrA[j]);
			}
			EtempAP(ii) = tempVal;
		}*/
		EtempAP = matA.matrix * EvecP;
		/* α計算 */
		//const double temp = EvecR.dot(EvecLDV);
		//alpha = temp / temp2;
		const double temp2 = EvecP.dot(EtempAP);		
		alpha = rur0 / temp2;
		

		/* 解ベクトルと残差計算 */
#ifdef OMP_USING_ICCG		
		#pragma omp parallel for
#endif
		for(slv_int i = 0 ; i < size ; i++){
			results[i] += alpha * EvecP(i);
		}
		EvecR = EvecR - alpha * EtempAP;
		const double pure_norm = EvecR.norm();
		const double normR = pure_norm / normalizer;

		/* フラグがonなら、残差保存 */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* 収束判定 */
		if(normR < conv_cri || pure_norm < abs_conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << ", " << temp << ", " << temp2  <<std::endl;
			is_conv = true;
			break;
		//}else if( It % 100 == 0 ){
		//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		
		if(normR < best_resi_value){
			best_resi_value = normR;
			/* 最良値の更新(フラグがonなら) */		
			if(is_save_best){
				for(slv_int i = 0; i < size; i++){
					best_results[i] = results[i];
				}
			}
		}
		/* 発散判定１ */
		if(diverge_judge_type == 1){
			/* 最良値×val以下なら、発散カウント初期化 */
			if(normR < best_resi_value * bad_div_val){
				bad_counter = 0;
			}
			/* 最良値×val以上なら、発散カウント＋ */
			if(normR >= best_resi_value * bad_div_val){
				bad_counter++;
			}
			/* 発散カウントが閾値オーバー＝発散扱いで終わる */
			if(bad_counter >= bad_div_count_thres){
				is_conv = false;
				break;
			}
		}

		/* v=(LDLT)-1rk　を計算 */
		IC_frbc_process(size, matL, matL_tr, diagD, EvecR, EvecLDV);
		/* β計算 */
		const double rur1 = EvecR.dot(EvecLDV);
		beta = rur1 / rur0;
		rur0 = rur1;
		//beta = EvecR.dot(EvecLDV);
		//beta /= temp;
		/* P計算 */
		EvecP = beta*EvecP + EvecLDV;
	}
	delete[] start_posA;
	delete[] end_posA;
	if(!is_conv){
		std::cout << "not Convergence!!! " << std::endl;
		/* 最良値を代入(フラグがonなら) */
		if(is_save_best){
			for(slv_int i = 0; i < size; i++){
				results[i] = best_results[i];
			}
			delete[] best_results;
		}
	}
	return is_conv ;
}

/* end of namespace */
};

