
#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
MRTRソルバ
//=======================================================
//=======================================================
*/


/*//=======================================================
// ● MRTRで解く(右辺ノルム内部計算パターン)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double *vecB, double *results, bool init){
	const slv_int size = size0;
	double norm = 0;
	for(slv_int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	return solveSGSMRTR(size, conv_cri, max_ite, norm, matA, vecB, results, init);
}


/*//=======================================================
// ● MRTRで解く(入力右辺がVector)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init){
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
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm, matA, vecBa, results_a, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}

/*//=======================================================
// ● MRTRで解く(入力右辺がEigen)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init){
	double* vecBa = new double[size0];
	double norm=0;
	for(slv_int i = 0 ; i < size0 ; i++){
		const double temp = vecB(i);
		vecBa[i] = temp;
		norm += temp*temp;
	}
	norm = sqrt(norm);
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm, matA, vecBa, results, init);
	delete[] vecBa;
	return bl;
}

/*//=======================================================
// ● MRTRで解く・外部実行本体
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init){
	/* 対角スケーリング */
	double* vecB2 = new double[size0];
	SparseMat matD = matA.diagScaling(vecB2, vecB);
	SparseMat matDAD = matD*matA*matD;

	double normB2 = 0;
	for(slv_int i = 0; i < size0; i++){
		normB2 += vecB2[i]*vecB2[i];
	}
	normB2 = sqrt(normB2);

	/* DADの下三角を取得 */
	SparseMat matL = matDAD.getMatLower();
	/* DADの上三角を取得 */
	SparseMat matL_tr = matL.trans();

	bool bl= solveSGSMRTR(size0, conv_cri, max_ite, normB, matDAD.matrix, matL.matrix, matL_tr.matrix, vecB2, results, init);
	delete[] vecB2;
	/* 元に戻す */
	double* result_true = matD*results;
	for(slv_int i = 0; i < size0; i++){
		results[i] = result_true[i];
	}
	delete[] result_true;	
	return bl;	
}

/*-----------------------------------------------*/



/*//=======================================================
// ● MRTRで解く(右辺ノルム内部計算パターン)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init){
	const slv_int size = size0;
	dcomplex norm = 0;
	for(slv_int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	double norm0 = abs(norm);
	double norm00 = sqrt(norm0);
	return solveSGSMRTR(size, conv_cri, max_ite, norm00, matA, vecB, results, init);
}


/*//=======================================================
// ● MRTRで解く(入力右辺がVector)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init){
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex* results_a = new dcomplex[size0];
	dcomplex norm=0;
	for(slv_int i = 0 ; i < size0 ; i++){
		const dcomplex temp = vecB[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = results[i];
	}
	double norm0 = abs(norm);
	double norm00 = sqrt(norm0);
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm00, matA, vecBa, results_a, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}

/*//=======================================================
// ● MRTRで解く(入力右辺がEigen)
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init){
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex norm=0;
	for(slv_int i = 0 ; i < size0 ; i++){
		const dcomplex temp = vecB[i];
		vecBa[i] = temp;
		norm += temp*temp;
	}
	double norm0 = abs(norm);
	double norm00 = sqrt(norm0);
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm00, matA, vecBa, results, init);
	delete[] vecBa;
	return bl;
}

/*//=======================================================
// ● MRTRで解く・外部実行本体
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init){
	/* 対角スケーリング */
	dcomplex* vecB2 = new dcomplex[size0];
	SparseMatC matD = matA.diagScaling(vecB2, vecB);
	SparseMatC matDAD = matD*matA*matD;

	dcomplex normB2 = 0;
	for(slv_int i = 0; i < size0; i++){
		normB2 += vecB2[i]*vecB2[i];
	}
	double norm0 = abs(normB2);
	double norm00 = sqrt(norm0);

	/* DADの下三角を取得 */
	SparseMatC matL = matDAD.getMatLower();
	/* DADの上三角を取得 */
	SparseMatC matL_tr = matL.trans();

	bool bl= solveSGSMRTR(size0, conv_cri, max_ite, norm00, matDAD.matrix, matL.matrix, matL_tr.matrix, vecB2, results, init);
	delete[] vecB2;
	/* 元に戻す */
	dcomplex* result_true = matD*results;
	for(slv_int i = 0; i < size0; i++){
		results[i] = result_true[i];
	}
	delete[] result_true;	
	return bl;	
}



/*========================================*/
/*========================================*/
/*//=======================================================
// ● MRTRで解く（本体）
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
	const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results, bool init){

	/* 要素確保 */
	Eigen::VectorXd EvecP = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd EvecRd(size);
	Eigen::VectorXd EvecU(size);
	Eigen::VectorXd EvecX(size);	
	Eigen::VectorXd EvecY(size);
	Eigen::VectorXd EvecARd(size);
	
	/* 初期設定 */
	if(init){
		for(slv_int i = 0 ; i < size ; i++){
			results[i] = 0;
			EvecX(i) = 0.0;
		}
	} else{
		for(int i = 0 ; i < size ; i++){
			EvecX(i) = results[i];
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
		/* 初期残差等計算*/
		EvecRd(ii) = vecB[ii] - ap_temp;
		EvecU(ii) = vecB[ii];
	}

	/* 絶対収束判定値をセット */
	const double abs_conv_cri = (normB*conv_cri*0.9 < small_abs_conv_val ? small_abs_conv_val : normB*conv_cri*0.9);

	/* 最初から答えだったら何もしない */
	const double first_normR = EvecRd.norm() / normB;
	if(first_normR < conv_cri*0.1 || first_normR*normB < abs_conv_cri*0.1){
		delete[] start_posA;
		delete[] end_posA;
		return true;
	}
	/* 収束判定用の補正|B| */
	Eigen::VectorXd tempVec = matL.matrix*EvecU;
	const double normBd = tempVec.norm();

	/* 残差正規化方法をセット */
	double normalizer = normBd;
	if(conv_normalize_type == 1){
		normalizer = first_normR;
	}else if(conv_normalize_type == 2){
		normalizer = conv_normalize_const;
	}

	/* 前処理rd=C^-1 * r */
	tempVec = EvecRd;
	fr_process(size, matL, tempVec, EvecRd);	
	/* y0 = -r0 */
	EvecY = -1.0*EvecRd;


	/* 最良結果の保存用（フラグがonなら） */
	double* best_results=nullptr;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = new double[size];
	}

	bool is_conv = false;
	int bad_counter = 0;

	double zeta=1.0;
	double zeta_old=1.0;
	double eta;
	double nu=1.0;
	/* 反復開始 */
	int It = 0;
	for(It = 0; It < max_ite; It++){
		/* u = C^tr*rd */
		bc_process(size, matL_tr, EvecRd, EvecU);	

		/* ARd = u + C-1 * (Rd-U) */
		tempVec = EvecRd - EvecU;
		fr_process(size, matL, tempVec, EvecARd);
		EvecARd += EvecU;		

		/* (Ard, rd) */
		double Ar_r = EvecARd.dot(EvecRd);
		/* (Ard, Ard) */
		double Ar_Ar =  EvecARd.dot(EvecARd);

		if(It == 0){
			/* ζ(0) */
			zeta = Ar_r / Ar_Ar;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (A*rd(k), y(k)) */
			double Ar_y = EvecARd.dot(EvecY);
			/* ζ(k), η(k) の式の分母 */
			double temp = 1.0 / (nu * Ar_Ar - Ar_y * Ar_y);
			/* ζ(k) */
			zeta = nu * Ar_r * temp;
			/* η(k) */
			eta = -1.0*Ar_y * Ar_r * temp;
		}

		/* ν(k + 1) */
		nu = zeta * Ar_r;

		/* p(k) */
		double temp2 = eta * zeta_old / zeta;
		EvecP = EvecU + temp2 * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY = eta*EvecY + zeta*EvecARd;
		/* rd(k + 1)*/
		EvecRd -= EvecY;

		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecRd.norm();
		/* 収束判定 */
		const double normR = norm_r / normalizer;
		/* フラグがonなら、残差保存 */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* 収束判定 */
		if(normR < conv_cri || norm_r < abs_conv_cri){
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
	}
	delete[] start_posA;
	delete[] end_posA;

	/* 結果代入 */
	for(slv_int i = 0; i < size; i++){
		results[i] = EvecX(i);
	}	
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


/*--------------------*/


/*========================================*/
/*========================================*/
/*//=======================================================
// ● CO-MRTRで解く（本体）
//=======================================================*/
bool MatSolvers::solveSGSMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
	const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *vecB, dcomplex *results, bool init){

	/* 要素確保 */
	Eigen::VectorXcd EvecP = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecRd(size);
	Eigen::VectorXcd EvecU(size);
	Eigen::VectorXcd EvecX(size);	
	Eigen::VectorXcd EvecY(size);
	Eigen::VectorXcd EvecARd(size);

	/* 初期設定 */
	if(init){
		for(slv_int i = 0 ; i < size ; i++){
			results[i] = 0;
			EvecX(i) = 0.0;
		}
	} else{
		for(slv_int i = 0 ; i < size ; i++){
			EvecX(i) = results[i];
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
		dcomplex ap_temp=0;
		for(slv_int j = start_posA[ii] ; j < c_size ; j++){
			ap_temp += val_ptrA[j] * results[col_ptrA[j]];
		}
		/* 初期残差等計算*/
		EvecRd(ii) = vecB[ii] - ap_temp;
		EvecU(ii) = vecB[ii];
	}

	/* 最初から答えだったら何もしない */
	const double first_normR = EvecRd.norm() / normB;
	if(first_normR < conv_cri*0.01){
		delete[] start_posA;
		delete[] end_posA;
		return true;
	}

	/* 収束判定用の補正|B| */
	Eigen::VectorXcd tempVec = matL.matrix*EvecU;
	const double normBd = tempVec.norm();

	/* 前処理rd=C^-1 * r */
	tempVec = EvecRd;
	fr_process(size, matL, tempVec, EvecRd);	
	/* y0 = -r0 */
	EvecY = -1.0*EvecRd;


	/* 最良結果の保存用（フラグがonなら） */
	dcomplex* best_results=nullptr;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = new dcomplex[size];
	}

	bool is_conv = false;
	int bad_counter = 0;

	dcomplex zeta=1.0;
	dcomplex zeta_old=1.0;
	dcomplex eta;
	dcomplex nu=1.0;
	/* 反復開始 */
	int It = 0;
	for(It = 0; It < max_ite; It++){
		/* u = C^tr*rd */
		bc_process(size, matL_tr, EvecRd, EvecU);	

		/* ARd = u + C-1 * (Rd-U) */
		tempVec = EvecRd - EvecU;
		fr_process(size, matL, tempVec, EvecARd);
		EvecARd += EvecU;		

		/* (Ard, rd) */
		//dcomplex Ar_r = EvecARd.dot(EvecRd);
		dcomplex Ar_r = EvecARd.transpose()*EvecRd;
		/* (Ard, Ard) */
		//dcomplex Ar_Ar =  EvecARd.dot(EvecARd);
		dcomplex Ar_Ar =  EvecARd.transpose()*EvecARd;

		if(It == 0){
			/* ζ(0) */
			zeta = Ar_r / Ar_Ar;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (A*rd(k), y(k)) */
			//dcomplex Ar_y = EvecARd.dot(EvecY);
			dcomplex Ar_y = EvecARd.transpose()*EvecY;
			/* ζ(k), η(k) の式の分母 */
			dcomplex temp = 1.0 / (nu * Ar_Ar - Ar_y * Ar_y);
			/* ζ(k) */
			zeta = nu * Ar_r * temp;
			/* η(k) */
			eta = -1.0*Ar_y * Ar_r * temp;
		}

		/* ν(k + 1) */
		nu = zeta * Ar_r;

		/* p(k) */
		dcomplex temp2 = eta * zeta_old / zeta;
		EvecP = EvecU + temp2 * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY = eta*EvecY + zeta*EvecARd;
		/* rd(k + 1)*/
		EvecRd -= EvecY;

		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecRd.norm();
		/* 収束判定 */
		const double normR = norm_r / normBd;
		/* フラグがonなら、残差保存 */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* 収束判定 */
		if(normR < conv_cri){
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
	}
	delete[] start_posA;
	delete[] end_posA;

	/* 結果代入 */
	for(slv_int i = 0; i < size; i++){
		results[i] = EvecX(i);
	}	
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

