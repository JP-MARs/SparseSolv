/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file MatSolverMRTR.cpp
 * @brief MRTR Linear matrix solvers for sparse matrices implemantation using Eigenstat's preconditioner
 * このファイルは疎行列のMRTR法ソルバの実装を提供します(Eigenstat前処理タイプ)
 */


#include <SparseSolv/MatSolverMRTR.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{



/*

//=======================================================
//=======================================================
//=======================================================
MRTRソルバ
//=======================================================
//=======================================================

*/

/*//=======================================================
// ● SGS+MRTRで解く・呼び出し元
//=======================================================*/
/** 
 * @brief MRTR solver - SGS preconditioner version
*/
bool MatSolverMRTR::solve_SGS_main(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){
	const slv_int size = matA.getSize();	

	/* 対角スケーリング */
	Eigen::VectorXd vecB2;
	Eigen::VectorXd matDiagD = SprsPreconditioner::DiagScaling(matA, vecB2, vecB);
	SparseMat matDAD = makeDiagMatTrans(matA, matDiagD);

	double normB2 = vecB2.norm();

	/* DADの下三角を取得 */
	SparseMat matL = SprsTrans::MatLower(matDAD);
	/* DADの上三角を取得 */
	SparseMat matL_tr = SprsTrans::MatUpper(matDAD);

	/* 解く */
	Eigen::VectorXd results_tmp = std::move(results);
	if(!is_init){
		results_tmp.array() /= matDiagD.array();
	}
	bool bl = solve_SGS_MRTR_main(size, conv_cri, max_ite, normB2, matDAD.mat(), matL, matL_tr, vecB2, results_tmp);

	/* 元に戻す */
	results = std::move(results_tmp);
	results.array() *= matDiagD.array();

	return bl;	
}



/*========================================*/
/*========================================*/
/*//=======================================================
// ● SGS+MRTRで解く（本体）
//=======================================================*/
/** 
 * @brief MRTR solver (main solver body)
*/
bool MatSolverMRTR::solve_SGS_MRTR_main(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& matA, const SparseMat& matL, const SparseMat& matL_tr, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){



	/* 要素確保 */
	Eigen::VectorXd EvecP = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd EvecRd(size);
	Eigen::VectorXd EvecU(size);
	Eigen::VectorXd EvecX = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd EvecY(size);
	Eigen::VectorXd EvecARd(size);

	/* 初期設定 */
	if(is_init){
		results.setZero();
	}else{
		EvecX = results;
	}

	if(is_save_residual_log){
		residual_log.clear();
	}
	/* 最良結果の保存用（フラグがonなら） */
	Eigen::VectorXd best_results;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = Eigen::VectorXd(size);
	}


	/* 初期化 */
	EvecRd = vecB;
	if(is_init){
		EvecRd -= matA*results;
	}
	//EvecU = vecB;
	Eigen::VectorXd tempVec = EvecRd;


	/* 前処理用の行列情報の事前取得 */
	MatrixInfoD mat_infoL;
	MatrixInfoD mat_infoL_tr;
	makeMatrixInfo(matL, mat_infoL);
	makeMatrixInfo(matL_tr, mat_infoL_tr);



	/* 収束判定用の補正|B| */
	/* 前処理rd=C^-1 * r */
	//Eigen::VectorXd tempVec;// = matL.mat()* EvecU;
	fr_process(size, mat_infoL, tempVec, EvecRd);	

	/*------------------------------------------------------------*/
	/* 収束判定系のセット */
	const double pure_normR = EvecRd.norm();
	const double rhs_eps = 1.0e-30;
	/* b = 0 系は絶対残差だけで判定する */
	if (normB <= rhs_eps) {
		if (pure_normR <= small_abs_conv_val) {
			if (is_save_residual_log) {
				residual_log.clear();
				residual_log.push_back(0.0);
			}
			return true;
		}
	}
	const double abs_conv_cri = (normB*conv_cri*0.9 < small_abs_conv_val ? small_abs_conv_val : normB*conv_cri*0.9);
	double first_normR = pure_normR / normB;
	/* 最初から答えだったら何もしない */
	if(first_normR < conv_cri*0.5 || pure_normR < abs_conv_cri*0.5){
		return true;
	}
	/* 残差正規化方法をセット */
	const double normalizer = setNormalizer(normB, pure_normR);

	if(is_save_residual_log){
		residual_log.push_back(pure_normR/normalizer);
	}
	/*------------------------------------------------------------*/


	/* 前処理rd=C^-1 * r */
	//tempVec = EvecRd;
	//fr_process(size, mat_infoL, tempVec, EvecRd);	
	
	/* y0 = -r0 */
	EvecY = -1.0*EvecRd;


	/* 初期値 */
	double zeta=1.0;
	double zeta_old=1.0;
	double eta;
	double nu=1.0;


	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* u = C^tr*rd */
		bc_process(size, mat_infoL_tr, EvecRd, EvecU);	

		/* ARd = u + C-1 * (Rd-U) */
		tempVec = EvecRd - EvecU;
		fr_process(size, mat_infoL, tempVec, EvecARd);	
		EvecARd += EvecU;		

		/* (Ard, rd) */
		double Ar_r = dot_no_conj(EvecARd, EvecRd);// EvecARd.dot(EvecRd);
		/* (Ard, Ard) */
		double Ar_Ar = dot_no_conj(EvecARd, EvecARd);// EvecARd.dot(EvecARd);

		if(It == 0){
			/* ζ(0) */
			zeta = Ar_r / Ar_Ar;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (A*rd(k), y(k)) */
			double Ar_y = dot_no_conj(EvecARd, EvecY);//EvecARd.dot(EvecY);
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
		EvecP *= temp2;
		EvecP += EvecU;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY *= eta;
		EvecY += zeta*EvecARd;
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
				best_results = EvecX;
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
	/* 結果代入 */
	results = EvecX;

	if(!is_conv){
		std::cout << "not Convergence!!! " << std::endl;
		/* 最良値を代入(フラグがonなら) */
		if(is_save_best){
			results = best_results;
		}
	}

	return is_conv ;
}



/*
//=======================================================
//=======================================================
//=======================================================
 複素　ICCOCGソルバ
//=======================================================
//=======================================================
*/


/*//=======================================================
// ● SGS+MRTRで解く・呼び出し元
//=======================================================*/
/** 
 * @brief MRTR solver - SGS preconditioner version
*/
bool MatSolverMRTR::solve_SGS_main(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){
	const slv_int size = matA.getSize();	

	/* 対角スケーリング */
	Eigen::VectorXcd vecB2;
	Eigen::VectorXcd matDiagD = SprsPreconditioner::DiagScaling(matA, vecB2, vecB);
	SparseMatC matDAD = makeDiagMatTrans(matA, matDiagD);

	double normB2 = vecB2.norm();

	/* DADの下三角を取得 */
	SparseMatC matL = SprsTrans::MatLower(matDAD);
	/* DADの上三角を取得 */
	SparseMatC matL_tr = SprsTrans::MatUpper(matDAD);

	/* 解く */
	Eigen::VectorXcd results_tmp = std::move(results);
	if(!is_init){
		results_tmp.array() /= matDiagD.array();
	}
	bool bl = solve_SGS_MRTR_main(size, conv_cri, max_ite, normB2, matDAD.mat(), matL, matL_tr, vecB2, results_tmp);

	/* 元に戻す */
	results = std::move(results_tmp);
	results.array() *= matDiagD.array();

	return bl;	
}


/*========================================*/
/*========================================*/
/*//=======================================================
// ● SGS+MRTRで解く（本体）
//=======================================================*/
/** 
 * @brief SGS+MRTR solver (main solver body)
*/
bool MatSolverMRTR::solve_SGS_MRTR_main(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
			const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& matA, const SparseMatC& matL, const SparseMatC& matL_tr, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){

	/* 要素確保 */
	Eigen::VectorXcd EvecP = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecRd(size);
	Eigen::VectorXcd EvecU(size);
	Eigen::VectorXcd EvecX = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecY(size);
	Eigen::VectorXcd EvecARd(size);


	/* 初期設定 */
	if(is_init){
		results.setZero();
	} else{
		EvecX = results;
	}

	if(is_save_residual_log){
		residual_log.clear();
	}
	/* 最良結果の保存用（フラグがonなら） */
	Eigen::VectorXcd best_results;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = Eigen::VectorXcd(size);
	}

	/* 初期化 */
	EvecRd = vecB;
	if(is_init){
		EvecRd -= matA * results;
	}
	//EvecU = vecB;
	Eigen::VectorXcd tempVec = EvecRd;//matL.mat()*EvecU;

	/* 前処理用の行列情報の事前取得 */
	MatrixInfoC mat_infoL;
	MatrixInfoC mat_infoL_tr;
	makeMatrixInfo(matL, mat_infoL);
	makeMatrixInfo(matL_tr, mat_infoL_tr);


	//Eigen::VectorXcd tempVec = matL.mat()*EvecU;
	//const double normBd = tempVec.norm();
	/* 収束判定用の補正|B| */
	/* 前処理rd=C^-1 * r */
	//Eigen::VectorXd tempVec;// = matL.mat()* EvecU;
	fr_process(size, mat_infoL, tempVec, EvecRd);	

	/*------------------------------------------------------------*/
	/* 収束判定系のセット */
	const double pure_normR = EvecRd.norm();
	const double rhs_eps = 1.0e-30;
	/* b = 0 系は絶対残差だけで判定する */
	if (normB <= rhs_eps) {
		if (pure_normR <= small_abs_conv_val) {
			if (is_save_residual_log) {
				residual_log.clear();
				residual_log.push_back(0.0);
			}
			return true;
		}
	}
	const double abs_conv_cri = (normB*conv_cri*0.9 < small_abs_conv_val ? small_abs_conv_val : normB*conv_cri*0.9);
	double first_normR = pure_normR / normB;
	/* 最初から答えだったら何もしない */
	if(first_normR < conv_cri*0.5 || pure_normR < abs_conv_cri*0.5){
		return true;
	}
	/* 残差正規化方法をセット */
	const double normalizer = setNormalizer(normB, pure_normR);

	if(is_save_residual_log){
		residual_log.push_back(pure_normR/normalizer);
	}
	/*------------------------------------------------------------*/

	/* 前処理rd=C^-1 * r */
	//tempVec = EvecRd;
	//fr_process(size, mat_infoL, tempVec, EvecRd);
	/* y0 = -r0 */
	EvecY = -1.0*EvecRd;


	/* 初期値 */
	dcomplex zeta=1.0;
	dcomplex zeta_old=1.0;
	dcomplex eta;
	dcomplex nu=1.0;

	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* u = C^tr*rd */
		bc_process(size, mat_infoL_tr, EvecRd, EvecU);

		/* ARd = u + C-1 * (Rd-U) */
		tempVec = EvecRd - EvecU;
		fr_process(size, mat_infoL, tempVec, EvecARd);
		EvecARd += EvecU;		

		/* (Ard, rd) */
		//dcomplex Ar_r = EvecARd.dot(EvecRd);
		dcomplex Ar_r = dot_no_conj(EvecARd, EvecRd);// EvecARd.transpose()*EvecRd;
		/* (Ard, Ard) */
		//dcomplex Ar_Ar =  EvecARd.dot(EvecARd);
		dcomplex Ar_Ar =  dot_no_conj(EvecARd, EvecARd);// EvecARd.transpose()*EvecARd;

		if(It == 0){
			/* ζ(0) */
			zeta = Ar_r / Ar_Ar;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (A*rd(k), y(k)) */
			//dcomplex Ar_y = EvecARd.dot(EvecY);
			dcomplex Ar_y = dot_no_conj(EvecARd, EvecY);// EvecARd.transpose()*EvecY;
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
		EvecP *= temp2;
		EvecP += EvecU;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY *= eta;
		EvecY += zeta*EvecARd;
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
				best_results = EvecX;
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
	/* 結果代入 */
	results = EvecX;

	if(!is_conv){
		std::cout << "not Convergence!!! " << std::endl;
		/* 最良値を代入(フラグがonなら) */
		if(is_save_best){
			results = best_results;
		}
	}

	return is_conv ;
}



/* end of namespace */
};
};

