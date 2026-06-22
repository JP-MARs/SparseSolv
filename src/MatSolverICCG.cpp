/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file MatSolverICCG.cpp
 * @brief ICCG Linear matrix solvers for sparse matrices implemantation
 * このファイルは疎行列のICCG法ソルバの実装を提供します
 */


#include <SparseSolv/MatSolverICCG.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{



/*

//=======================================================
//=======================================================
//=======================================================
ICCGソルバ
//=======================================================
//=======================================================

*/


/*//=======================================================
// ● ICCGで解く(右辺ノルム内部計算パターン)
//=======================================================*/
/** 
 * @brief ICCG solver (calculate norm of rhs internally)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const SparseMat& matA, const double *vecB, double *results) {
	const slv_int size = matA.getSize();
	Eigen::VectorXd vecB_eigen(size);
	Eigen::VectorXd results_eigen(size);
	double norm = 0;
	for (slv_int i = 0; i < size; i++) {
		vecB_eigen(i) = vecB[i];
		results_eigen(i) = results[i];
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	bool bl = solve(conv_cri, max_ite, norm, matA, vecB_eigen, results_eigen);
	Eigen::Map<Eigen::VectorXd>(results, size) = results_eigen;
	return bl;
}
/*//=======================================================
// ● ICCGで解く(右辺ノルム内部計算パターン / Eigen)
//=======================================================*/
/**
 * @brief ICCG solver (calculate norm of rhs internally, Eigen ver)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) {
	double norm = vecB.norm();
	return solve(conv_cri, max_ite, norm, matA, vecB, results);
}

/*//=======================================================
// ● ICCGで解く - 右辺ポインタ、Eigenに変換する
//=======================================================*/
/**
 * @brief ICCG solver (ptr ver >> change Eigen::Vector)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double* vecB, double* results) {
	const slv_int size = matA.getSize();
	Eigen::VectorXd vecB_eigen(size);
	Eigen::VectorXd results_eigen(size);
	for (slv_int i = 0; i < size; i++) {
		vecB_eigen(i) = vecB[i];
		results_eigen(i) = results[i];
	}
	bool bl = solve(conv_cri, max_ite, normB, matA, vecB_eigen, results_eigen);
	Eigen::Map<Eigen::VectorXd>(results, size) = results_eigen;
	return bl;
}



/*//=======================================================
// ● ICCGで解く・外部実行本体
//=======================================================*/
/** 
 * @brief ICCG solver (with autmatic accel. factor determine)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){
	
	/* 対角スケーリングあり */
	if(is_diag_scale){
		return solve_diag_scale(conv_cri, max_ite, normB, matA, vecB, results);
	}
	
	/* 通常 */
	const slv_int size = matA.getSize();
	/* コレスキー用スパース行列作成 */
	double* diagD = new double[size];
	SparseMat matL;
	/* 加速係数が負なら自動決定モードへ */
	if(accel_factor < 0.0){
		auto_accel_determine(size, accel_factor, matA, diagD, matL);
	}else{
		matL = SprsPreconditioner::IC_decomp(matA, diagD, accel_factor);
	}
	/* 確定 */
	SparseMat matL_tr = SprsTrans::trans(matL);

	bool bl = solve_main(size, conv_cri, max_ite, normB, diagD, matA.mat(), matL, matL_tr, vecB, results);
	delete[] diagD;
	return bl;
}

/*  */

/*//=======================================================
// ● ICCG専用内部処理(対角スケーリングあり)
//=======================================================*/
/** 
 * @brief ICCG solver (main solver body with diag scaling)
*/
bool MatSolverICCG::solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){
	
	const slv_int size = matA.getSize();
	/* 対角スケーリング */
	Eigen::VectorXd vecB2;
	Eigen::VectorXd matDiagD = SprsPreconditioner::DiagScaling(matA, vecB2, vecB);
	SparseMat matDAD = makeDiagMatTrans(matA, matDiagD);

	double normB2 = vecB2.norm();

	/* コレスキー用スパース行列作成 */
	double* diagD = new double[size];
	SparseMat matL;
	/* 加速係数が負なら自動決定モードへ */
	if(accel_factor < 0.0){
		auto_accel_determine(size, accel_factor, matDAD, diagD, matL);
	}else{
		matL = SprsPreconditioner::IC_decomp(matDAD, diagD, accel_factor);
	}
	/* 確定 */
	SparseMat matL_tr = SprsTrans::trans(matL);

	/* 解く */
	Eigen::VectorXd results_tmp = std::move(results);
	if(!is_init){
		results_tmp.array() /= matDiagD.array();
	}
	bool bl = solve_main(size, conv_cri, max_ite, normB2, diagD, matDAD.mat(), matL, matL_tr, vecB2, results_tmp);
	delete[] diagD;

	/* 元に戻す */
	results = std::move(results_tmp);
	results.array() *= matDiagD.array();
	
	return bl;
}

/*========================================*/
/*========================================*/
/*//=======================================================
// ● ICCGで解く（本体）
//=======================================================*/
/** 
 * @brief ICCG solver (main solver body)
*/
bool MatSolverICCG::solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const double* diagD, const Eigen::SparseMatrix<double, Eigen::RowMajor>& matA, const SparseMat& matL, const SparseMat& matL_tr, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){

	/* 要素確保 */
	double alpha;
	double beta;
	Eigen::VectorXd EvecP(size);
	Eigen::VectorXd EvecR(size);
	Eigen::VectorXd EvecLDV(size);
	Eigen::VectorXd EtempAP(size);


	/* 初期設定 */
	if(is_init){
		results.setZero();
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
	EvecR = vecB;
	if(!is_init){
		EtempAP = matA * results;
		EvecR -= EtempAP;
	}


	/*------------------------------------------------------------*/
	/* 収束判定系のセット */
	const double pure_normR = EvecR.norm();
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

	/* 前処理用の行列情報の事前取得 */
	MatrixInfoD mat_infoL;
	MatrixInfoD mat_infoL_tr;
	makeMatrixInfo(matL, mat_infoL);
	makeMatrixInfo(matL_tr, mat_infoL_tr);
	/* 前処理 */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecP);
	EvecLDV = EvecP;

	/* 初期値 */
	double rur0 = dot_no_conj(EvecR, EvecLDV);


	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* AP計算 */
		EtempAP = matA * EvecP;
		/* α計算 */
		const double temp2 = dot_no_conj(EvecP, EtempAP);
		alpha = rur0 / temp2;

		/* 解ベクトルと残差計算 */
		results += alpha * EvecP;
		EvecR -= (alpha * EtempAP);
		const double pure_norm = EvecR.norm();
		const double normR = pure_norm / normalizer;

		/* フラグがonなら、残差保存 */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* 収束判定 */
		if(normR < conv_cri || pure_norm < abs_conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << ", " << normR <<std::endl;
			is_conv = true;
			break;
		//}else if( It % 10 == 0 ){
		//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		
		if(normR < best_resi_value){
			best_resi_value = normR;
			/* 最良値の更新(フラグがonなら) */		
			if(is_save_best){
				best_results = results;
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
		IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecLDV);
		/* β計算 */
		const double rur1 = dot_no_conj(EvecR, EvecLDV);
		beta = rur1 / rur0;
		rur0 = rur1;
		//beta = EvecR.dot(EvecLDV);
		//beta /= temp;
		/* P計算 */
		EvecP *= beta;
		EvecP += EvecLDV;
	}
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
// ● ICCGで解く(右辺ノルム内部計算パターン)
//=======================================================*/
/** 
 * @brief ICCOCG solver (calculate norm of rhs internally)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex *vecB, dcomplex* results) {
	const slv_int size = matA.getSize();
	Eigen::VectorXcd vecB_eigen(size);
	Eigen::VectorXcd results_eigen(size);
	for (slv_int i = 0; i < size; i++) {
		vecB_eigen(i) = vecB[i];
		results_eigen(i) = results[i];
	}
	double normB2 = vecB_eigen.norm();
	bool bl = solve(conv_cri, max_ite, normB2, matA, vecB_eigen, results_eigen);
	Eigen::Map<Eigen::VectorXcd>(results, size) = results_eigen;
	return bl;
}
/*//=======================================================
// ● ICCGで解く(右辺ノルム内部計算パターン / Eigen)
//=======================================================*/
/**
 * @brief ICCG solver (calculate norm of rhs internally, Eigen ver)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) {
	double norm = vecB.norm();
	return solve(conv_cri, max_ite, norm, matA, vecB, results);
}

/*//=======================================================
// ● ICCGで解く - 右辺ポインタ、Eigenに変換する
//=======================================================*/
/**
 * @brief ICCG solver (ptr ver >> change Eigen::Vector)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) {
	const slv_int size = matA.getSize();
	Eigen::VectorXcd vecB_eigen(size);
	Eigen::VectorXcd results_eigen(size);
	for (slv_int i = 0; i < size; i++) {
		vecB_eigen(i) = vecB[i];
		results_eigen(i) = results[i];
	}
	bool bl = solve(conv_cri, max_ite, normB, matA, vecB_eigen, results_eigen);
	Eigen::Map<Eigen::VectorXcd>(results, size) = results_eigen;
	return bl;
}



/*//=======================================================
// ● ICCGで解く・外部実行本体
//=======================================================*/
/** 
 * @brief ICCOCG solver ((with autmatic accel. factor determine)
*/
bool MatSolverICCG::solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){
	/* 対角スケーリングあり */
	if(is_diag_scale){
		return solve_diag_scale(conv_cri, max_ite, normB, matA, vecB, results);
	}

	/* 通常 */
	const slv_int size = matA.getSize();
	/* コレスキー用スパース行列作成 */
	dcomplex* diagD = new dcomplex[size];
	/* 加速係数が負なら自動決定モードへ */
	SparseMatC matL;
	if(accel_factor < 0.0){
		auto_accel_determine(size, accel_factor, matA, diagD, matL);
	}else{
		matL = SprsPreconditioner::IC_decomp(matA, diagD, accel_factor);
	}
	/* 確定 */
	SparseMatC matL_tr = SprsTrans::trans(matL);

	bool bl= solve_main(size, conv_cri, max_ite, normB, diagD, matA.mat(), matL, matL_tr, vecB, results);
	delete[] diagD;
	return bl;
}

/*//=======================================================
// ● ICCGで解く・外部実行本体
//=======================================================*/
/** 
 * @brief ICCOCG solver (main solver body with diag scaling)
*/
bool MatSolverICCG::solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){
	
	const slv_int size = matA.getSize();
	/* 対角スケーリング */
	Eigen::VectorXcd vecB2;
	Eigen::VectorXcd matDiagD = SprsPreconditioner::DiagScaling(matA, vecB2, vecB);
	SparseMatC matDAD = makeDiagMatTrans(matA, matDiagD);

	double normB2 = vecB2.norm();

	/* コレスキー用スパース行列作成 */
	dcomplex* diagD = new dcomplex[size];
	SparseMatC matL;
	/* 加速係数が負なら自動決定モードへ */
	if(accel_factor < 0.0){
		auto_accel_determine(size, accel_factor, matDAD, diagD, matL);
	}else{
		matL = SprsPreconditioner::IC_decomp(matDAD, diagD, accel_factor);
	}
	/* 確定 */
	SparseMatC matL_tr = SprsTrans::trans(matL);

	/* 解く */
	Eigen::VectorXcd results_tmp = std::move(results);
	if(!is_init){
		results_tmp.array() /= matDiagD.array();
	}
	bool bl = solve_main(size, conv_cri, max_ite, normB2, diagD, matDAD.mat(), matL, matL_tr, vecB2, results_tmp);
	delete[] diagD;

	/* 元に戻す */
	results = std::move(results_tmp);
	results.array() *= matDiagD.array();
	
	return bl;
}

/*========================================*/
/*========================================*/
/*//=======================================================
// ● ICCGで解く（本体）
//=======================================================*/
/** 
 * @brief ICCOCG solver (main solver body)
*/
bool MatSolverICCG::solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const dcomplex* diagD, const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& matA, const SparseMatC& matL, const SparseMatC& matL_tr, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){

	/* 要素確保 */
	dcomplex alpha;
	dcomplex beta;
	Eigen::VectorXcd EvecP(size);
	Eigen::VectorXcd EvecR(size);
	Eigen::VectorXcd EvecLDV(size);
	Eigen::VectorXcd EtempAP(size);


	/* 初期設定 */
	if(is_init){
		results.setZero();
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
	EvecR = vecB;
	if(!is_init){
		EtempAP = matA * results;
		EvecR -= EtempAP;
	}

	/*------------------------------------------------------------*/
	/* 収束判定系のセット */
	const double pure_normR = EvecR.norm();
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

	/* 前処理用の行列情報の事前取得 */
	MatrixInfoC mat_infoL;
	MatrixInfoC mat_infoL_tr;
	makeMatrixInfo(matL, mat_infoL);
	makeMatrixInfo(matL_tr, mat_infoL_tr);
	/* 前処理 */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecP);
	EvecLDV = EvecP;

	/* 初期値 */
	dcomplex rur0 = dot_no_conj(EvecR, EvecLDV);

	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* AP計算 */
		EtempAP = matA * EvecP;
		/* α計算 */
		const dcomplex temp2 = dot_no_conj(EvecP, EtempAP);//EvecP.transpose()*EtempAP;
		alpha = rur0 / temp2;

		/* 解ベクトルと残差計算 */
		results += alpha * EvecP;
		EvecR -= alpha * EtempAP;
		const double pure_norm = EvecR.norm();
		const double normR = pure_norm / normalizer;


		/* フラグがonなら、残差保存 */
		if(is_save_residual_log){
			residual_log.push_back(normR);
		}
		/* 収束判定 */
		if(normR < conv_cri || pure_norm < abs_conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << std::endl;
			is_conv = true;
			break;
			//}else if( It % 100 == 0 ){
			//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		if(normR < best_resi_value){
			best_resi_value = normR;
			/* 最良値の更新(フラグがonなら) */		
			if(is_save_best){
				best_results = results;
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
		IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecLDV);
		/* β計算 */
		//const dcomplex rur1 = EvecR.dot(EvecLDV);
		const dcomplex rur1 = dot_no_conj(EvecR, EvecLDV);//EvecR.transpose()*EvecLDV;
		beta = rur1 / rur0;
		rur0 = rur1;
		//beta = EvecR.dot(EvecLDV);
		//beta /= temp;
		/* P計算 */
		EvecP *= beta;
		EvecP += EvecLDV;
	}
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

