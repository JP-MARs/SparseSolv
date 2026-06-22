/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file MatSolverMRTR.cpp
 * @brief MRTR Linear matrix solvers for sparse matrices implemantation using IC preconditioner
 * このファイルは疎行列のMRTR法ソルバの実装を提供します(IC前処理タイプ)
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
// MRTRソルバ
//=======================================================
//=======================================================

*/


/*//=======================================================
// ● MRTRで解く(右辺ノルム内部計算パターン)
//=======================================================*/
/** 
 * @brief MRTR solver (calculate norm of rhs internally)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const SparseMat& matA, const double *vecB, double *results) {
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
// ● MRTRで解く(右辺ノルム内部計算パターン / Eigen)
//=======================================================*/
/**
 * @brief MRTR solver (calculate norm of rhs internally, Eigen ver)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) {
	double norm = vecB.norm();
	return solve(conv_cri, max_ite, norm, matA, vecB, results);
}

/*//=======================================================
// ● MRTRで解く - 右辺ポインタ、Eigenに変換する
//=======================================================*/
/**
 * @brief MRTR solver (ptr ver >> change Eigen::Vector)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double* vecB, double* results) {
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
// ● MRTRで解く・外部実行本体
//=======================================================*/
/** 
 * @brief MRTR solver (with autmatic accel. factor determine)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){

	/* 前処理がEigenstatの方法の場合、別処理へ飛んで終わる */
	if(precondition_type == 1){
		return solve_SGS_main(conv_cri, max_ite, normB, matA, vecB, results);
	}
	
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
// ● MRTR専用内部処理(対角スケーリングあり)
//=======================================================*/
/** 
 * @brief MRTR solver (main solver body with diag scaling)
*/
bool MatSolverMRTR::solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){
	
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
// ● MRTRで解く（本体）
//=======================================================*/
/** 
 * @brief MRTR solver (main solver body)
*/
bool MatSolverMRTR::solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const double* diagD, const Eigen::SparseMatrix<double, Eigen::RowMajor>& matA, const SparseMat& matL, const SparseMat& matL_tr, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){

	/* 要素確保 */
	Eigen::VectorXd EvecP = Eigen::VectorXd::Zero(size);
	Eigen::VectorXd EvecR(size);
	Eigen::VectorXd EvecU(size);
	Eigen::VectorXd EvecX = Eigen::VectorXd::Zero(size);	
	Eigen::VectorXd EvecY(size);
	Eigen::VectorXd EvecZ(size);


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
	Eigen::VectorXd best_results;
	double best_resi_value = 1.0e+6;
	if(is_save_best){
		best_results = Eigen::VectorXd(size);
	}

	/* 初期化 */
	EvecR = vecB;
	if(is_init){
		EvecR -= matA * results;
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
	/* 前処理u=M^-1 * r */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecU);
	/* y0 = -r0 */
	EvecY = -1.0*EvecR;
	/* 前処理z=M^-1 * y */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecY, EvecZ);


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
		/* v = Au */
		Eigen::VectorXd vecV = matA * EvecU;
		/* w = M^-1 * v */
		Eigen::VectorXd vecW(size);
		IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, vecV, vecW);

		/* (w, r(k)) */
		double w_r = dot_no_conj(vecW, EvecR);// vecW.dot(EvecR);
		/* (v, w) */
		double v_w = dot_no_conj(vecV, vecW);// vecV.dot(vecW);

		if(It == 0){
			/* ζ(0) */
			zeta = w_r / v_w;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (w, y(k)) */
			double w_y = dot_no_conj(vecW, EvecY);// vecW.dot(EvecY);
			/* ζ(k), η(k) の式の分母 */
			const double temp0 = 1.0 / (nu * v_w - w_y * w_y);
			/* ζ(k) */
			zeta = nu * w_r * temp0;
			/* η(k) */
			eta = -1.0*w_y * w_r * temp0;
		}

		/* ν(k + 1) */
		nu = zeta * w_r;
		/* p(k) */
		double temp = eta * zeta_old / zeta;
		EvecP = EvecU;
		EvecP += temp * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY *= eta;
		EvecY += zeta*vecV;
		/* r(k + 1)*/
		EvecR -= EvecY;


		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecR.norm();
		/* 収束判定 */
		const double normR = norm_r / normalizer;;
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
		/* z(k + 1) */
		EvecZ *= eta;
		EvecZ += zeta*vecW;
		/* u(k + 1) */
		EvecU -= EvecZ;
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
// ● MRTRで解く(右辺ノルム内部計算パターン)
//=======================================================*/
/** 
 * @brief ICCOCG solver (calculate norm of rhs internally)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex *vecB, dcomplex* results) {
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
// ● MRTRで解く(右辺ノルム内部計算パターン / Eigen)
//=======================================================*/
/**
 * @brief MRTR solver (calculate norm of rhs internally, Eigen ver)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) {
	double norm = vecB.norm();
	return solve(conv_cri, max_ite, norm, matA, vecB, results);
}

/*//=======================================================
// ● MRTRで解く - 右辺ポインタ、Eigenに変換する
//=======================================================*/
/**
 * @brief MRTR solver (ptr ver >> change Eigen::Vector)
 * @param init Initialize solution vector by zero or not (if false, user must initialize primaly)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) {
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
// ● MRTRで解く・外部実行本体
//=======================================================*/
/** 
 * @brief ICCOCG solver ((with autmatic accel. factor determine)
*/
bool MatSolverMRTR::solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){

	/* 前処理がEigenstatの方法の場合、別処理へ飛んで終わる */
	if(precondition_type == 1){
		return solve_SGS_main(conv_cri, max_ite, normB, matA, vecB, results);
	}


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
// ● MRTRで解く・外部実行本体
//=======================================================*/
/** 
 * @brief ICCOCG solver (main solver body with diag scaling)
*/
bool MatSolverMRTR::solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){
	
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
// ● MRTRで解く（本体）
//=======================================================*/
/** 
 * @brief ICCOCG solver (main solver body)
*/
bool MatSolverMRTR::solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const dcomplex* diagD, const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& matA, const SparseMatC& matL, const SparseMatC& matL_tr, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){

	/* 要素確保 */
	Eigen::VectorXcd EvecP = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecR(size);
	Eigen::VectorXcd EvecU(size);
	Eigen::VectorXcd EvecX = Eigen::VectorXcd::Zero(size);
	Eigen::VectorXcd EvecY(size);
	Eigen::VectorXcd EvecZ(size);


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
	EvecR = vecB;
	if(is_init){
		EvecR -= matA * results;
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
	/* 前処理u=M^-1 * r */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecR, EvecU);
	/* y0 = -r0 */
	EvecY = -1.0*EvecR;
	/* 前処理z=M^-1 * y */
	IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, EvecY, EvecZ);

	/* 初期値 */
	dcomplex zeta;
	dcomplex zeta_old;
	dcomplex eta;
	dcomplex nu=1.0;

	bool is_conv = false;
	int bad_counter = 0;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* v = Au */
		Eigen::VectorXcd vecV = matA * EvecU;
		/* w = M^-1 * v */
		Eigen::VectorXcd vecW(size);
		IC_frbc_process(size, mat_infoL, mat_infoL_tr, diagD, vecV, vecW);

		/* (w, r(k)) */
		//dcomplex w_r = vecW.dot(EvecR);
		dcomplex w_r = dot_no_conj(vecW, EvecR);// vecW.transpose()* EvecR;
		/* (v, w) */
		//dcomplex v_w = vecV.dot(vecW);
		dcomplex v_w = dot_no_conj(vecV, vecW);// vecV.transpose()* vecW;

		if(It == 0){
			/* ζ(0) */
			zeta = w_r / v_w;
			zeta_old = zeta;
			/* η(0) */
			eta = 0.0;
		} else{
			/* (w, y(k)) */
			//dcomplex w_y = vecW.dot(EvecY);
			dcomplex w_y = dot_no_conj(vecW, EvecY);//vecW.transpose()* EvecY;
			/* ζ(k), η(k) の式の分母 */
			const dcomplex temp0 = 1.0 / (nu * v_w - w_y * w_y);
			/* ζ(k) */
			zeta = nu * w_r * temp0;
			/* η(k) */
			eta = -1.0*w_y * w_r * temp0;
		}

		/* ν(k + 1) */
		nu = zeta * w_r;


		/* p(k) */
		dcomplex temp = eta * zeta_old / zeta;
		EvecP = EvecU;
		EvecP += temp * EvecP;
		zeta_old = zeta;

		/* x(k + 1) */
		EvecX += zeta*EvecP;
		/* y(k + 1) */
		EvecY *= eta;
		EvecY += zeta*vecV;
		/* r(k + 1)*/
		EvecR -= EvecY;

		/* (r(k + 1), r(k + 1)) */
		double norm_r = EvecR.norm();
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

		/* z(k + 1) */
		EvecZ *= eta;
		EvecZ += zeta*vecW;
		/* u(k + 1) */
		EvecU -= EvecZ;
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

