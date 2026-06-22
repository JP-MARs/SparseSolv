/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
* @file MatSolverBase.cpp
* @brief Linear matrix solvers for sparse matrices implementation of base class
*/
#include <SparseSolv/MatSolverBase.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/*
//=======================================================
// ■ スパース行列用ソルバの基底クラス
//=======================================================*/

/*//=======================================================
// ● 設定コンストラクタ
//=======================================================*/
/** 
 * @brief Solver base class constructor (set default settings)
*/
MatSolverBase::MatSolverBase(){
	is_init = true;
	is_diag_scale = false;
	is_save_best = false;
	is_save_residual_log = false;
	residual_log.clear();
	/* 発散判定：
	0=最大反復までやる
	1：最良×bad_div_valより大きい値がbad_div_countだけ続いたら発散として終わる */
	diverge_judge_type = 0;
	bad_div_val = 1000.0;
	bad_div_count_thres = 1000;
	/* 収束判定の正規化タイプ(0:右辺、1:初期残差、2:外部指定) */
	conv_normalize_type = 0;
	conv_normalize_const = 1.0;
}

/*//=======================================================
// ● ログ取得
//=======================================================*/
/** 
 * @brief Getter of residual log
*/
void MatSolverBase::getResidualLog(std::vector<double>& log){
	log.clear();
	const size_t the_size = residual_log.size();
	log.resize(the_size);
	if(is_save_residual_log){
		for(size_t i = 0; i < the_size; i++){
			log[i] = residual_log[i];
		}
	}
}

/*//=======================================================
// ● 収束判定正規化定数の計算
//=======================================================*/
/** 
 * @brief set normalizer for convergence check
*/
double MatSolverBase::setNormalizer(const double normB, const double pure_normR){
	/* 収束判定の正規化タイプ(0:右辺、1:初期残差、2:外部指定) */
	/*
	* 0: normalize by norm of rhs
	* 1: normalize by initial residual norm
	* 2: normalize by user specified value
	*/
	if(conv_normalize_type == 0){
		return normB;
	}else if(conv_normalize_type == 1){
		return pure_normR;
	}else if(conv_normalize_type == 2){
		return conv_normalize_const;
	}
	return normB;
}


/*--------------------------------------------------------------------------*/


/*//=======================================================
// ● 対角スケーリング処理時の行列変換
//=======================================================*/
/** 
 * @brief Trans. MatA for diagonal scaling
*/
SparseMat MatSolverBase::makeDiagMatTrans(const SparseMat& matA, const Eigen::VectorXd& matDiagD){
	SparseMat out = matA;

    std::vector<slv_int> start_pos;
    std::vector<slv_int> end_pos;
    out.getCols(start_pos, end_pos);

    const slv_int* col_ptr = out.getColPtr();
    double* val_ptr = out.getValuePtr();

    const slv_int size = out.getSize();
    for (slv_int i = 0; i < size; ++i) {
        const double di = matDiagD(i);
        for (slv_int k = start_pos[i]; k < end_pos[i]; ++k) {
            const slv_int j = col_ptr[k];
            val_ptr[k] *= di * matDiagD(j);
        }
    }
    return out;
}
/** 
 * @brief Trans. MatA for diagonal scaling (complex)
*/
SparseMatC MatSolverBase::makeDiagMatTrans(const SparseMatC& matA, const Eigen::VectorXcd& matDiagD){
	SparseMatC out = matA;

    std::vector<slv_int> start_pos;
    std::vector<slv_int> end_pos;
    out.getCols(start_pos, end_pos);

    const slv_int* col_ptr = out.getColPtr();
    dcomplex* val_ptr = out.getValuePtr();

    const slv_int size = out.getSize();
    for (slv_int i = 0; i < size; ++i) {
        const dcomplex di = matDiagD(i);
        for (slv_int k = start_pos[i]; k < end_pos[i]; ++k) {
            const slv_int j = col_ptr[k];
            val_ptr[k] *= di * matDiagD(j);
        }
    }
    return out;
}



/*//=======================================================
// ● 疎行列積用の対角情報の事前抜き取り
//=======================================================*/
/** 
 * @brief matrix info extraction for sparse matrix-vector product
*/
void MatSolverBase::makeMatrixInfo(const SparseMat& matA, struct MatrixInfoD& mat_info){
	const slv_int size = matA.getSize();

	mat_info.start_pos.resize(size);
	mat_info.end_pos.resize(size);
	matA.getCols(mat_info.start_pos, mat_info.end_pos);
	mat_info.col_ptr = const_cast<slv_int*>(matA.getColPtr());
	mat_info.val_ptr = const_cast<double*>(matA.getValuePtr());
}
void MatSolverBase::makeMatrixInfo(const SparseMatC& matA, struct MatrixInfoC& mat_info){
	const slv_int size = matA.getSize();

	mat_info.start_pos.resize(size);
	mat_info.end_pos.resize(size);
	matA.getCols(mat_info.start_pos, mat_info.end_pos);
	mat_info.col_ptr = const_cast<slv_int*>(matA.getColPtr());
	mat_info.val_ptr = const_cast<dcomplex*>(matA.getValuePtr());
}


/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
/** 
 * @brief for IC decomposition forward and backward substitution
*/
void MatSolverBase::IC_frbc_process(const slv_int size0, const struct MatrixInfoD& mat_infoL, const struct MatrixInfoD& mat_infoL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	Eigen::VectorXd vec_temp = Eigen::VectorXd(size);
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = mat_infoL.end_pos[i] - 1;
		for (slv_int jj = mat_infoL.start_pos[i]; jj < c_size; jj++) {
			slv_int j = mat_infoL.col_ptr[jj];
			s -= mat_infoL.val_ptr[jj] * vec_temp(j);
		}
		vec_temp(i) = s / mat_infoL.val_ptr[c_size];
	}
	for (slv_int i = size; i-- > 0;) {	
		double s = 0;
		const slv_int c_size = mat_infoL_tr.end_pos[i];
		for (slv_int j = mat_infoL_tr.start_pos[i]+1 ; j < c_size; j++) {
			slv_int the_row = mat_infoL_tr.col_ptr[j];
			s -= mat_infoL_tr.val_ptr[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s + vec_temp(i);
	}
}


/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
void MatSolverBase::IC_frbc_process(const slv_int size0, const struct MatrixInfoC& mat_infoL, const struct MatrixInfoC& mat_infoL_tr, const dcomplex *diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;

	/* 第一方程式計算 */
	Eigen::VectorXcd vec_temp = Eigen::VectorXcd(size);
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = mat_infoL.end_pos[i] - 1;
		for (slv_int jj = mat_infoL.start_pos[i]; jj < c_size; jj++) {
			slv_int j = mat_infoL.col_ptr[jj];
			s -= mat_infoL.val_ptr[jj] * vec_temp(j);
		}
		vec_temp(i) = s / mat_infoL.val_ptr[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
		dcomplex s = 0;
		const slv_int c_size = mat_infoL_tr.end_pos[i];
		for (slv_int j = mat_infoL_tr.start_pos[i]+1 ; j < c_size; j++) {
			slv_int the_row = mat_infoL_tr.col_ptr[j];
			s -= mat_infoL_tr.val_ptr[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s + vec_temp(i);
	}
}


/*--------------------------------------------------------------------------*/


/*//=======================================================
// ● 前進代入
//=======================================================*/
/** 
 * @brief forward substitution
*/
void MatSolverBase::fr_process(const slv_int size0, const struct MatrixInfoD& mat_infoL, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;

	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = mat_infoL.end_pos[i] - 1;
		for (slv_int jj = mat_infoL.start_pos[i]; jj < c_size; jj++) {
			slv_int j = mat_infoL.col_ptr[jj];
			s -= mat_infoL.val_ptr[jj] * vec(j);
		}
		vec(i) = s / mat_infoL.val_ptr[c_size];
	}
}



/*//=======================================================
// ● 前進代入
//=======================================================*/
void MatSolverBase::fr_process(const slv_int size0, const struct MatrixInfoC& mat_infoL, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = mat_infoL.end_pos[i] - 1;
		for (slv_int jj = mat_infoL.start_pos[i]; jj < c_size; jj++) {
			slv_int j = mat_infoL.col_ptr[jj];
			s -= mat_infoL.val_ptr[jj] * vec(j);
		}
		vec(i) = s / mat_infoL.val_ptr[c_size];
	}
}


/*----------------------------------------------*/



/*//=======================================================
// ● 後退代入
//=======================================================*/
/** 
 * @brief backward substitution
*/
void MatSolverBase::bc_process(const slv_int size0, const struct MatrixInfoD& mat_infoL_tr, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;

	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = EvecR(i);
		const slv_int c_size = mat_infoL_tr.end_pos[i];
		for (slv_int j = mat_infoL_tr.start_pos[i]+1 ; j < c_size; j++) {
			slv_int the_row = mat_infoL_tr.col_ptr[j];
			s -= mat_infoL_tr.val_ptr[j] * vec(the_row);
		}
		vec(i) = s  / mat_infoL_tr.val_ptr[mat_infoL_tr.start_pos[i]];
	}
}


/*//=======================================================
// ● 後退代入
//=======================================================*/
void MatSolverBase::bc_process(const slv_int size0, const struct MatrixInfoC& mat_infoL_tr, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = EvecR(i);
		const slv_int c_size = mat_infoL_tr.end_pos[i];
		for (slv_int j = mat_infoL_tr.start_pos[i]+1 ; j < c_size; j++) {
			slv_int the_row = mat_infoL_tr.col_ptr[j];
			s -= mat_infoL_tr.val_ptr[j] * vec(the_row);
		}
		vec(i) = s  / mat_infoL_tr.val_ptr[mat_infoL_tr.start_pos[i]];
	}
}


/*----------------------------------------------*/
/*----------------------------------------------*/
/*----------------------------------------------*/
/*----------------------------------------------*/
/*----------------------------------------------*/



/*//=======================================================
// ● 加速係数の自動決定
//=======================================================*/
/** 
 * @brief automatic determination of acceleration coefficient for IC decomposition (make diagonal positive)
*/
void MatSolverBase::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMat& matA, double* diagD, SparseMat& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* 対角が正になるまで実施 */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMat matL_temp = SprsPreconditioner::IC_decomp(matA, diagD, accela_val);
		bool ok = true;
		for(slv_int i = 0; i < size0; i++){
			ok &= (diagD[i] > 0);
		}
		if(ok){
			matL = std::move(matL_temp);
			return;
		}
		accela_val += 0.01;
	}
	/* 失敗時は適当に返すか */
	SparseMat matL_temp0 = SprsPreconditioner::IC_decomp(matA, diagD, 1.8);
	matL = std::move(matL_temp0);
}

/*//=======================================================
// ● 加速係数の自動決定(複素)
//=======================================================*/
void MatSolverBase::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMatC& matA, dcomplex* diagD, SparseMatC& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* 対角が正になるまで実施 */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMatC matL_temp = SprsPreconditioner::IC_decomp(matA, diagD, accela_val); ;//matA.IC_decomp(diagD, accela_val);
		bool ok = true;
		for(slv_int i = 0; i < size0; i++){
			ok &= (diagD[i].real() > 0);
		}
		if(ok){
			matL = std::move(matL_temp);
			return;
		}
		accela_val += 0.01;
	}
	/* 失敗時は適当に返すか */
	SparseMatC matL_temp0 = SprsPreconditioner::IC_decomp(matA, diagD, 1.8);
	matL = std::move(matL_temp0);
}

/* end of namespace */
};
};

