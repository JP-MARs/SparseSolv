/**
* @file MatSolvers_base.cpp
* @brief define basic processes for matrix solvers
*/
#include <SparseSolve/MatSolvers.hpp>
#include <SparseSolve/SparseMatOperators.hpp>

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
ソルバ基本処理
//=======================================================
//=======================================================
*/

/*//=======================================================
// ● 設定コンストラクタ
//=======================================================*/
MatSolvers::MatSolvers(){
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
}

/*//=======================================================
// ● ログ取得
//=======================================================*/
void MatSolvers::getResidualLog(std::vector<double>& log){
	log.clear();
	const size_t the_size = residual_log.size();
	log.resize(the_size);
	if(is_save_residual_log){
		for(size_t i = 0; i < the_size; i++){
			log[i] = residual_log[i];
		}
	}
}

/*--------------------------------------------------------------------------*/


/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const double *vecR, double *vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	/* 第一方程式計算 */
	Eigen::VectorXd vec_temp = Eigen::VectorXd(size);
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] = s + vec_temp(i);
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#ifdef AAAAAAAAAAAAA
	const slv_int size = size0;
	//double s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}

/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	Eigen::VectorXd vec_temp = Eigen::VectorXd(size);
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	for (slv_int i = size; i-- > 0;) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s + vec_temp(i);
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;

#ifdef AAAAAAAAAAAAAAAAAA
	const slv_int size = size0;
	//double s;
	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}


/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	Eigen::VectorXcd vec_temp = Eigen::VectorXcd(size);
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] = s + vec_temp(i);
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#ifdef AAAAAAAAAAAAAAAAAA
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}

/*//=======================================================
// ● IC 前進後退代入
//=======================================================*/
void MatSolvers::IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	Eigen::VectorXcd vec_temp = Eigen::VectorXcd(size);
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec_temp(j);
		}
		vec_temp(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) = s + vec_temp(i);
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;

#ifdef AAAAAAAAAAAAAAA
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
#endif
}


/*--------------------------------------------------------------------------*/



/*//=======================================================
// ● 前進代入
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseD& matL, const double *vecR, double *vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}

/*//=======================================================
// ● 前進代入
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseD& matL, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;
	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		double s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}


/*//=======================================================
// ● 前進代入
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseC& matL, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}

/*//=======================================================
// ● 前進代入
//=======================================================*/
void MatSolvers::fr_process(const slv_int size0, const SparseMatBaseC& matL, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	delete[] start_posL1;
	delete[] end_posL1;
}


/*----------------------------------------------*/



/*//=======================================================
// ● 後退代入
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const double *vecR, double *vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();

	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = vecR[i];
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		vec[i] = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// ● 後退代入
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	//double s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		double s = EvecR(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		vec(i) = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}


/*//=======================================================
// ● 後退代入
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const dcomplex *vecR, dcomplex *vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = vecR[i];
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec[the_row];
		}
		vec[i] = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// ● 後退代入
//=======================================================*/
void MatSolvers::bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec){
	const slv_int size = size0;
	//dcomplex s;

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第二方程式計算 */
	for (slv_int i = size; i-- > 0;) {	
	//for (slv_int i = size-1; i >= 0; i--) {	
		dcomplex s = EvecR(i);
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s -= val_ptrL2[j] * vec(the_row);
		}
		vec(i) = s  / val_ptrL2[start_posL2[i]];
	}
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// ● 加速係数の自動決定
//=======================================================*/
void MatSolvers::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMat& matA, double* diagD, SparseMat& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* 対角が正になるまで実施 */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMat matL_temp = matA.IC_decomp(diagD, accela_val);
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
}

/*//=======================================================
// ● 加速係数の自動決定(複素)
//=======================================================*/
void MatSolvers::auto_accel_determine(const slv_int size0, double accel_ini, const SparseMatC& matA, dcomplex* diagD, SparseMatC& matL){
	double accela_val = fabs(accel_ini);
	if(accela_val < 0.9 || accela_val > 1.8){
		accela_val = 1;
	}
	/* 対角が正になるまで実施 */
	for(int kkk = 0; kkk < 80; kkk++){
		SparseMatC matL_temp = matA.IC_decomp(diagD, accela_val);
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
}

/* end of namespace */
};

