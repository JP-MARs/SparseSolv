/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file SprsPreconditioner.cpp
 * @brief Matrix Preconditioner Operators for Sparse Matrix Solvers implementation
 * 
 * このファイルは ソルバで使う前処理の実装です
 */

#include <cfloat>
#include <SparseSolv/SprsPreconditioner.hpp>


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{
/* 疎行列の前処理用の専用名前空間 */
namespace SprsPreconditioner{


/* 疎行列前処理用の内部空間間 */
/**
* @namespace SprsPreconditioner::detai
* @brief Original namespace for internal use of Matrix Preconditioner
*/
namespace detail {
	
/*//=======================================================
// ● 不完全コレスキー分解 サブルーチン
//=======================================================*/
/** 
 * @brief IC decomposition operator with accearation factor subrutin
 * @param accela acceleration factor
*/
template<typename DType1>
void IC_decomp_sub(SparseMatTMPL<DType1>& mat1, DType1* diagD, const double accela){

	const slv_int the_size = mat1.getSize();
	for (slv_int i = 0; i < the_size; i++) {
		diagD[i] = 0;
	}

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	mat1.getCols(start_pos1, end_pos1);
	auto col_ptr1 = mat1.getColPtr();
	auto val_ptr1 = mat1.getValuePtr();

	/* LとDを求める */
	DType1 s;
	/* LとDを求める */
	for (slv_int i = 0; i < the_size; i++) {
		/* i行目にある非ゼロの数をゲット */
		const slv_int c_size = end_pos1[i];
		/* i行目の非ゼロの回数だけ列ループをまわす */
		for (slv_int jj = start_pos1[i]; jj < c_size; jj++) {
			/* 現在の列番号jを取得 */
			const slv_int j = col_ptr1[jj];
			if(j >= i) break;
			s = val_ptr1[jj];
			/* i行目の列の非ゼロ列をサーチ */
			const slv_int L_ksize = end_pos1[i];
			for(slv_int kk = start_pos1[i] ; kk < L_ksize ; kk++){
				if(col_ptr1[kk] <= j-1){
					const slv_int L_Jsize = end_pos1[j];;//tempMat->column_size[j];
					for(slv_int ll = start_pos1[j] ; ll < L_Jsize ; ll++){
						if(col_ptr1[ll] == col_ptr1[kk]){
							s -= val_ptr1[kk] * val_ptr1[ll] * diagD[col_ptr1[ll]];
							break;
						}else if(col_ptr1[kk] < col_ptr1[ll]){
							break;
						}
					}
				}else{
					break;
				}
			}
			/* 値を[i][j]に代入 */
			val_ptr1[jj] = s;
		}
		const slv_int last = end_pos1[i] - 1;
		//std::cout << "diag " << i << ", " << last << ", "<< col_ptr1[last] << std::endl;
		/* もし対角が無ければエラーで落とす */
		if(col_ptr1[last] != i){
			std::cout << "there is not diagonal element! "<< std::endl;
			exit(1);
		}
		s = val_ptr1[last] * accela;
		for (slv_int kk = start_pos1[i]; kk < last; kk++) {
			if(col_ptr1[kk] <= i-1){
				s -= val_ptr1[kk] * val_ptr1[kk] * diagD[col_ptr1[kk]];
			}else{
				break;
			}
		}
		/* 対角に代入 */
		val_ptr1[last] = s;
		if(abs(s) < 1.0e-12){
			diagD[i] = 1.0 / (DBL_EPSILON);
		}else{
			diagD[i] = 1.0 / s;
		}
	}
}

/* end of subnamespace */
};


/*==================================================================*/
/*==================================================================*/
/* 実装本体*/
/*==================================================================*/
/*==================================================================*/
/*==================================================================*/



/*//=======================================================
// ● 不完全コレスキー分解
//=======================================================*/
/** 
 * @brief IC decomposition operator with accearation factor
 * @param accela acceleration factor
*/
SparseMat IC_decomp(const SparseMat& mat1, double* diagD, const double accela){
	SparseMat mat_lower = SprsTrans::MatLower(mat1);
	detail::IC_decomp_sub<double>(mat_lower, diagD, accela);
	return mat_lower;
}
SparseMatC IC_decomp(const SparseMatC& mat1, dcomplex* diagD, const double accela){
	SparseMatC mat_lower = SprsTrans::MatLower(mat1);
	detail::IC_decomp_sub<dcomplex>(mat_lower, diagD, accela);
	return mat_lower;
}

/*//=======================================================
// ● 対角スケーリング
//=======================================================*/
/** 
 * @brief diagonal scaling operator
*/
Eigen::VectorXd DiagScaling(const SparseMat& mat1, Eigen::VectorXd& trans_vec, const Eigen::VectorXd& ori_vec){

	const slv_int the_size = mat1.getSize();
	trans_vec = Eigen::VectorXd(the_size);

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	mat1.getCols(start_pos1, end_pos1);
	auto col_ptr1 = mat1.getColPtr();
	auto val_ptr1 = mat1.getValuePtr();

	Eigen::VectorXd diagD(the_size);

	for (slv_int i = 0; i < the_size; i++) {
		/* i行目にある非ゼロの数をゲット */
		const slv_int c_size = end_pos1[i];
		/* i行目の非ゼロの回数だけ列ループをまわす */
		for(slv_int jj = start_pos1[i]; jj < c_size; jj++) {
			/* 対角位置なら処理 */
			if(col_ptr1[jj] == i){
				double temp1 = abs(val_ptr1[jj]);
				double temp1b = sqrt(temp1);
				double temp2;
				if (temp1b < 1.0e-12) {
					temp2 = 1.0 / (DBL_EPSILON);
				}else{
					temp2 = 1.0 / temp1b;
				}
				diagD(i) = temp2;
				/* 右辺も更新 */
				trans_vec(i) = ori_vec(i) * temp2;
				break;
			}
		}
	}
	return diagD;
}

/*//=======================================================
// ● 対角スケーリング（複素）
//=======================================================*/
/** 
 * @brief diagonal scaling operator for complex matrix
*/
Eigen::VectorXcd DiagScaling(const SparseMatC& mat1, Eigen::VectorXcd& trans_vec, const Eigen::VectorXcd& ori_vec){

	const slv_int the_size = mat1.getSize();
	trans_vec = Eigen::VectorXcd(the_size);

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	mat1.getCols(start_pos1, end_pos1);
	auto col_ptr1 = mat1.getColPtr();
	auto val_ptr1 = mat1.getValuePtr();

	Eigen::VectorXcd diagD(the_size);

	for (slv_int i = 0; i < the_size; i++) {
		/* i行目にある非ゼロの数をゲット */
		const slv_int c_size = end_pos1[i];
		/* i行目の非ゼロの回数だけ列ループをまわす */
		for(slv_int jj = start_pos1[i]; jj < c_size; jj++) {
			/* 対角位置なら処理 */
			if(col_ptr1[jj] == i){
				double temp1 = std::norm(val_ptr1[jj]);
				double temp1b = sqrt(temp1);
				double temp2;
				if (temp1b < 1.0e-12) {
					temp2 = 1.0 / (DBL_EPSILON);
				}else{
					temp2 = 1.0 / temp1b;
				}
				dcomplex temp3(temp2, 0.0);
				diagD(i) = temp3;
				/* 右辺も更新 */
				trans_vec(i) = ori_vec(i) * temp3;
				break;
			}
		}
	}
	return diagD;
}

/* end of namespace */
};
};
};


