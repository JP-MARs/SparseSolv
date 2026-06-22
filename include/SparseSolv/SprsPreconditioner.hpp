/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 
 /**
 * @file SprsPreconditioner.hpp
 * @brief Matrix Preconditioner Operators for Sparse Matrix Solvers
 * 
 * このファイルは ソルバで使う前処理を提供します
 */
#ifndef DEF_SPARSE_MAT_PRECONDITIONER_OPERATOR_DEFINES
#define DEF_SPARSE_MAT_PRECONDITIONER_OPERATOR_DEFINES

#include "SparseMat.hpp"
#include "SprsOperator.hpp"
#include "SprsTrans.hpp"


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/* 疎行列の前処理用の専用名前空間 */
/**
* @namespace SprsPreconditioner
* @brief Original namespace for sparse matrix Matrix Preconditioner Operators
*/
namespace SprsPreconditioner{
	/* 不完全コレスキー分解 */
	SparseMat IC_decomp(const SparseMat& mat1, double* diagD, const double accela);
	SparseMatC IC_decomp(const SparseMatC& mat1, dcomplex* diagD, const double accela);

	/* 対角スケーリング */
	Eigen::VectorXd DiagScaling(const SparseMat& mat1, Eigen::VectorXd& trans_vec, const Eigen::VectorXd& ori_vec);
	Eigen::VectorXcd DiagScaling(const SparseMatC& mat1, Eigen::VectorXcd& trans_vec, const Eigen::VectorXcd& ori_vec);
}


/* end of namespace */
};
};


#endif


