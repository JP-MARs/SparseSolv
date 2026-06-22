/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 /**
 * @file MatSolversMRTR.hpp
 * @brief MRTR Linear matrix solvers for sparse matrices
 * このファイルは疎行列のMRTR法ソルバを提供します
 */
#ifndef DEF_MAT_SOLVER_MRTR_MYDEF_CLASS
#define DEF_MAT_SOLVER_MRTR_MYDEF_CLASS


#include <SparseSolv/MatSolverBase.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/*
//=======================================================
// ■ スパース行列用ソルバのICCGクラス
//=======================================================*/
/**
 * @class MatSolverMRTR
 * @brief MRTR Linear matrix solvers for sparse matrices
 * MRTRを用いた疎行列の線形ソルバクラスです。
 * MatSolverBaseを継承しており、内部状態があるのでインスタンス化して使ってください
 */
class MatSolverMRTR : public MatSolverBase {
private:
	/* 設定 */
	/** @brief Accel factor */
	double accel_factor = 1.05;							/* 加速係数 */
	/* 設定 */
	/** @brief PreconditionType */
	int precondition_type = 1;							/* 前処理タイプ:0=IC, 1=SGS */

	/* 対角スケーリングありの場合の分岐 */
	bool solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results);
	/* MRTRソルバ本体 */
	bool solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const double* diagD, const Eigen::SparseMatrix<double, Eigen::RowMajor>& matA, const SparseMat& matL, const SparseMat& matL_tr, const Eigen::VectorXd& vecB, Eigen::VectorXd& results);

	/* 対角スケーリングありの場合の分岐 */
	bool solve_diag_scale(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results);
	/* MRTRソルバ本体 */
	bool solve_main(const slv_int size, const double conv_cri, const int max_ite, const double normB,
					const dcomplex* diagD, const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& matA, const SparseMatC& matL, const SparseMatC& matL_tr, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results);

	/* SGS前処理バージョンのMRTR・本体 */
	bool solve_SGS_main(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results);
	bool solve_SGS_MRTR_main(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const Eigen::SparseMatrix<double, Eigen::RowMajor>& matA, const SparseMat& matL, const SparseMat& matL_tr, const Eigen::VectorXd& vecB, Eigen::VectorXd& results);
	/* SGS前処理バージョンのMRTR・本体(複素) */
	bool solve_SGS_main(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results);
	bool solve_SGS_MRTR_main(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& matA, const SparseMatC& matL, const SparseMatC& matL_tr, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results);

public:
	MatSolverMRTR()=default;
	~MatSolverMRTR() override = default;
	/*--------------------------------------------------*/
	void setAccelFactor(double accel_factor0) { accel_factor = accel_factor0; };					/* 加速係数の設定 */
	void setPreconditionType(int precondition_type0) { precondition_type = precondition_type0; };	/* 前処理タイプの設定 */
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	bool solve(const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results) override;
	bool solve(const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) override;
	bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double* vecB, double* results) override;
	bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) override;
	/*--------------------------------------------------*/
	bool solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) override;
	bool solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) override;
	bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) override;
	bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) override;
};

/* end of namespace */
};
};

#endif
