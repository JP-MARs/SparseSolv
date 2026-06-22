/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 /**
 * @file MatSolver.hpp
 * @brief Linear matrix solvers for sparse matrices
 * このファイルは疎行列の線形ソルバ用の基底クラスを提供します
 * 加速係数付きのICCG法やMRTR法の基底クラスを実装しています。
 */
#ifndef DEF_MAT_SOLVER_MYDEF_BASECLASS
#define DEF_MAT_SOLVER_MYDEF_BASECLASS


#include <omp.h>

#include <SparseSolv/SparseMat.hpp>
#include <SparseSolv/SprsOperator.hpp>
#include <SparseSolv/SprsTrans.hpp>
#include <SparseSolv/SprsPreconditioner.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/*
//=======================================================
// ■ スパース行列用ソルバの基底クラス
//=======================================================*/
/**
 * @class MatSolverBase
 * @brief Base Class of Linear matrix solvers for sparse matrices
 * このファイルは疎行列の線形ソルバの基底クラスを提供します。
 * 収束判定用の内部情報や前進後退代入処理、加速係数の自動決定処理などを実装しています。
 */
class MatSolverBase{
protected:
	struct MatrixInfoD{
		std::vector<slv_int> start_pos;	/* 各行の列のスタート位置 */
		std::vector<slv_int> end_pos;	/* 各行の列の終了位置 */
		slv_int* col_ptr;				/* 各列番号の生ポインタ */
		double* val_ptr;				/* 各値の生ポインタ */
	};
	struct MatrixInfoC{
		std::vector<slv_int> start_pos;	/* 各行の列のスタート位置 */
		std::vector<slv_int> end_pos;	/* 各行の列の終了位置 */
		slv_int* col_ptr;				/* 各列番号の生ポインタ */
		dcomplex* val_ptr;				/* 各値の生ポインタ */
	};

	/*--------------------------------------------------*/
	bool is_init;									/* 解ベクトルをゼロクリアするかどうか */
	/** @brief save best solution? */
	/** @brief appy diagonal scaling? */
	bool is_diag_scale;								/* 対角スケーリングを行うか */	
	bool is_save_best;								/* 最良結果を保存して失敗時に利用するか */
	/** @brief save residual log? */
	bool is_save_residual_log;						/* 残差履歴を残すか */
	/** @brief residual log data */
	std::vector<double> residual_log;				/* 残差履歴 */
	/* 発散判定：
	  0=最大反復までやる
	  1：最良×bad_div_valより大きい値がbad_div_countだけ続いたら発散として終わる */
	/** @brief divergence check type
	* 0: no divergence check
	* 1: divergence check by bad divergence value and count threshold
	*/
	int diverge_judge_type;
	/** @brief bad divergence value  */
	double bad_div_val;
	/** @brief divergence check thres */
	int bad_div_count_thres;
	/* 収束判定の正規化タイプ(0:右辺、1:初期残差、2:外部指定) */
	/** @brief  convergence check type 
	* 0: normalize by norm of rhs
	* 1: normalize by initial residual norm
	* 2: normalize by user specified value
	*/
	int conv_normalize_type;
	/** @brief user specified normalization constant */
	double conv_normalize_const;	/* 外部指定する正規化係数 */
	/* 絶対収束の最小値 */  
	/** @brief small absolute convergence value */
	static constexpr double small_abs_conv_val = 1.0e-20;
	/**/
	/*--------------------------------------------------*/
	/* 収束判定正規化定数の計算 */
	double setNormalizer(const double normB, const double pure_normR);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* 対角スケーリング処理時の行列変換 */
	SparseMat makeDiagMatTrans(const SparseMat& matA, const Eigen::VectorXd& matDiagD);
	SparseMatC makeDiagMatTrans(const SparseMatC& matA, const Eigen::VectorXcd& matDiagD);
	/*--------------------------------------------------*/
	/* 疎行列積用の対角情報の事前抜き取り */
	void makeMatrixInfo(const SparseMat& matA, struct MatrixInfoD& mat_info);
	void makeMatrixInfo(const SparseMatC& matA, struct MatrixInfoC& mat_info);
	/* 単純な下三角の前進代入処理 */
	void fr_process(const slv_int size0, const struct MatrixInfoD& mat_infoL, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void fr_process(const slv_int size0, const struct MatrixInfoC& mat_infoL, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/* 単純な下三角(転置)の後退代入処理 */
	void bc_process(const slv_int size0, const struct MatrixInfoD& mat_infoL_tr, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void bc_process(const slv_int size0, const struct MatrixInfoC& mat_infoL_tr, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/* IC分解の前進後退代入処理 */
	void IC_frbc_process(const slv_int size0, const struct MatrixInfoD& mat_infoL, const struct MatrixInfoD& mat_infoL_tr, const double* diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void IC_frbc_process(const slv_int size0, const struct MatrixInfoC& mat_infoL, const struct MatrixInfoC& mat_infoL_tr, const dcomplex* diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* 加速係数の自動決定 */
	void auto_accel_determine(const slv_int size0, double accel_ini, const SparseMat& matA, double* diagD, SparseMat& matL);
	void auto_accel_determine(const slv_int size0, double accel_ini, const SparseMatC& matA, dcomplex* diagD, SparseMatC& matL);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* 内積の独自ラッパー */
	double dot_no_conj(const Eigen::VectorXd& vec1, const Eigen::VectorXd& vec2) { return vec1.dot(vec2); };
	dcomplex dot_no_conj(const Eigen::VectorXcd& vec1, const Eigen::VectorXcd& vec2) { return (vec1.array() * vec2.array()).sum(); };

public:
	MatSolverBase();
	virtual ~MatSolverBase()=default;
	/**/
	void setInitSolution(bool bl){ is_init = bl; };
	void setDiagScale(bool bl){is_diag_scale=bl;};
	void setSaveBest(bool bl){ is_save_best = bl; };
	void setSaveLog(bool bl){ is_save_residual_log = bl; };
	void getResidualLog(std::vector<double>& log);
	void setDirvegeType(int x){ diverge_judge_type = x; };
	void setBadDivVal(double x){bad_div_val=x;};
	void setBadDivCount(int x){bad_div_count_thres=x;};
	void setConvNormalizeType(int x){ conv_normalize_type = x; };
	void setConvNormalizeConst(double x){ conv_normalize_const = x; };
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	virtual bool solve(const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results) = 0;		
	virtual bool solve(const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) = 0;
	virtual bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double* vecB, double* results) = 0;
	virtual bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results) = 0;
	/*--------------------------------------------------*/
	virtual bool solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) = 0;		
	virtual bool solve(const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) = 0;
	virtual bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results) = 0;
	virtual bool solve(const double conv_cri, const int max_ite, const double normB, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results) = 0;
	
};

/* end of namespace */
};
};

#endif
