#ifndef DEF_MAT_SOLVER_ICCG_MYDEF
#define DEF_MAT_SOLVER_ICCG_MYDEF

/* 単純なICCG並列化(行列ベクトル積と内積の単純なOMP化)するとき、on */
//#define OMP_USING_ICCG

#include "SparseMat.hpp"
#include "SparseMatC.hpp"
#include <cfloat>
#include <omp.h>


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
namespace py = pybind11;

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列用ソルバ
//=======================================================*/
class MatSolvers{
private:
	/*--------------------------------------------------*/
	/* 設定 */
	bool is_diag_scale;								/* 対角スケーリングを行うか */
	bool is_save_best;								/* 最良結果を保存して失敗時に利用するか */
	bool is_save_residual_log;						/* 残差履歴を残すか */
	std::vector<double> residual_log;				/* 残差履歴 */
	/* 発散判定：
	  0=最大反復までやる
	  1：最良×bad_div_valより大きい値がbad_div_countだけ続いたら発散として終わる */
	int diverge_judge_type;
	double bad_div_val;
	int bad_div_count_thres;
	/* 収束判定の正規化タイプ(0:右辺、1:初期残差、2:外部指定) */
	int conv_normalize_type;
	double conv_normalize_const;	/* 外部指定する正規化係数 */
	/* 絶対収束の最小値 */  
	static constexpr double small_abs_conv_val = 1.0e-12;
	/**/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* 単純な下三角の前進代入処理 */
	void fr_process(const slv_int size0, const SparseMatBaseD& matL, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void fr_process(const slv_int size0, const SparseMatBaseD& matL, const double* vecR, double* vec);	
	void fr_process(const slv_int size0, const SparseMatBaseC& matL, const dcomplex* vecR, dcomplex* vec);	
	void fr_process(const slv_int size0, const SparseMatBaseC& matL, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/* 単純な下三角(転置)の後退代入処理 */
	void bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void bc_process(const slv_int size0, const SparseMatBaseD& matL_tr, const double* vecR, double* vec);	
	void bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const dcomplex* vecR, dcomplex* vec);	
	void bc_process(const slv_int size0, const SparseMatBaseC& matL_tr, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/* IC分解の前進後退代入処理 */
	void IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double* diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void IC_frbc_process(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double* diagD, const double* vecR, double* vec);	
	void IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex* diagD, const dcomplex* vecR, dcomplex* vec);	
	void IC_frbc_process(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex* diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* 加速係数の自動決定 */
	void auto_accel_determine(const slv_int size0, double accel_ini, const SparseMat& matA, double* diagD, SparseMat& matL);
	void auto_accel_determine(const slv_int size0, double accel_ini, const SparseMatC& matA, dcomplex* diagD, SparseMatC& matL);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/**/
	/* ICCG専用内部処理(対角スケーリングあり) */
	bool solveICCG_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	/* ICCG専用内部処理 */
	bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double* vecB, double* results, bool init=false);
	/**/
	/* ICCOCG専用内部処理(対角スケーリングあり) */
	bool solveICCG_diag(slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	/* ICCOCG専用内部処理 */
	bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
						  const dcomplex* diagD, const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex* vecB, dcomplex* results, bool init=false);
	/**/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--By Shingo Hiruma-------------------------------------------*/
	/* ABMC ordering付きICCG内部処理*/
	//
	// Sort by Algebraic Multi-Color Ordering
	void sortAlgebraicMultiColor(slv_int size0, const SparseMatBaseD &A, SparseMatBaseD &PA, const double *b, double *Pb,
		std::vector<std::vector<slv_int>> &color_list, std::vector<slv_int> &ordering, std::vector<slv_int> &reverse_ordering,
		int num_colors);
	
	// Make Algebraic Block
	void makeAlgebraicBlock(slv_int size0, const SparseMatBaseD &A, SparseMatBaseD &Mb,
		std::vector<std::vector<slv_int>> &block_list, int num_blocks);
	
	// Sort by Algebraic Block Multi-Color Ordering
	void sortAlgebraicBlockMultiColor(slv_int size0, const SparseMatBaseD &A, SparseMatBaseD &PA, const double *b, double *Pb,
		std::vector<std::vector<slv_int>> &block_list, std::vector<std::vector<slv_int>> &color_list, std::vector<slv_int> &ordering, std::vector<slv_int> &reverse_ordering,
		int num_blocks, int num_colors);

	// ICCG Method using ABMC ordering
	bool parallelIccgSolv(slv_int size0, const SparseMatBaseD &A, const SparseMatBaseD &L, const SparseMatBaseD &Lt, 
		std::vector<std::vector<slv_int>> &block_list, std::vector<std::vector<slv_int>> &color_list,
		const double *diagD, const double *b, double *x, double eps, int numItr);

	// solve [L][D][Lt]{x} = {b} for ICCG solver using ABMC ordering
	void parallelIcSolv(slv_int size0, const SparseMatBaseD &L, const SparseMatBaseD &Lt, 
		std::vector<std::vector<slv_int>> &block_list, std::vector<std::vector<slv_int>> &color_list,
		const double *diagD, const Eigen::VectorXd &b, Eigen::VectorXd &x);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/**/
	/* MRTR専用内部処理(対角スケーリングあり) */
	bool solveICMRTR_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	/* MRTR専用内部処理 */
	bool solveICMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double* vecB, double* results, bool init=false);
	/**/
	/* CO-MRTR専用内部処理(対角スケーリングあり) */
	bool solveICMRTR_diag(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	/* CO-MRTR専用内部処理 */
	bool solveICMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const dcomplex* diagD, const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex* vecB, dcomplex* results, bool init=false);
	/**/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* SGS-MRTR専用内部処理 */
	bool solveSGSMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double* vecB, double* results, bool init=false);
	/* SGS--CO-MRTR専用内部処理 */
	bool solveSGSMRTR(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex* vecB, dcomplex* results, bool init=false);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
public:
	MatSolvers();		/* コンストラクタ */

	void setDiagScale(bool bl){ is_diag_scale = bl; };
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
	/* ICCG法 */
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init=false);	
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init=false);
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false);
	/* ICCOCG法 */
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init=false);
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--By Shingo Hiruma-------------------------------------------*/
	/* ABMC ICCG法*/
	bool solveICCGwithABMC(const slv_int size0, const double conv_cri,  const int max_ite, const double accera, const double normB, const SparseMat &matA, 
		const double *vec_b, double *vec_x, int num_blocks, int num_colors, bool init=false);
	bool solveICCGwithABMC(const slv_int size0, const double conv_cri,  const int max_ite, const double accera, const SparseMat &matA, 
						   const double *vec_b, double *vec_x, int num_blocks, int num_colors, bool init=false);
	bool solveICCGwithABMC(const slv_int size0, const double conv_cri,  const int max_ite, const double accera, const SparseMat &matA, 
						   const std::vector<double>& vec_b, std::vector<double>& vec_x, int num_blocks, int num_colors, bool init=false);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/* Based on Tsuburaya Tomonori */
	/* IC MRTR法 */
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init=false);	
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init=false);
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false);
	/* ICCO MRTR法 */
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init=false);
	bool solveICMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
	/**/
	/* SGS MRTR法 */
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double *vecB, double *results, bool init=false);	
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init=false);
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false);
	/* SGS CO MRTR法 */
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init=false);
	bool solveSGSMRTR(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	/*--------------------------------------------------*/
	std::vector<double> getResidualLog_py();
	bool solveICCG_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, py::array_t<double> vecB, py::array_t<double> results, bool init=false);
	bool solveICCG_py(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init=false);
	bool solveICCGwithABMC_py(const slv_int size0, const double conv_cri,  const int max_ite, const double accera, const SparseMat &matA, 
						   py::array_t<double> vec_b, py::array_t<double> results, int num_blocks, int num_colors, bool init=false);
	bool solveICMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, py::array_t<double> vecB, py::array_t<double> results, bool init=false);
	bool solveICMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init=false);
	bool solveSGSMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA,  py::array_t<double> vecB, py::array_t<double> results, bool init=false);	
	bool solveSGSMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init=false);

};

/* end of namespace */
};

#endif
