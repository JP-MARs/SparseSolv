/**
 * @file SparseMatC
 * @brief Sparse Matrix for complex<double>
 *
 * このファイルはcomplex<double>用スパース行列クラス SparseMatC を定義しています。
 * 内部でSparseMatBaseCを保持し、行列の生成・操作・演算を行います。行列の実態はSparseMatBaseCが管理します。
 */
#ifndef DEF_SPAR_MAT_C
#define DEF_SPAR_MAT_C


#include "SparseMatTMPL.hpp"
#include "SparseMat.hpp"
#include "SparseMatOperators.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列（実数）
//=======================================================*/
/**
 * @class SparseMatC
 * @brief Main Sparse Matrix for complex<double>
 * このファイルは本ライブラリが提供する複素実数の疎行列 SparseMatC を定義します。
 * 内部でSparseMatBaseCを保持し、行列の生成・操作・演算を行います。行列の実態はSparseMatBaseCが管理します。
 * このクラスはSparseMatBaseCに対する演算等を定義しています。
 */
class SparseMatC{	
friend class SparseMat;
friend class SparseMatOperators;
friend class MatSolvers;
friend class MatSolversEigenMKL;
private:
	/** @brief innner sparse matrix data */
	SparseMatBaseC matrix;
public:	
	SparseMatC()=default;
	~SparseMatC()=default;
	SparseMatC(slv_int x);
	SparseMatC(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<dcomplex>& vals);/* 作成済み疎行列データから初期化 */
	SparseMatC(const Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>& Mat0);
	SparseMatC(const SparseMatC& mat);
	SparseMatC(SparseMatC&& mat) noexcept;
	SparseMatC(const SparseMatBaseC& mat);
	SparseMatC(SparseMatBaseC&& mat) noexcept;
	SparseMatC& operator=(const SparseMatC& mat);
	SparseMatC& operator=(SparseMatC&& mat) noexcept;
	/**/
	bool isFixed() const{return matrix.is_fix;};												/* 確定済みかどうか */
	bool isEmpty() const {return matrix.isEmpty();};											/* 行列の中身が空かどうか */
	void tempInitialize() {matrix.tempInitialize();};											/* 一時vector行列を作成 */
	void all_initialize(slv_int x);																/* 行列の完全初期化 */
	void fix(bool toSquare=false) { matrix.fix(toSquare); };									/* 一時vector配列を確定させる */
	void refresh(){matrix.refresh();};
	void resetMat() {matrix.resetMat();};																								/* 確定済み行列の値を０に再セット */
	bool isInclude(slv_int gyo, slv_int target_r, slv_int& result_retu)const{return matrix.isInclude(gyo, target_r, result_retu);};		/* i行目にtarget_r列があるかどうか（あったらそのindexを返す） */	
	void add(slv_int gyo, slv_int retu, dcomplex val){matrix.add(gyo, retu, val);};					/* 一時配列にpush */
	void getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<dcomplex>& row_val)const{matrix.getTargetRowVal(target, row_pos, row_val);};			/* 指定した列の非ゼロの行位置と値をvectorに書き出す */
	void getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<dcomplex>& col_val)const{matrix.getTargetColVal(target, col_pos, col_val);};			/* 指定した行の非ゼロの列位置と値をvectorに書き出す */
	slv_int getMaxCol()const{return matrix.getMaxCol();};											/* スパース内の最大の列位置を返す */
	void puruned(double thres = 1.0e-12){matrix.pruned(thres);};
	void back_unfixed(){ matrix.back_unfixed(); };
	void delFlagPosition(const bool* flag){matrix.delFlagPosition(flag);};							/* フラグがTrueの位置の列・行のデータを削除（FEM用） */
	void DiagFlagPosition(const bool* flag){matrix.DiagFlagPosition(flag);};						/* フラグがTrueの位置の列・行のデータをゼロにして対角だけ１に（FEM用） */
	void printMat(const std::string& str="Mat.csv") {matrix.printMat(str);};;
	void print() { matrix.print(); };
	void readMat(const std::string& filename){ matrix.readMat(filename); };
	/* */
	/* オペレータ群 */
	dcomplex* operator*(const double* vec) const;	
	Eigen::VectorXcd operator*(const Eigen::VectorXd& vec) const;	
	Eigen::VectorXcd operator*(const Eigen::VectorXcd& vec) const;	
	dcomplex* operator*(const dcomplex* vec) const;
	void operator*=(const double x) {matrix *=(x);};
	void operator*=(const dcomplex x) {matrix *=(x);};
	SparseMatC operator*(const double x) const;
	SparseMatC operator*(const dcomplex x) const;
	SparseMatC operator*(const SparseMat& mat) const;
	SparseMatC operator*(const SparseMatC& mat) const;
	SparseMatC operator+(const SparseMat& mat) const;
	SparseMatC operator+(const SparseMatC& mat) const;
	/**/
	/* その他オペレータ */
	/**/
	SparseMatC trans() const;		/* 転置 */
	SparseMatC makeSubMat(slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b);	/* 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、Bに渡す */
	void MatDiv(SparseMatC& matK11, SparseMatC& matK12, slv_int rangeA, slv_int rangeB);	/* 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する */
	SparseMatC getMatLower() const;				/* 下三角行列の取得 */
	SparseMatC getMatUpper() const;				/* 上三角行列の取得 */
	SparseMatC inv() const;						/* 逆行列（密アルゴリズム・大行列でやるとメモリが吹き飛ぶ！） */
	SparseMatC makePrsdInv(double eps) const;	/* 疑似逆行列用の行列（A^T*A+epsI）を作成 */
	/* 不完全コレスキー分解 */
	SparseMatC IC_decomp(dcomplex* diagD, const double accela) const;
	/* 対角スケーリング */
	SparseMatC diagScaling(dcomplex* trans_vec, const dcomplex* ori_vec) const;
	/**/
};

/* end of namespace */
};

#endif
