#ifndef DEF_SPAR_MAT_D
#define DEF_SPAR_MAT_D


#include "SparseMatTMPL.hpp"
#include "SparseMatC.hpp"
#include "SparseMatOperators.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列（実数）
//=======================================================*/
class SparseMat{	
friend class SparseMatC;
friend class SparseMatOperators;
friend class MatSolvers;
friend class MatSolversEigenMKL;
private:
	SparseMatBaseD matrix;
public:	
	SparseMat(){ ; };
	~SparseMat(){ ; };
	SparseMat(slv_int x);
	SparseMat(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<double>& vals);/* 作成済み疎行列データから初期化 */
	SparseMat(const Eigen::SparseMatrix<double, Eigen::RowMajor>& Mat0);
	SparseMat(const SparseMat& mat);
	SparseMat(SparseMat&& mat) noexcept;
	SparseMat(const SparseMatBaseD& mat);
	SparseMat(SparseMatBaseD&& mat) noexcept;
	SparseMat& operator=(const SparseMat& mat);
	SparseMat& operator=(SparseMat&& mat) noexcept;
	/**/
	bool isFixed() const{return matrix.is_fix;};												/* 確定済みかどうか */
	bool isEmpty() const {return matrix.isEmpty();};											/* 行列の中身が空かどうか */
	void tempInitialize() {matrix.tempInitialize();};											/* 一時vector行列を作成 */
	void fix(bool toSquare=false) { matrix.fix(toSquare); };									/* 一時vector配列を確定させる */
	void refresh(){matrix.refresh();};
	void resetMat() {matrix.resetMat();};																		/* 確定済み行列の値を０に再セット */
	slv_int isInclude(slv_int gyo, slv_int target_r)const{return matrix.isInclude(gyo, target_r);};			/* i行目にtarget_r列があるかどうか（あったらそのindexを返す） */	
	void add(slv_int gyo, slv_int retu, double val){matrix.add(gyo, retu, val);};					/* 一時配列にpush */
	void getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<double>& row_val)const{matrix.getTargetRowVal(target, row_pos, row_val);};			/* 指定した列の非ゼロの行位置と値をvectorに書き出す */
	void getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<double>& col_val)const{matrix.getTargetColVal(target, col_pos, col_val);};			/* 指定した行の非ゼロの列位置と値をvectorに書き出す */
	slv_int getMaxCol()const{return matrix.getMaxCol();};											/* スパース内の最大の列位置を返す */
	void delFlagPosition(const bool* flag){matrix.delFlagPosition(flag);};							/* フラグがTrueの位置の列・行のデータを削除（FEM用） */
	void DiagFlagPosition(const bool* flag){matrix.DiagFlagPosition(flag);};						/* フラグがTrueの位置の列・行のデータをゼロにして対角だけ１に（FEM用） */
	void printMat(const std::string& str="Mat.csv") {matrix.printMat(str);};;
	void print() { matrix.print(); };
	void readMat(const std::string& filename){ matrix.readMat(filename); };
	/* */
	/* オペレータ群 */
	double* operator*(const double* vec) const;	
	dcomplex* operator*(const dcomplex* vec) const;
	Eigen::VectorXd operator*(const Eigen::VectorXd& vec) const;	
	Eigen::VectorXcd operator*(const Eigen::VectorXcd& vec) const;	
	std::vector<double> operator*(const std::vector<double>& vec) const;
	std::vector<dcomplex> operator*(const std::vector<dcomplex>& vec) const;
	void operator*=(const double x) {matrix *= x;};
	SparseMat operator*(const double x) const;
	SparseMatC operator*(const dcomplex x) const;
	SparseMat operator*(const SparseMat& mat) const;
	SparseMatC operator*(const SparseMatC& mat) const;
	SparseMat operator+(const SparseMat& mat) const;
	SparseMatC operator+(const SparseMatC& mat) const;
	/**/
	/* その他オペレータ */
	/**/
	SparseMat trans() const;		/* 転置 */
	SparseMat makeSubMat(slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b);	/* 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、Bに渡す */
	void MatDiv(SparseMat& matK11, SparseMat& matK12, slv_int rangeA, slv_int rangeB);			/* 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する */
	SparseMat getMatLower() const;				/* 下三角行列の取得 */
	SparseMat getMatUpper() const;				/* 上三角行列の取得 */
	SparseMat inv() const;						/* 逆行列（密アルゴリズム・大行列でやるとメモリが吹き飛ぶ！） */
	SparseMat makePrsdInv(double eps) const;	/* 疑似逆行列用の行列（A^T*A+epsI）を作成 */
	void round();								/* 行列値を丸めて整数にする */
	/**/
	/* 不完全コレスキー分解 */
	SparseMat IC_decomp(double* diagD, const double accela) const;
	/* 対角スケーリング */
	SparseMat diagScaling(double* trans_vec, const double* ori_vec) const;
};

/* end of namespace */
};


#endif
