/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 /**
 * @file SparseMatTMPL.hpp
 * @brief Sparse Matrix Template Class
 * このファイルはテンプレートクラス SparseMatTMPL を提供します。
 * 内部ではSparseBuilderで構築した疎行列情報を読み取り、
 * Eigen::SparseMatrix を構築します。
 * 行列が完成したあとは Eigen の高速な演算が利用可能になります。
 */
#ifndef DEF_SPR_MAT_DEFINE_USING_TMPL
#define DEF_SPR_MAT_DEFINE_USING_TMPL

#include <SparseSolv/SparseBuilder.hpp>
#include <string>
#include <fstream>



/* 行列積などをOpenMP並列化するとき、onしてください */
//#define OMP_USING_MAT_SOL


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/*
//=======================================================
// ■ スパース行列テンプレート
//=======================================================*/
/**
 * @class SparseMatTMPL
 * @brief Sparse matrix template (this is a main matrix class)
 * このファイルはテンプレートクラス SparseMatTMPL を提供します。
 * 内部ではSparseBuilderで構築した疎行列情報を読み取り、
 * Eigen::SparseMatrix を構築します。
 * 行列が完成したあとは Eigen の高速な演算が利用可能になります。
 * @tparam DType （double, dcomplex）
 */
template<typename DType> class SparseMatTMPL{
private:
	/** @brief Eigen sparse matrix */
	Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix;											/* Eigenスパース行列 */
	/*  */
public:
	SparseMatTMPL()=default;
	~SparseMatTMPL()=default;
	SparseMatTMPL(const SparseMatTMPL<DType>& Mat)=default;										/* コピーコンストラクタ */
	SparseMatTMPL(SparseMatTMPL<DType>&& Mat) noexcept = default;								/* ムーブコンストラクタ */
	SparseMatTMPL<DType>& operator=(const SparseMatTMPL<DType>& Mat) = default;					/* 代入オペレータ */
	SparseMatTMPL<DType>& operator=(SparseMatTMPL<DType>&& Mat) noexcept = default;				/* 代入オペレータ（ムーブ） */

	SparseMatTMPL(SparseBuilderTMPL<DType>& builder, bool to_square=false){this->build(builder, to_square);};						/* コンストラクタ(Builderから作成):to_square=長方形の場合に正方形にするか */
	SparseMatTMPL(const Eigen::SparseMatrix<DType, Eigen::RowMajor>& Mat0){matrix=Mat0;};
	SparseMatTMPL(Eigen::SparseMatrix<DType, Eigen::RowMajor>&& Mat0){matrix=std::move(Mat0);};
	SparseMatTMPL(const std::vector<Eigen::Triplet<DType>>& sparse_data);																	/* Eigenタプルから初期化 */
	SparseMatTMPL(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<DType>& vals);	/* 作成済み疎行列データから初期化 */
	/**/
	void build(SparseBuilderTMPL<DType>& builder, bool to_square=false);							/* Builderから作成: to_square=長方形の場合に正方形にするか */
	/**/
	slv_int getSize() const;																	/* 行数を返す(サイズ範囲チェックありで) */
	Eigen::Index getTrueSize() const{return matrix.rows();};									/* 行数を返す */
	slv_int getMaxCol() const;																	/* スパース内の最大の列位置を返す */
	void add(slv_int gyo, slv_int retu, DType val){matrix.coeffRef(gyo, retu) += val;};			/* 確定済み行列にpush */
	void compressed() { matrix.makeCompressed(); };												/* 行列を圧縮する（非ゼロ位置が変わった場合は呼ぶ。それ以外は不要 */

	Eigen::SparseMatrix<DType, Eigen::RowMajor>& mat() noexcept{ return matrix; };				/* 行列を返す */
	const Eigen::SparseMatrix<DType, Eigen::RowMajor>& mat() const noexcept{ return matrix; };	/* 行列を返す */
	void setMatrix(Eigen::SparseMatrix<DType, Eigen::RowMajor>&& mat) { matrix = std::move(mat); };

	void getCols(slv_int* cols_data1, slv_int* cols_data2) const;								/* 各行の列のスタート・終了位置を返す */
	void getCols(std::vector<slv_int>& cols_data1, std::vector<slv_int>& cols_data2) const;		/* 各行の列のスタート・終了位置を返す */
	slv_int* getColPtr(){ return matrix.innerIndexPtr(); };										/* 列のポインタを返す */
	DType* getValuePtr(){ return matrix.valuePtr(); };											/* 値のポインタを返す */
	const slv_int* getColPtr()const{ return matrix.innerIndexPtr(); };							/* 列のポインタを返す(const) */
	const DType* getValuePtr()const{ return matrix.valuePtr(); };								/* 値のポインタを返す(const) */
	void resetMat();																			/* 値だけゼロにする（位置は保持したまま） */

	bool isInclude(slv_int gyo, slv_int target_r, slv_int& result_retu) const;								/*  i行目にtarget_r列があるかどうか（あったらそのindexを返す） */
	void getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<DType>& row_val) const; /* 指定した列の非ゼロの行位置と値をvectorに書き出す  */
	void getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<DType>& col_val)const;	/* 指定した行の非ゼロの列位置と値をvectorに書き出す */

	void prune(double thres);																				/* ほぼゼロになっている非ゼロ位置の削除 */

	void printMat(const std::string& str="Mat.csv") const;
	void print() const { std::cout << matrix << std::endl;};
	void readMat(const std::string& filename);

};

/**
 * @page sparsemat_aliases Sparse Matrix Aliases
 * @brief Sparse matrix type aliases used in this library.
 *
 * This library provides the following type aliases:
 *
 * - `SparseMat` : alias of `SparseMatTMPL<double>`
 * - `SparseMatC`: alias of `SparseMatTMPL<dcomplex>`
 *
 * @section sparsemat_alias_real SparseMat
 * `SparseMat` is the double precision sparse matrix type.
 *
 * @code{.cpp}
 * using SparseMat = SparseMatTMPL<double>;
 * @endcode
 *
 * @section sparsemat_alias_complex SparseMatC
 * `SparseMatC` is the complex sparse matrix type.
 *
 * @code{.cpp}
 * using SparseMatC = SparseMatTMPL<dcomplex>;
 * @endcode
 */

/** 
 * @brief Base Sparse matrix (double precision)
 *
 * inner double sparse matrix class
 */
using SparseMat = SparseMatTMPL<double>;

/** 
 * @brief Base Sparse matrix (double complex)
 *
 * inner complex<double> sparse matrix class
 */
using SparseMatC = SparseMatTMPL<dcomplex>;


/* ========== */
extern template class SparseMatTMPL<double>;
extern template class SparseMatTMPL<dcomplex>;





/* end of namespace */
};
};

#endif
