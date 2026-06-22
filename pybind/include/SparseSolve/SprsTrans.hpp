/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 
 /**
 * @file SprsTrans.hpp
 * @brief Matrix transform Operators for Sparse Matrix
 * 
 * このファイルは SparseMatTMPL の行列構造を変換する処理を提供します。
 */
#ifndef DEF_SPARSE_MATTRANSFORM_OPERATOR_DEFINES
#define DEF_SPARSE_MATTRANSFORM_OPERATOR_DEFINES

#include "SparseMat.hpp"


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/* 疎行列変換処理用の専用名前空間 */
/**
* @namespace SprsTrans
* @brief Original namespace for sparse matrix transform operators
*/
namespace SprsTrans {
	/* 行列の値を四捨五入） */
	void round(SparseMat& matA);

	/* 行列Aにa1をかける */
	void SelfDot(SparseMat& matA, double a1);
	void SelfDot(SparseMatC& matA, double a1);

	/* ほぼゼロになっている非ゼロ位置の削除 */
	void prune(SparseMat& matA, double thres);

	/* 転置 */
	SparseMat trans(const SparseMat& matA);
	SparseMatC trans(const SparseMatC& matA);

	/* 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作る */
	SparseMat makeSubMat(const SparseMat& matA, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b);
	SparseMatC makeSubMat(const SparseMatC& matA, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b);

	/* 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する */
	void MatDiv(SparseMat& matK11, SparseMat& matK12, const SparseMat& matA, slv_int rangeA, slv_int rangeB);
	void MatDiv(SparseMatC& matK11, SparseMatC& matK12, const SparseMatC& matA, slv_int rangeA, slv_int rangeB);

	/* 下三角行列の取得 */
	SparseMat MatLower(const SparseMat& matA);
	SparseMatC MatLower(const SparseMatC& matA);
	/* 上三角行列の取得 */
	SparseMat MatUpper(const SparseMat& matA);
	SparseMatC MatUpper(const SparseMatC& matA);


	/* フラグがTrueの位置の列・行のデータを削除（FEM用） */
	void DelFlagPosition(SparseMat& matA, const bool* flag);
	void DelFlagPosition(SparseMatC& matA, const bool* flag);

	/* フラグがTrueの位置の列・行のデータをゼロにして対角だけ１に（FEM用） */
	void DiagFlagPosition(SparseMat& matA, const bool* flag);
	void DiagFlagPosition(SparseMatC& matA, const bool* flag);
}


/* end of namespace */
};
};


#endif


