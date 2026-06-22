/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 /**
 * @file SparseBuilder.hpp
 * @brief Sparse Matrix Builder Template Class
 * このファイルはテンプレートクラス SparseBuilder を提供します。
 * 内部では一時保存用の map[] を用いて要素を構築し、
 * 疎行列に必要な情報を構築します。
 * このクラスを疎行列クラスに渡すことで疎行列を完成させます
 */


#ifndef DEF_SPR_MAT_BUILDER_TMPL
#define DEF_SPR_MAT_BUILDER_TMPL

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

//#include <BasicDefines.hpp>
#include "BasicDefinesC.hpp"
#include "BasicDefinesEigenSparse.hpp"

#include "SparseTypeDefine.hpp"

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{

/*
//=======================================================
// ■ スパース行列構築補助クラステンプレート（自作形式）
//=======================================================*/
/**
 * @class SparseBuilderTMPL
 * @brief Sparse matrix builder template
 * このファイルはテンプレートクラス SparseBuilderTMPL を提供します。
 * 内部では一時保存用の map[] を用いて要素を構築し、
 * 疎行列に必要な情報を構築します。
 * このクラスを疎行列クラスに渡すことで疎行列を完成させます
 * @tparam DType （double, dcomplex）
 */
template<typename DType> class SparseBuilderTMPL{
private:
	/** @brief row size */
	slv_int size;																/* 行列の行数 */
	/** @brief temporary storage for constructing the matrix */
	std::vector<std::map<slv_int, DType>> tempMat;								/* 一時保存行列（キーが列位置） */	
	/*  */
public:
	SparseBuilderTMPL()=default;
	SparseBuilderTMPL(slv_int size0){this->initialize(size0);};											/* コンストラクタ */
	~SparseBuilderTMPL()=default;	
	SparseBuilderTMPL(const SparseBuilderTMPL<DType>& Mat)=default;										/* コピーコンストラクタ */
	SparseBuilderTMPL(SparseBuilderTMPL<DType>&& Mat) noexcept = default;								/* ムーブコンストラクタ */
	SparseBuilderTMPL<DType>& operator=(const SparseBuilderTMPL<DType>& Mat) = default;					/* 代入オペレータ */
	SparseBuilderTMPL<DType>& operator=(SparseBuilderTMPL<DType>&& Mat) noexcept = default;				/* 代入オペレータ（ムーブ） */
	/**/
	slv_int getSize() const{return size;};
	/**/
	void initialize(slv_int ss);																/* 一時行列を作成 */
	void initialize();																			/* 一時行列を作成 */
	void reInitialize(slv_int new_size){this->initialize(new_size);};							/* 一時行列を作成&初期化 */
	void resetMat();																			/* 構築済み部分の値を０に再セット */
	void add(slv_int gyo, slv_int retu, DType val){tempMat[gyo][retu] += val;};					/* 一時配列にpush */

	slv_int build(std::vector<Eigen::Triplet<DType>>& tripletList, bool to_square);				/* 一次配列からtripletを作り、一次配列を削除 */

	void printTempMat()const;																	/* 一次配列の状況をprint */
};

/* ===== usingで型隠し ===== */

/**
 * @page SparseBuilderTMPL_aliases SparseBuilderTMPL Aliases
 * @brief Sparse matrix builder aliases used in this library.
 *
 * This library provides the following type aliases:
 *
 * - `SparseBuilder` : alias of `SparseBuilderTMPL<double>`
 * - `SparseBuilderC`: alias of `SparseBuilderTMPL<dcomplex>`
 *
 * @section sparsebuilder_alias_real SparseBuilder
 * `SparseBuilder` is the double precision sparse matrix builder type.
 *
 * @code{.cpp}
 * using SparseBuilder = SparseBuilderTMPL<double>;
 * @endcode
 *
 * @section sparsebuilder_alias_complex SparseBuilderC
 * `SparseBuilderC` is the complex sparse matrix builder type.
 *
 * @code{.cpp}
 * using SparseBuilderC = SparseBuilderTMPL<dcomplex>;
 * @endcode
 */


/** 
 * @typedef SparseBuilder
 * @brief Base Sparse matrix builder for double precision
 *
 */
using SparseBuilder = SparseBuilderTMPL<double>;

/** 
 * @typedef SparseBuilderC
 * @brief Base Sparse matrix builder for complex<double> precision
 *
 */
using SparseBuilderC = SparseBuilderTMPL<dcomplex>;


/* ========== */
extern template class SparseBuilderTMPL<double>;
extern template class SparseBuilderTMPL<dcomplex>;



/* end of namespace */
};
};

#endif
