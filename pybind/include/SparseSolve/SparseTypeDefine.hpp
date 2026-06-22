/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 
 /**
 * @file SparseTypeDefine.hpp
 * @brief Sparse Matrix Type Definitions
 * このファイルはSparseSolvで使う型を定義します
 */
#ifndef DEF_SPR_MAT_DEFINE_TYPE_USING_DEF
#define DEF_SPR_MAT_DEFINE_TYPE_USING_DEF


/* オリジナル名前空間(JPMARsライブラリ) */
/**
* @namespace JPMRspace
* @brief Original namespace for JP-MARs
*/
namespace JPMRspace{

	
/* 疎行列ソルバ用オリジナル名前空間 */
/**
* @namespace SparseSolv
* @brief Original namespace for sparse matrix solvers
*/
namespace SparseSolv{



/* ソルバ内のint定義 */
/**
* @typedef slv_int
* @brief Alias for integer type used in solvers
*/
using slv_int = int;





/* end of namespace */
};
};

#endif
