/*
 * <SparseSolv>
 * Copyright (c) 2026 Takahiro Sato
 *
 * This source code is licensed under the MPL2 License.
 * See the LICENSE file in the root directory for details.
 */

 
 /**
 * @file SprsOperator.hpp
 * @brief Operators for Sparse Matrix
 * 
 * このファイルは SparseMatTMPL への基本的な四則演算を提供します。
 */
#ifndef DEF_SPARSE_OPERATOR_DEFINES
#define DEF_SPARSE_OPERATOR_DEFINES

#include "SparseMat.hpp"


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/* 疎行列演算用の専用名前空間 */
/**
* @namespace SprsOperator
* @brief Original namespace for sparse matrix operators
*/
namespace SprsOperator {
    /* 行列同士の加算 */
    SparseMat add(const SparseMat& matA, const SparseMat& matB);
    SparseMatC add(const SparseMatC& matA, const SparseMat& matB);
    SparseMatC add(const SparseMat& matA, const SparseMatC& matB);
    SparseMatC add(const SparseMatC& matA, const SparseMatC& matB);
    /* 行列同士の引き算 */
    SparseMat substruct(const SparseMat& matA, const SparseMat& matB);
    SparseMatC substruct(const SparseMatC& matA, const SparseMat& matB);
    SparseMatC substruct(const SparseMat& matA, const SparseMatC& matB);
    SparseMatC substruct(const SparseMatC& matA, const SparseMatC& matB);
    /* 行列同士の掛け算 */
    SparseMat dot(const SparseMat& matA, const SparseMat& matB);
    SparseMatC dot(const SparseMatC& matA, const SparseMat& matB);
    SparseMatC dot(const SparseMat& matA, const SparseMatC& matB);
    SparseMatC dot(const SparseMatC& matA, const SparseMatC& matB);

    /*-------------------------------*/
    /*-------------------------------*/
    /*-------------------------------*/
    /* ベクトル掛け算オペレータたち(c++デフォルト) */
	double* dot(const SparseMat& matA, const double* vec);
	dcomplex* dot(const SparseMat& matA, const dcomplex* vec);
	dcomplex* dot(const SparseMatC& matA, const double* vec);
	dcomplex* dot(const SparseMatC& matA, const dcomplex* vec);
	std::vector<double> dot(const SparseMat& matA, const std::vector<double>& vec);
	std::vector<dcomplex> dot(const SparseMatC& matA, const std::vector<double>& vec);
	std::vector<dcomplex> dot(const SparseMat& matA, const std::vector<dcomplex>& vec);
	std::vector<dcomplex> dot(const SparseMatC& matA, const std::vector<dcomplex>& vec);
    /* ベクトル掛け算オペレータたち(Eigen) */
    Eigen::VectorXd dot(const SparseMat& matA, const Eigen::VectorXd& vec);
    Eigen::VectorXcd dot(const SparseMat& matA, const Eigen::VectorXcd& vec);
    Eigen::VectorXcd dot(const SparseMatC& matA, const Eigen::VectorXd& vec);
    Eigen::VectorXcd dot(const SparseMatC& matA, const Eigen::VectorXcd& vec);


    /*-------------------------------*/
    /*-------------------------------*/
    /*-------------------------------*/
    /* 足し算オペレータ(matAにmatBを加える(a1*A+a2*B)。Bの位置を行方向にpos1、列方向にpos2だけずらす。) */
    SparseMat ShiftPus(const SparseMat& matA, const SparseMat& matB, double a1, double a2, slv_int pos1, slv_int pos2);
    SparseMatC ShiftPus(const SparseMatC& matA, const SparseMatC& matB, double a1, double a2, slv_int pos1, slv_int pos2);

    /* matAにa1*matBを加える。 開始位置はposだけずらす */
    void MatSelfPlus(SparseMat& matA, const SparseMat& matB, const double a1, const slv_int pos1, const slv_int pos2);
    void MatSelfPlus(SparseMatC& matA, const SparseMatC& matB, const dcomplex a1, const slv_int pos1, const slv_int pos2);

    /* (mat1 + mat2)*vecBを計算 */
    double* dot2(const SparseMat& matA, const SparseMat& matB, const double* vec);
    dcomplex* dot2(const SparseMatC& matA, const SparseMatC& matB, const double* vec);
    dcomplex* dot2(const SparseMat& matA, const SparseMat& matB, const dcomplex* vec);
    dcomplex* dot2(const SparseMatC& matA, const SparseMatC& matB, const dcomplex* vec);
    std::vector<double> dot2(const SparseMat& matA, const SparseMat& matB, const std::vector<double>& vec);
    std::vector<dcomplex> dot2(const SparseMatC& matA, const SparseMatC& matB, const std::vector<double>& vec);
    std::vector<dcomplex> dot2(const SparseMat& matA, const SparseMat& matB, const std::vector<dcomplex>& vec);
    std::vector<dcomplex> dot2(const SparseMatC& matA, const SparseMatC& matB, const std::vector<dcomplex>& vec);
    /* ベクトル掛け算オペレータたち(Eigen) */
    Eigen::VectorXd dot2(const SparseMat& matA, const SparseMat& matB, const Eigen::VectorXd& vec);
    Eigen::VectorXcd dot2(const SparseMatC& matA, const SparseMatC& matB, const Eigen::VectorXd& vec);
    Eigen::VectorXcd dot2(const SparseMat& matA, const SparseMat& matB, const Eigen::VectorXcd& vec);
    Eigen::VectorXcd dot2(const SparseMatC& matA, const SparseMatC& matB, const Eigen::VectorXcd& vec);

    /* 行列AとBとCをかけて結果を返す */
	SparseMat dot3(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC);
    SparseMatC dot3(const SparseMat& matA, const SparseMatC& matB, const SparseMat& matC);
    SparseMatC dot3(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC);

    /* 行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす */
    SparseMat MatPlusMinusShift(const SparseMat& matA, const SparseMat& matB, double a1, double a2, slv_int pos1, slv_int pos2);
    SparseMatC MatPlusMinusShift(const SparseMatC& matA, const SparseMatC& matB, double a1, double a2, slv_int pos1, slv_int pos2);

    /* 行列A±B±Cを加える(a1*A+a2*B+a3*C) */
    SparseMat add3(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC, double a1, double a2, double a3);
    SparseMatC add3(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC, double a1, double a2, double a3);

    /* 行列Aの逆行列を求めてInvに渡す(密用) */
	SparseMat inv(const SparseMat& targetMat);
	SparseMatC inv(const SparseMatC& targetMat);
}


/* end of namespace */
};
};


#endif


