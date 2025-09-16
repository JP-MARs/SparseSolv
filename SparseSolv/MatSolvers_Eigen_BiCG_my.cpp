#include "MatSolversEigenMKL.hpp"
#include <000_thirdparty/Eigen/IterativeLinearSolvers>
#include <000_thirdparty/Eigen/SparseCholesky>	

#include "IncompleteLUT_my.h"


/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
Eigenソルバ for BiCG自作修正（加速係数付き）
//=======================================================
//=======================================================

*/



/*//=======================================================
// ● BiCGstabで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenBiCGstab_accel(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init){
	/* 設定 */
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT_my<double>> bicg;
	bicg.setMaxIterations(max_ite);
	bicg.setTolerance(conv_cri);

	bicg.compute(matA.matrix.matrix);
	if(init){
		results = bicg.solve(vecB);
	}else{
		results = bicg.solveWithGuess(vecB, results);
	}
	double err = bicg.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenBiCGstab_accel(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init){
	Eigen::VectorXd vecB2(size0);
	Eigen::VectorXd results2(size0);
	for(slv_int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenBiCGstab(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}

/*//=======================================================
// ● BiCGstabで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenBiCGstab_accel(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init){
	/* 設定 */
	Eigen::BiCGSTAB<Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>, Eigen::IncompleteLUT_my<dcomplex>> bicg;
	bicg.setMaxIterations(max_ite);
	bicg.setTolerance(conv_cri);

	bicg.compute(matA.matrix.matrix);
	if(init){
		results = bicg.solve(vecB);
	}else{
		results = bicg.solveWithGuess(vecB, results);
	}
	double err = bicg.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenBiCGstab_accel(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init){
	Eigen::VectorXcd vecB2(size0);
	Eigen::VectorXcd results2(size0);
	for(slv_int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenBiCGstab(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}


/* end of namespace */
};



