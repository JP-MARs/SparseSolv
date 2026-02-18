/**
 * @file MatSolvers_Eigen_unsupported.cpp
 * @brief Linear matrix solvers's wrapper for Eigen in unsupported module
 */

#include <SparseSolve/MatSolversEigenMKL.hpp>
#include <000_thirdparty/Eigen_unsupported/Eigen/IterativeSolvers>


/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
Eigen unsupportedソルバ
//=======================================================
//=======================================================

*/
	
	
/*//=======================================================
// ● GMRESで解く(Eigen unsupported ソルバ)
//=======================================================*/
/** 
 * @brief Eigen's unsupported GMRES solver (preconditioner=ILUT, set_restart=30)
 */
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init){
	/* 設定 */
	Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> gmres;
	gmres.setMaxIterations(max_ite);
	gmres.setTolerance(conv_cri);

	/* トラブル対策2：リスタート設定（メモリ爆発を防ぐ！）
    // これをしないと、計算がどんどん重くなるよ ・・・らしい */
    gmres.set_restart(30); 

	gmres.compute(matA.matrix.matrix);
	if(init){
		results = gmres.solve(vecB);
	}else{
		results = gmres.solveWithGuess(vecB, results);
	}
	double err = gmres.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init){
	Eigen::VectorXd vecB2(size0);
	Eigen::VectorXd results2(size0);
	for(slv_int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
		results2(i) = results[i];
	}
	bool bl = solveEigenGMRES(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const std::vector<double>& vecB, std::vector<double>& results, bool init){
	bool bl = solveEigenGMRES(size0, conv_cri, max_ite, matA, vecB.data(), results.data(), init);
	return bl;
}


	
/*//=======================================================
// ● GMRESで解く(Eigen unsupported ソルバ) for complex
//=======================================================*/
/** 
 * @brief Eigen's unsupported GMRES solver (preconditioner=ILUT, set_restart=30) for complex
 */
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init){
	/* 設定 */
	Eigen::GMRES<Eigen::SparseMatrix<dcomplex>, Eigen::IncompleteLUT<dcomplex>> gmres;
	gmres.setMaxIterations(max_ite);
	gmres.setTolerance(conv_cri);

	/* トラブル対策2：リスタート設定（メモリ爆発を防ぐ！）
    // これをしないと、計算がどんどん重くなるよ ・・・らしい */
    gmres.set_restart(30); 

	gmres.compute(matA.matrix.matrix);
	if(init){
		results = gmres.solve(vecB);
	}else{
		results = gmres.solveWithGuess(vecB, results);
	}
	double err = gmres.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init){
	Eigen::VectorXcd vecB2(size0);
	Eigen::VectorXcd results2(size0);
	for(slv_int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
		results2(i) = results[i];
	}
	bool bl = solveEigenGMRES(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(slv_int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}
bool MatSolversEigenMKL::solveEigenGMRES(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const std::vector<dcomplex>& vecB, std::vector<dcomplex>& results, bool init){
	bool bl = solveEigenGMRES(size0, conv_cri, max_ite, matA, vecB.data(), results.data(), init);
	return bl;
}


/* end of namespace */
};



