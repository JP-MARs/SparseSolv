
#include "MatSolversEigenMKL.hpp"
#include <000_thirdparty/Eigen/IterativeLinearSolvers>
#include <000_thirdparty/Eigen/SparseCholesky>	


/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
Eigenソルバ
//=======================================================
//=======================================================

*/
	
	
/*//=======================================================
// ● ICCGで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init){
	/* 設定 */
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double, Eigen::Lower, Eigen::NaturalOrdering<int>>> iccg;
	iccg.setMaxIterations(max_ite);
	iccg.setTolerance(conv_cri);

	iccg.compute(matA.matrix.matrix);
	if(init){
		results = iccg.solve(vecB);
	}else{
		results = iccg.solveWithGuess(vecB, results);
	}
	double err = iccg.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init){
	Eigen::VectorXd vecB2(size0);
	Eigen::VectorXd results2(size0);
	for(int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenICCG(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}

/*//=======================================================
// ● BiCGstabで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init){
	/* 設定 */
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>> bicg;
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
bool MatSolversEigenMKL::solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init){
	Eigen::VectorXd vecB2(size0);
	Eigen::VectorXd results2(size0);
	for(int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenBiCGstab(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}


/*//=======================================================
// ● ICCGで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init){
	/* 設定 */
	Eigen::ConjugateGradient<Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<dcomplex, Eigen::Lower, Eigen::NaturalOrdering<int>>> iccg;
	iccg.setMaxIterations(max_ite);
	iccg.setTolerance(conv_cri);

	iccg.compute(matA.matrix.matrix);
	if(init){
		results = iccg.solve(vecB);
	}else{
		results = iccg.solveWithGuess(vecB, results);
	}
	double err = iccg.error();
	//std::cout << "err " << err << std::endl;
	return(err <= conv_cri);
}
bool MatSolversEigenMKL::solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init){
	Eigen::VectorXcd vecB2(size0);
	Eigen::VectorXcd results2(size0);
	for(int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenICCG(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}


/*//=======================================================
// ● BiCGstabで解く(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init){
	/* 設定 */
	Eigen::BiCGSTAB<Eigen::SparseMatrix<dcomplex, Eigen::RowMajor>, Eigen::IncompleteLUT<dcomplex>> bicg;
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
bool MatSolversEigenMKL::solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init){
	Eigen::VectorXcd vecB2(size0);
	Eigen::VectorXcd results2(size0);
	for(int i = 0 ; i < size0 ; i++){
		vecB2(i) = vecB[i];
	}
	bool bl = solveEigenBiCGstab(size0, conv_cri, max_ite, matA, vecB2, results2, init);
	for(int i = 0 ; i < size0 ; i++){
		results[i] = results2(i);
	}
	return bl;
}


/*//=======================================================
// ● SimplicialLDLT(Eigenソルバ)
//=======================================================*/
bool MatSolversEigenMKL::solveEigenSimplicialLDLT(const slv_int size0, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results){	
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
	solver.compute(matA.matrix.matrix);
	if(solver.info() != Eigen::Success) {
		return false;
	}
	results = solver.solve(vecB);
	if(solver.info() != Eigen::Success) {
		return false;;
	}
	return true;
}
bool MatSolversEigenMKL::solveEigenSimplicialLDLT(const slv_int size0, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results){	
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<dcomplex> > solver;
	solver.compute(matA.matrix.matrix);
	if(solver.info() != Eigen::Success) {
		return false;
	}
	results = solver.solve(vecB);
	if(solver.info() != Eigen::Success) {
		return false;;
	}
	return true;
}

/* end of namespace */
};



