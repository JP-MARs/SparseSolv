#define _CRT_SECURE_NO_WARNINGS

/* Eigen自体の並列化の無効化 */
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;

#include <SparseSolve/SparseMatTMPL.hpp>
#include <SparseSolve/SparseMatC.hpp>
#include <SparseSolve/MatSolvers.hpp>
#include <SparseSolve/MatSolversEigenMKL.hpp>

#include <chrono>

// タイマー
double timeit(std::function<void()> f) {
    auto start = chrono::high_resolution_clock::now();
    f();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> diff = end - start;
    return diff.count();
}

/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char* argv[]){
	double matA[5][5] ={
		{ 61,    0,   3.1,  -0.1,   1},
		{  0,   42.5, 2.0,   0.5,   4},
		{ 3.1,  2.0, 10.5,     0,   0},
		{-0.1, 0.5,     0,  80.4, 1.5},
		{   1,   4,     0,   1.5,  1000.0}
	};
	const int total_size = 5;


	/* 疎行列を定義 */
	SRLfem::SparseMat matAs1(total_size);
	/* 行列A,B,Cの値をセット */
	for(int i = 0 ; i < 5 ; i++){
		for(int j = 0 ; j < 5 ; j++){
			if( fabs(matA[i][j]) > 1.0e-12 ){
				matAs1.add(i, j, matA[i][j]);
			}
		}
	}
	/* 疎行列位置を確定させる */
	matAs1.refresh();


	/* 適当なベクトル */
	double vecB[] = {1, 2, 3, 4, 5};
	SRLfem::dcomplex* vec0C = new SRLfem::dcomplex[total_size];
	for(int i = 0 ; i < total_size ; i++){
		vec0C[i] = i;//rand() / (1.0+RAND_MAX);
	}

	/**/
	/**/
	/* ソルバ */
	/**/
	/**/

	/* 解を初期化 */
	double* results00 = new double[total_size];
	for(int i = 0 ; i < total_size ; i++){
		results00[i] = 0;
	}

	cout <<"start" << endl;
	SRLfem::MatSolvers solvers;
	//solvers.setSaveLog(true);
	//solvers.setDiagScale(true);
	/* ICCG */
	//bool bl1 = solvers.solveSGSMRTR(total_size, 1.0e-8, 10000, matAs1, vecB, results00);
	//bool bl1 = solvers.solveICCGwithABMC(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00, 5, 5);

	vector<double> vecB_vec(vecB, vecB+5);
	vector<double> results00_vec;
	results00_vec.resize(total_size, 0.0);
	//bool bl1 = solvers.solveICCG(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);
	bool bl1 = solvers.solveICCG(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB_vec, results00_vec);
	// 
	//bool bl1 = solvers.solveICMRTR(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);

	vector<double> log1;
	solvers.getResidualLog(log1);

	/* ICCOCG */
	//bool bl2 = solvers.solveSGSMRTR(total_size, 1.0e-8, 10000, -1.05, matCs1, vec0C, results01);

	for(auto itr : log1){
		cout << itr << endl;
	}


	cout << "D norm" << endl;
	/*
	double* norm1 = matAs1 * results00;
	for(int i = 0 ; i < total_size ; i++){
		cout << norm1[i] << ", " << results00[i] << endl;
	}
	delete[] norm1;
	*/
	double* norm1 = matAs1 * results00_vec.data();
	for(int i = 0 ; i < total_size ; i++){
		cout << norm1[i] << ", " << results00_vec[i] << endl;
	}
	delete[] norm1;


	cout << endl << "BiCGstab" << endl;
	results00_vec.clear();
	results00_vec.resize(total_size, 0.0);
	bl1 = SRLfem::MatSolversEigenMKL::solveEigenBiCGstab_accel(total_size, 1.0e-8, 10000, matAs1, vecB_vec, results00_vec);
	for(int i = 0 ; i < total_size ; i++){
		cout << results00_vec[i] << endl;
	}

	cout << endl << "GMRES" << endl;
	results00_vec.clear();
	results00_vec.resize(total_size, 0.0);
	bl1 = SRLfem::MatSolversEigenMKL::solveEigenGMRES(total_size, 1.0e-8, 10000, matAs1, vecB_vec, results00_vec);
	for(int i = 0 ; i < total_size ; i++){
		cout << results00_vec[i] << endl;
	}


	getchar();

	return 1;
}


#ifdef AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char* argv[]){
	cout << sizeof(vector<map<int, double>>) << endl;
	cout << sizeof(Eigen::SparseMatrix<double>) << endl;
	cout << sizeof(SRLfem::SparseMatBaseD) << endl;
	cout << sizeof(SRLfem::SparseMat) << endl;

	/* Eigenの並列化を無効化 */
	Eigen::initParallel();

	Eigen::setNbThreads(1);
	int n = Eigen::nbThreads( );
	cout << "Thred Eigen " << n << endl;

	/*size_t test_size = 10;
	for(size_t i = test_size; i-- > 0;) {
		cout << "test : " << i << endl;
	}
	return 1;*/



#ifdef AAAAAAAAAAAA

	//EIGEN_NO_DEBUG 



	int dim, data_size;
	fstream fp1("./mat.csv", std::ios::in);
	fp1 >> dim >> data_size;
	cout << dim << ", " << data_size << endl;
	SRLfem::SparseMatC matA(dim);
	for(int i = 0; i < data_size; i++){
		int gyo, retu;
		double val_r, val_i;
		fp1 >> gyo >> retu >> val_r >> val_i;
		SRLfem::dcomplex val(val_r, val_i);
		matA.add(gyo, retu, val);
	}
	matA.fix();
	fp1.close();

	SRLfem::dcomplex *vec = new SRLfem::dcomplex[dim];
	Eigen::VectorXcd vecE = Eigen::VectorXcd::Zero(dim);

	fstream fp2("./vec.csv", std::ios::in);
	fp2 >> dim;
	for(int i = 0; i < dim; i++){
		double val_r, val_i;
		fp2 >> val_r >> val_i;
		SRLfem::dcomplex val(val_r, val_i);
		vec[i] = val;
		vecE(i) = val;
	}	
	fp2.close();


	cout <<"start" << endl;
	/* 解を初期化 */
	SRLfem::dcomplex* results01 = new SRLfem::dcomplex[dim];
	Eigen::VectorXcd resultE = Eigen::VectorXcd::Zero(dim);
	for(int i = 0 ; i < dim ; i++){
		results01[i] = 0;
	}

	SRLfem::MatSolvers solvers;
	solvers.setSaveLog(true);
	//solvers.setDiagScale(true);
	/* ICCG */
	//bool bl1 = solvers.solveSGSMRTR(dim, 1.0e-6, 10000, matA, vec, results01);
	//bool bl1 = solvers.solveICCG(dim, 1.0e-6, 1000, 1.2, matA, vec, results01);
	bool bl1 = solvers.solveICMRTR(dim, 1.0e-6, 1000, 1.2, matA, vec, results01);
	bool bl2 = false;//solvers.solveSGSMRTR(dim, 1.0e-6, 1000, matA, vec, results01);
	
	//bool bl2 = false;//SRLfem::MatSolversEigenMKL::solveEigenICCG(dim, 1.0e-6, 1000, matA, vecE, resultE);

	std::vector<double> log;
	fp1.open("log.csv", std::ios::out);
	solvers.getResidualLog(log);
	int c=1;
	for(auto val : log){
		cout << val << endl;
		fp1 << c << ", " << val << endl;
		c++;
	}
	fp1.close();

	cout << "conv = " << bl1 << ", " << bl2 << endl;
	getchar();
	return 1;

#endif





	/*
	const int test_size = 1500;
	SRLfem::SparseMatBaseD mat1(test_size);
	SRLfem::SparseMatBaseD mat2(test_size);
	mat2.fix();

	auto t_start = std::chrono::system_clock::now();
	for(int i = 0; i < test_size; i++) {
		for(int j = 0; j < test_size; j++) {
			mat1.add(i, j, 1);
		}
	}
	mat1.fix();
	cout << "end1"<<endl;
	auto t_end = std::chrono::system_clock::now();

	for(int i = 0; i < test_size; i++) {
		for(int j = 0; j < test_size; j++) {
			mat2.add(i, j, 1);
		}
	}
	auto t_end2 = std::chrono::system_clock::now();
	cout << "end2"<<endl;

	auto t_dur = t_end - t_start;
	auto t_dur2 = t_end2 - t_end;
	auto t_sec = std::chrono::duration_cast<std::chrono::seconds>(t_dur).count();
	auto t_sec2 = std::chrono::duration_cast<std::chrono::seconds>(t_dur2).count();
	std::cout << "Original insert  time: " << t_sec << " sec \n";
	std::cout << "Eigen insert time: " << t_sec2 << " sec \n";
	*/

	double matA[5][5] ={
		{ 61,    0,   3.1,  -0.1,   1},
		{  0,   42.5, 2.0,   0.5,   4},
		{ 3.1,  2.0, 10.5,     0,   0},
		{-0.1, 0.5,     0,  80.4, 1.5},
		{   1,   4,     0,   1.5,  1000.0}
	};
	/* 複素行列C */
	SRLfem::dcomplex matC[5][5] ={
		{15.5, 2, 3i, 0, 1.5},
		{2, 10, 0, 0, 0},
		{3i, 0, 8, 0, 1},
		{0, 0, 0, 4.8, 0.5i},
		{1.5, 0, 1, 0.5i, 20.5}
	};

	const int total_size = 5;


	/* 疎行列を定義 */
	SRLfem::SparseMat matAs1(total_size);
	SRLfem::SparseMat matBs1(total_size);
	SRLfem::SparseMatC matCs1(total_size);
	/* 行列A,B,Cの値をセット */
	for(int i = 0 ; i < 5 ; i++){
		for(int j = 0 ; j < 5 ; j++){
			if( fabs(matA[i][j]) > 1.0e-12 ){
				matAs1.add(i, j, matA[i][j]);
				matBs1.add(i, j, matA[i][j]);
			}
			if( fabs(abs(matC[i][j])) > 1.0e-12 ){
				matCs1.add(i, j, matC[i][j]);
			}
		}
	}
	/* 疎行列位置を確定させる */
	matAs1.refresh();
	matBs1.refresh();
	matCs1.refresh();

	SRLfem::SparseMatC matDot = matAs1*matBs1*matCs1;
	matDot.print();

	double time_eigen = timeit([&]() {
		for(int i = 0; i < 1000000; i++){
			matDot = matAs1*matBs1*matCs1;
			//cout << "--------------" << endl;
			//matDot.print();
		}
	});

	double time_manual = timeit([&]() {
		for(int i = 0; i < 1000000; i++){
			SRLfem::SparseMatOperators::dotFix(matDot, matAs1, matBs1, matCs1);
			//SRLfem::SparseMatC matDot0 = SRLfem::SparseMatOperators::dotMats(matAs1, matBs1, matCs1);
			//cout << "--------------" << endl;
			//matDot.print();
		}
	});

    std::cout << "Eigen product time: " << time_eigen << " sec" << std::endl;
    std::cout << "Manual product time: " << time_manual << " sec" << std::endl;

	return 1;


	Eigen::VectorXd vec0 = Eigen::VectorXd::Zero(5);
	//while(true){
		Eigen::VectorXd vec1 = matAs1 * vec0;
		SRLfem::SparseMat testX;
		testX = matAs1 * matBs1;
	//}

	matAs1.printMat("./testMat.csv");


	/*cout << "cols " << matAs1.getMaxCol() << endl;
	SRLfem::SparseMat tempMat0 = matAs1.makePrsdInv(1.0e-1);
	tempMat0.printMat();

	SRLfem::SparseMat tempMat1;
	string str1 = "Mat.csv";
	tempMat1.readMat(str1);
	tempMat1.printMat("test2.csv");*/


	SRLfem::SparseMat tempMat01 = matAs1 + matAs1;
	SRLfem::SparseMatOperators::plusFix(tempMat01, matAs1, matAs1);
	//tempMat01.printMat();
	//getchar();


	double diagD[100];
	SRLfem::SparseMat matLL = matAs1.IC_decomp(diagD, 1.0);
	//matLL.print();


	/*double vec[5] ={1,2,3,4,5};

	SRLfem::dcomplex* vec2 = matC*vec;
	for(int i = 0 ; i < 5 ; i++){
		cout << vec2[i] <<endl;
	}
	getchar();*/






	/* 適当なベクトル */
	double vecB[] = {1, 2, 3, 4, 5};
	SRLfem::dcomplex* vec0C = new SRLfem::dcomplex[total_size];
	for(int i = 0 ; i < total_size ; i++){
		vec0C[i] = i;//rand() / (1.0+RAND_MAX);
	}

	/* 行列・ベクトル積 */
	/*while(true){
		SRLfem::SparseMat tempMat01;
		SRLfem::SparseMat tempMat02;
		matAs1.MatDiv(tempMat01, tempMat02, 1, 2);
		SRLfem::SparseMat tempMat03 = matAs1.getMatUpper();
	}*/



	/**/
	/**/
	/* ソルバ */
	/**/
	/**/

	/* 解を初期化 */
	double* results00 = new double[total_size];
	for(int i = 0 ; i < total_size ; i++){
		results00[i] = 0;
	}
	SRLfem::dcomplex* results01 = new SRLfem::dcomplex[total_size];
	for(int i = 0 ; i < total_size ; i++){
		results01[i] = 0;
	}

	cout <<"start" << endl;
	SRLfem::MatSolvers solvers;
	solvers.setSaveLog(true);
	solvers.setDiagScale(true);
	/* ICCG */
	//bool bl1 = solvers.solveSGSMRTR(total_size, 1.0e-8, 10000, matAs1, vecB, results00);
	//bool bl1 = solvers.solveICCGwithABMC(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00, 5, 5);
	//bool bl1 = solvers.solveICCG(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);
	//bool bl1 = solvers.solveICMRTR(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);
	SRLfem::MatSolversEigenMKL solver2;
	bool bl1 = solver2.solveEigenBiCGstab_accel(total_size, 1.0e-8, 10000, matAs1, vecB, results00);

	vector<double> log1;
	solvers.getResidualLog(log1);

	/* ICCOCG */
	//bool bl2 = solvers.solveSGSMRTR(total_size, 1.0e-8, 10000, -1.05, matCs1, vec0C, results01);

	for(auto itr : log1){
		cout << itr << endl;
	}


	cout << "D norm" << endl;
	double* norm1 = matAs1 * results00;
	for(int i = 0 ; i < total_size ; i++){
		cout << norm1[i] << ", " << results00[i] << endl;
	}
	cout << "C norm" << endl;
	SRLfem::dcomplex* norm2 = matCs1 * results01;
	for(int i = 0 ; i < total_size ; i++){
		cout << norm2[i] << ", " << results01[i] << endl;
	}

	delete[] norm1;
	delete[] norm2;

	return 1;



	Eigen::VectorXd vecBe(total_size);
	for(int i = 0 ; i < total_size ; i++){
		vecBe(i) = vecB[i];
	}
	Eigen::VectorXd results_e(total_size);

	matAs1.print();
	cout << "---------" << endl;
	cout << "start...." << endl;

	bool bl = SRLfem::MatSolversEigenMKL::solveEigenICCG(total_size, 1.0e-8, 10000, matAs1, vecBe, results_e, true);
	cout << results_e << endl;
	cout << "D norm3" << endl;
	Eigen::VectorXd norm1b = matAs1 * results_e;
	cout << norm1b <<endl;
	bool bl2x = SRLfem::MatSolversEigenMKL::solveEigenICCG(total_size, 1.0e-8, 10000, matAs1, vecBe, results_e, false);
	cout << results_e << endl;
	cout << "D norm4" << endl;
	norm1b = matAs1 * results_e;
	cout << norm1b <<endl;

	
	//matAs1.add(4, 2, 11);
	/*vector<int> aa;
	vector<double> bb;
	matAs1.getTargetColVal(0, aa, bb);
	cout << "results, "  << endl;
	for(auto itr : aa) {
		cout << itr << endl;
	}
	matAs1.print();


	SRLfem::SparseMatBaseD matB2 = matAs1;
	matB2 *= 1.5;
	matB2.print();

	cout << &matAs1 << ", " << &matB2 << endl;*/

	return 1;

}

#endif
