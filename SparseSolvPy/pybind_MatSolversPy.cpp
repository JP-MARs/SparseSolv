#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* 専用名前空間 */
namespace SRLfem{


/*//=======================================================
// ● 残差履歴取得forPython
//=======================================================*/
std::vector<double> MatSolvers::getResidualLog_py(){
	return residual_log;
}

/*//=======================================================
// ● ICCGで解く
//=======================================================*/
bool MatSolvers::solveICCG_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, py::array_t<double> vecB, py::array_t<double> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    double *b_ptr = (double *)b_buf.ptr;
	auto x_buf = results.request();
    double *x_ptr = (double *)x_buf.ptr;
	/* */

	double* vecBa = new double[size0];
	double* results_a = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const double temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
		std::cout << temp << std::endl;
	}
	norm = sqrt(norm);
	bool bl = solveICCG(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results_a, init);

	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}


/*//=======================================================
// ● ICCGで解く
//=======================================================*/
bool MatSolvers::solveICCG_py(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    std::complex<double> *b_ptr = (std::complex<double> *)b_buf.ptr;
	auto x_buf = results.request();
    std::complex<double> *x_ptr = (std::complex<double> *)x_buf.ptr;
	/* */
	
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex* results_a = new dcomplex[size0];
	dcomplex norm = 0;
	for(int i = 0 ; i < size0 ; i++){
		const dcomplex temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
	}
	double normB2 = abs(norm);
	normB2 = sqrt(normB2);
	bool bl = solveICCG(size0, conv_cri, max_ite, accera, normB2, matA, vecBa, results_a, init);

	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}

/*//=======================================================
// ● ABMC- ICCGで解く
//=======================================================*/
bool MatSolvers::solveICCGwithABMC_py(const slv_int size0, const double conv_cri,  const int max_ite, const double accera, const SparseMat &matA, 
						   py::array_t<double> vec_b, py::array_t<double> results, int num_blocks, int num_colors, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vec_b.request();
    double *b_ptr = (double *)b_buf.ptr;
	auto x_buf = results.request();
    double *x_ptr = (double *)x_buf.ptr;
	/* */

    double* vecBa = new double[size0];
    double* results_a = new double[size0];
    double norm=0;
    for(int i = 0 ; i < size0 ; i++){
        const double temp = b_ptr[i];
        vecBa[i] = temp;
        norm += temp*temp;
        results_a[i] = x_ptr[i];
    }
    norm = sqrt(norm);
    bool bl = solveICCGwithABMC(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results_a, num_blocks, num_colors, init);
    for(int i = 0 ; i < size0 ; i++){
        x_ptr[i] = results_a[i];
    }
    delete[] vecBa;
    delete[] results_a;
    return bl;
}

/*//=======================================================
// ● IC-MRTRで解く
//=======================================================*/
bool MatSolvers::solveICMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, py::array_t<double> vecB, py::array_t<double> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    double *b_ptr = (double *)b_buf.ptr;
	auto x_buf = results.request();
    double *x_ptr = (double *)x_buf.ptr;
	/* */

	double* vecBa = new double[size0];
	double* results_a = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const double temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
	}
	norm = sqrt(norm);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results_a, init);
	
	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}
/*//=======================================================
// ● IC-MRTRで解く
//=======================================================*/
bool MatSolvers::solveICMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    std::complex<double> *b_ptr = (std::complex<double> *)b_buf.ptr;
	auto x_buf = results.request();
    std::complex<double> *x_ptr = (std::complex<double> *)x_buf.ptr;
	/* */
	dcomplex* vecBa = new dcomplex[size0];
	dcomplex* results_a = new dcomplex[size0];
	dcomplex norm = 0;
	for(int i = 0 ; i < size0 ; i++){
		const dcomplex temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
	}
	double normB2 = abs(norm);
	normB2 = sqrt(normB2);
	bool bl = solveICMRTR(size0, conv_cri, max_ite, accera, normB2, matA, vecBa, results_a, init);
	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] results_a;
	delete[] vecBa;
	return bl;
}

/*//=======================================================
// ● SGS-MRTRで解く
//=======================================================*/
bool MatSolvers::solveSGSMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, py::array_t<double> vecB, py::array_t<double> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    double *b_ptr = (double *)b_buf.ptr;
	auto x_buf = results.request();
    double *x_ptr = (double *)x_buf.ptr;
	/* */

	double* vecBa = new double[size0];
	double* results_a = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const double temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
	}
	norm = sqrt(norm);
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm, matA, vecBa, results_a, init);
	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}
/*//=======================================================
// ● SGS-MRTRで解く
//=======================================================*/
bool MatSolvers::solveSGSMRTR_py(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, py::array_t<std::complex<double>> vecB, py::array_t<std::complex<double>> results, bool init){
	/* Python側へのポインタ設定 */
	auto b_buf = vecB.request();
    std::complex<double> *b_ptr = (std::complex<double> *)b_buf.ptr;
	auto x_buf = results.request();
    std::complex<double> *x_ptr = (std::complex<double> *)x_buf.ptr;
	/* */

	dcomplex* vecBa = new dcomplex[size0];
	dcomplex* results_a = new dcomplex[size0];
	dcomplex norm=0;
	for(int i = 0 ; i < size0 ; i++){
		const dcomplex temp = b_ptr[i];
		vecBa[i] = temp;
		norm += temp*temp;
		results_a[i] = x_ptr[i];
	}
	double norm0 = abs(norm);
	double norm00 = sqrt(norm0);
	bool bl = solveSGSMRTR(size0, conv_cri, max_ite, norm00, matA, vecBa, results_a, init);
	for(int i = 0 ; i < size0 ; i++){
		x_ptr[i] = results_a[i];
	}
	delete[] vecBa;
	delete[] results_a;
	return bl;
}


/* end of namespace */
}



