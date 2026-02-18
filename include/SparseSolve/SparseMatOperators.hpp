/**
 * @file SparseMatOperators.hpp
 * @brief Operator define for Sparse Matrix
 * 
 * このファイルは SparseMatTMPL への演算を提供します。
 */
#ifndef DEF_SPR_MAT_OPERATORS
#define DEF_SPR_MAT_OPERATORS

#include "SparseMatTMPL.hpp"
#include <000_thirdparty/Eigen/Dense>
#include <cfloat>

/* 行列積などをOpenMP並列化するとき、onしてください */
//#define OMP_USING_MAT_SOL


/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列　操作staticクラス
//=======================================================*/
/**
 * @class SparseMatOperators
 * @brief Operator static class for Sparse matrices
 * このファイルは SparseMatTMPL への演算を提供するstaticクラスです。
 */
class SparseMatOperators {
friend class SparseMat;
friend class SparseMatC;
friend class MatSolvers;
private:
	/* 行列ベクトル積 */
	template<typename Mat1, typename DType1, typename DType2, typename DType3> 
	static void mat_vec_product(DType3** ans, const Mat1& matA, const DType2* const vec);
	/* 行列ベクトル積(std::vectorとの) */
	template<typename Mat1, typename DType3> 
	static void mat_vec_product(std::vector<DType3>& ans, const Mat1& matA, const std::vector<double>& vec);
	template<typename Mat1> 
	static void mat_vec_product(std::vector<dcomplex>& ans, const Mat1& matA, const std::vector<dcomplex>& vec);
	/* 行列ベクトル積(std::vectorと。結果はEigen) */
	template<typename Mat1, typename DType3> 
	static void mat_vec_product(DType3& ans, const Mat1& matA, const std::vector<double>& vec);
	template<typename Mat1> 
	static void mat_vec_product(Eigen::VectorXcd& ans, const Mat1& matA, const std::vector<dcomplex>& vec);
	/* 足し算オペレータ(matAにmatBを加える(a1*A+a2*B)。Bの位置を行方向にpos1、列方向にpos2だけずらす。) */
	template<typename Mat1, typename Mat2, typename Mat3, typename DType1, typename DType2> 
	static void plus_operators(Mat3& mat3, const Mat1& matA, const Mat2& matB, double a1, double a2, slv_int pos1, slv_int pos2);
	/* matAにa1*matBを加える。 開始位置はposだけずらす */
	template<typename Mat1, typename Mat2, typename DType1, typename DType2, typename DType3> 
	static void MatSelfPlus(Mat1& matA, const Mat2& matB, const DType3 a1, const slv_int pos1, const slv_int pos2);
	/**/
	/* その他 */
	/**/
	/* 実数→複素化 */
	static void copyDtoC(SparseMatBaseC& mat_ans, const SparseMatBaseD& mat);
	/* 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、Bに渡す */
	template<typename Mat1, typename DType>
	static void makeSubMat(Mat1& mat_ans, const Mat1& targetMat, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b);
	/* 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する */
	template<typename Mat1, typename DType>
	static void MatDiv(const Mat1& targetMat, Mat1& matK11, Mat1& matK12, slv_int rangeA, slv_int rangeB);	
	/* 下三角行列を取得 */
	template<typename Mat1, typename DType>
	static void getMatLower(Mat1& ans_mat, const Mat1& targetMat);
	/* 上三角行列を取得 */
	template<typename Mat1, typename DType>	
	static void getMatUpper(Mat1& ans_mat, const Mat1& targetMat);
	/* 行列Aの逆行列を求めてInvに渡す(密用) */
	template<typename Mat1, typename DType, typename MType>	
	static void MatInv(Mat1& ans_mat, const Mat1& targetMat);
	/**/
	/* 不完全コレスキー分解 */
	template<typename Mat1, typename DType1, typename DTypeD>
	static void IC_decomp(Mat1& mat_ans, const Mat1& mat1, DTypeD* diagD, const double accela);
	/**/
	/* 対角スケーリング */
	static void diagScaling(SparseMatBaseD& mat_ans, const SparseMatBaseD& mat1, double* trans_vec, const double* ori_vec);
	/**/
	/* 対角スケーリング(複素) */
	static void diagScalingComplex(SparseMatBaseC& mat_ans, const SparseMatBaseC& mat1, dcomplex* trans_vec, const dcomplex* ori_vec);
	/**/
	/* A^T*A + epsIを作る（疑似逆行列用） */
	template<typename Mat1, typename DType1>
	static void AtA_eps(Mat1& mat_ans, const Mat1& mat1, double eps);
	/**/
	/*===========================================*/
	/* Staticで使う専用の内部処理 */
	/*===========================================*/
	/**/
	/* (mat1 + mat2)*vecBを計算 */
	template<typename Mat1, typename Mat2, typename DType1, typename DType2, typename DType0, typename DTypeA> 
	static void productVecMat2(DTypeA** ans, const Mat1& mat1, const Mat2& mat2, const DType0* vecB);
	/* 行列AとBとCをかけて結果を返す*/
	template<typename Mat1, typename Mat2, typename Mat3, typename Mat0, typename DType1, typename DType2, typename DType3, typename DType0> 
	static void MatProduct(Mat0& mat_ans, const Mat1& matA, const Mat2& matB, const Mat3& matC);		
	/**/
	/**/
	/* nonzero既知固定用 */
	/**/
	/* (配列位置確定時)自身に、行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす */
	template<typename Mat1, typename Mat2, typename Mat0, typename DType1, typename DType2, typename DType0> 
	static void PlusMinusShiftFix(Mat0& matAB, const Mat1& matA, const Mat2& matB, double a1, double a2, slv_int pos1, slv_int pos2);	
	/* (配列位置確定時)自身に、行列A±B±Cを加える(a1*A+a2*B+a3*C)。 */
	template<typename Mat1, typename Mat2, typename Mat3, typename Mat0, typename DType1, typename DType2, typename DType3, typename DType0> 
	static void PlusMinusShiftFix(Mat0& matABC, const Mat1& matA, const Mat2& matB, const Mat3& matC, double a1, double a2, double a3);	
	/* (配列位置確定時)行列AとBをかけて、自身ABに格納する */
	template<typename Mat1, typename Mat2, typename Mat0, typename DType1, typename DType2, typename DType0> 
	static void MatProductFix(Mat0& matAB, const Mat1& matA, const Mat2& matB);
	/**/
public:
	/* (a1*mat1 + a2*mat2)で、mat2に足す位置をpos1,pos2ずらすを計算 */
	static SparseMat plusShift(const SparseMat& mat1, const SparseMat& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2);
	static SparseMatC plusShift(const SparseMatC& mat1, const SparseMatC& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2);
	/* (mat1 + mat2)*vecBを計算 */
	static double* dotVecMat2(const SparseMat& mat1, const SparseMat& mat2, const double* vec1);
	static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMatC& mat2, const double* vec1);
	static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMat& mat2, const double* vec1);
	static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMatC& mat2, const double* vec1);
	static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMat& mat2, const dcomplex* vec1);
	static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMatC& mat2, const dcomplex* vec1);
	static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMat& mat2, const dcomplex* vec1);
	static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMatC& mat2, const dcomplex* vec1);
	/* 行列AとBとCをかけて結果を返す*/
	static SparseMat dotMats(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC);
	static SparseMatC dotMats(const SparseMat& matA, const SparseMat& matB, const SparseMatC& matC);
	static SparseMatC dotMats(const SparseMat& matA, const SparseMatC& matB, const SparseMatC& matC);
	static SparseMatC dotMats(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC);
	/**/
	/* (配列位置確定時)自身に、行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす */
	static void plusFix(SparseMat& matAB, const SparseMat& matA, const SparseMat& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
	static void plusFix(SparseMatC& matAB, const SparseMat& matA, const SparseMatC& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
	static void plusFix(SparseMatC& matAB, const SparseMatC& matA, const SparseMatC& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
	/* (配列位置確定時)自身に、行列A±B±Cを加える(a1*A+a2*B+a3*C) */
	static void plusFix(SparseMat& matABC, const SparseMat& matA, const SparseMat& matB, const SparseMat& matC, double a1=1.0, double a2=1.0, double a3=1.0);
	static void plusFix(SparseMatC& matABC, const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC, double a1=1.0, double a2=1.0, double a3=1.0);

	/* **こいつらは自作掛け算なので、たぶん普通に上で定義している掛け算を使ったほうが速い */
	/* (配列位置確定時)行列AとBをかけて、自身ABに格納する */
	static void dotFix(SparseMat& matAB, const SparseMat& matA, const SparseMat& matB);
	static void dotFix(SparseMatC& matAB, const SparseMat& matA, const SparseMatC& matB);
	static void dotFix(SparseMatC& matAB, const SparseMatC& matA, const SparseMatC& matB);
	/* (配列位置確定時)行列AとBとCをかけて、自身ABCに格納する */
	static void dotFix(SparseMat& matABC, const SparseMat& matA, const SparseMat& matB, const SparseMat& matC);
	static void dotFix(SparseMatC& matABC, const SparseMat& matA, const SparseMat& matB, const SparseMatC& matC);
	static void dotFix(SparseMatC& matABC, const SparseMat& matA, const SparseMatC& matB, const SparseMatC& matC);
	static void dotFix(SparseMatC& matABC, const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC);
};


/*//=======================================================
// ● 行列ベクトル積
//=======================================================*/
template<typename Mat1, typename DType1, typename DType2, typename DType3> 
void SparseMatOperators::mat_vec_product(DType3** ans, const Mat1& matA, const DType2* const vec){	
	const slv_int the_size = matA.size;
	auto col_ptr = matA.getColPtr();
	auto val_ptr = matA.getValuePtr();
	/* 列の開始と終わりのPtrを取得 */
	slv_int* start_pos = new slv_int[the_size];
	slv_int* end_pos = new slv_int[the_size];
	matA.getCols(start_pos, end_pos);		
	/* 掛け算開始 */
	DType3* new_vec = new DType3[the_size];
	for(slv_int i = 0 ; i < the_size ; i++){
		DType3 temp = 0;
		const slv_int ss = end_pos[i];
		for(slv_int j = start_pos[i] ; j < ss ; j++){
			const slv_int pos = col_ptr[j];
			temp += vec[pos]*val_ptr[j];
		}
		new_vec[i] = temp;
	}
	*ans = std::move(new_vec);
	new_vec = nullptr;
	delete[] start_pos;
	delete[] end_pos;
}

/*//=======================================================
// ● 行列ベクトル積(std::vectorとの)
//=======================================================*/
template<typename Mat1, typename DType3> 
void SparseMatOperators::mat_vec_product(std::vector<DType3>& ans, const Mat1& matA, const std::vector<double>& vec){
	const slv_int the_size = matA.size;
	ans.resize(the_size);
	auto col_ptr = matA.getColPtr();
	auto val_ptr = matA.getValuePtr();
	/* 列の開始と終わりのPtrを取得 */
	slv_int* start_pos = new slv_int[the_size];
	slv_int* end_pos = new slv_int[the_size];
	matA.getCols(start_pos, end_pos);		
	/* 掛け算開始 */
	for(slv_int i = 0 ; i < the_size ; i++){
		DType3 temp = 0;
		const slv_int ss = end_pos[i];
		for(slv_int j = start_pos[i] ; j < ss ; j++){
			const slv_int pos = col_ptr[j];
			temp += vec[pos]*val_ptr[j];
		}
		ans[i] = temp;
	}
	delete[] start_pos;
	delete[] end_pos;

	/*const slv_int the_size = vec.size();
	Eigen::VectorXd vecE(the_size);
	for(slv_int i = 0; i < the_size; i++){
		vecE(i) = vec[i];
	}*/
	//Eigen::VectorXd vec_new = matA.matrix*vecE;	
	/*ans.resize(the_size);
	for(slv_int i = 0; i < the_size; i++){
		ans[i] = vec_new(i);
	}*/
}
template<typename Mat1> 
void SparseMatOperators::mat_vec_product(std::vector<dcomplex>& ans, const Mat1& matA, const std::vector<dcomplex>& vec){
	const slv_int the_size = matA.size;
	ans.resize(the_size);
	auto col_ptr = matA.getColPtr();
	auto val_ptr = matA.getValuePtr();
	/* 列の開始と終わりのPtrを取得 */
	slv_int* start_pos = new slv_int[the_size];
	slv_int* end_pos = new slv_int[the_size];
	matA.getCols(start_pos, end_pos);		
	/* 掛け算開始 */
	for(slv_int i = 0 ; i < the_size ; i++){
		dcomplex temp = 0;
		const slv_int ss = end_pos[i];
		for(slv_int j = start_pos[i] ; j < ss ; j++){
			const slv_int pos = col_ptr[j];
			temp += vec[pos]*val_ptr[j];
		}
		ans[i] = temp;
	}
	delete[] start_pos;
	delete[] end_pos;
	
	/*const slv_int the_size = vec.size();
	Eigen::VectorXcd vecE(the_size);
	for(slv_int i = 0; i < the_size; i++){
		vecE(i) = vec[i];
	}
	Eigen::VectorXcd vec_new = matA.matrix*vecE;
	ans.resize(the_size);
	for(slv_int i = 0; i < the_size; i++){
		ans[i] = vec_new(i);
	}*/
}


/*//=======================================================
// ● 行列ベクトル積(std::vectorと。結果はEigen)
//=======================================================*/
template<typename Mat1, typename DType3> 
void SparseMatOperators::mat_vec_product(DType3& ans, const Mat1& matA, const std::vector<double>& vec){
	const slv_int the_size = vec.size();
	Eigen::VectorXd vecE(the_size);
	for(slv_int i = 0; i < the_size; i++){
		vecE(i) = vec[i];
	}
	ans = matA.matrix*vecE;
}

template<typename Mat1> 
void SparseMatOperators::mat_vec_product(Eigen::VectorXcd& ans, const Mat1& matA, const std::vector<dcomplex>& vec){
	const slv_int the_size = vec.size();
	Eigen::VectorXcd vecE(the_size);
	for(slv_int i = 0; i < the_size; i++){
		vecE(i) = vec[i];
	}
	ans = matA.matrix*vecE;
}


/*//=======================================================
// ● 足し算オペレータ(matAにmatBを加える(a1*A+a2*B)。Bの位置を行方向にpos1、列方向にpos2だけずらす。)
//=======================================================*/
template<typename Mat1, typename Mat2, typename Mat3, typename DType1, typename DType2> 
void SparseMatOperators::plus_operators(Mat3& mat3, const Mat1& matA, const Mat2& matB, double a1, double a2, slv_int pos1, slv_int pos2){
	/* 最大の行数を計算 */
	const slv_int size1 = matA.size;
	const slv_int size2 = matB.size;
	const slv_int max = (size1 > size2+pos1 ? size1 : size2+pos1);

	/* 足した後の行列 */
	Mat3 matAB(max);

	/* まずAを足す */
	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr1[aj];
			DType1 avalue = val_ptr1[aj];
			avalue *= a1;
			matAB.add(i, col, avalue);
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;

	/* Bをずらしながら足す */
	slv_int* start_pos2 = new slv_int[size2];
	slv_int* end_pos2 = new slv_int[size2];
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr2[aj];
			DType2 avalue = val_ptr2[aj];
			avalue *= a2;
			matAB.add(i+pos1, col+pos2, avalue);
		}
	}
	delete[] start_pos2;
	delete[] end_pos2;

	/* 確定させ、ムーブする */
	matAB.fix();
	mat3 = std::move(matAB);
}

/*//=======================================================
// ● matAにa1*matBを加える。 開始位置はposだけずらす
//=======================================================*/
template<typename Mat1, typename Mat2, typename DType1, typename DType2, typename DType3> 
void SparseMatOperators::MatSelfPlus(Mat1& matA, const Mat2& matB, const DType3 a1, const slv_int pos1, const slv_int pos2){
	if(matA.is_fix){
		const slv_int max_size = (matA.size > matB.size+pos1 ? matA.size : matB.size+pos1);
		slv_int* start_pos1 = new slv_int[matA.size];
		slv_int* end_pos1 = new slv_int[matA.size];
		matA.getCols(start_pos1, end_pos1);
		auto col_ptr1 = matA.getColPtr();
		auto val_ptr1 = matA.getValuePtr();

		Mat1 tempA(max_size);
		/* まず自身をコピー */
		for(int i=0 ; i < matA.size ; i++){
			const slv_int the_size0 = end_pos1[i];
			for(slv_int aj = start_pos1[i] ; aj < the_size0 ; aj++){
				slv_int row = col_ptr1[aj];
				DType1 avalue = val_ptr1[aj];
				tempA.add(i, row, avalue);
			}
		}
		delete[] start_pos1;
		delete[] end_pos1;

		/* 足す行列をコピー */
		const slv_int size2 = matB.size;
		slv_int* start_pos2 = new slv_int[size2];
		slv_int* end_pos2 = new slv_int[size2];
		matB.getCols(start_pos2, end_pos2);
		auto col_ptr2 = matB.getColPtr();
		auto val_ptr2 = matB.getValuePtr();
		for(slv_int i = 0 ; i < size2 ; i++){
			const slv_int the_size = end_pos2[i];
			/* i行目のAが持つ列を探索 */
			for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
				slv_int row = col_ptr2[aj];
				DType2 avalue = val_ptr2[aj];
				avalue *= a1;
				tempA.add(i+pos1, row+pos2, avalue);
			}
		}
		delete[] start_pos2;
		delete[] end_pos2;
		/* 確定 */
		tempA.refresh();
		matA = std::move(tempA);
	}else{
		const slv_int size2 = matB.size;
		slv_int* start_pos2 = new slv_int[size2];
		slv_int* end_pos2 = new slv_int[size2];
		matB.getCols(start_pos2, end_pos2);
		auto col_ptr2 = matB.getColPtr();
		auto val_ptr2 = matB.getValuePtr();
		for(slv_int i = 0 ; i < size2 ; i++){
			const slv_int the_size = end_pos2[i];
			/* i行目のAが持つ列を探索 */
			for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
				slv_int row = col_ptr2[aj];
				DType2 avalue = val_ptr2[aj];
				avalue *= a1;
				matA.add(i+pos1, row+pos2, avalue);
			}
		}
		delete[] start_pos2;
		delete[] end_pos2;
	}
}


/*//=======================================================
// ● (mat1 + mat2)*vecBを計算
//=======================================================*/
template<typename Mat1, typename Mat2, typename DType1, typename DType2, typename DType0, typename DTypeA> 
void SparseMatOperators::productVecMat2(DTypeA** ans, const Mat1& mat1, const Mat2& mat2, const DType0* vecB){
	const slv_int size1 = mat1.size;
	const slv_int size2 = mat2.size;
	const slv_int min_size = (size1 > size2 ? size2 : size1);
	const slv_int max_size = (size1 > size2 ? size1 : size2);
	DTypeA* results = new DTypeA[max_size];

	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	mat1.getCols(start_pos1, end_pos1);
	auto col_ptr1 = mat1.getColPtr();
	auto val_ptr1 = mat1.getValuePtr();

	slv_int* start_pos2 = new slv_int[size2];
	slv_int* end_pos2 = new slv_int[size2];
	mat2.getCols(start_pos2, end_pos2);
	auto col_ptr2 = mat2.getColPtr();
	auto val_ptr2 = mat2.getValuePtr();

#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < min_size ; i++){
		const slv_int the_size = end_pos1[i];
		DTypeA temp=0;		
		for(slv_int j = start_pos1[i] ; j < the_size ; j++){
			temp += val_ptr1[j] * vecB[col_ptr1[j]];
		}
		const slv_int the_size2 = end_pos2[i];
		for(slv_int j = 0 ; j < the_size2 ; j++){
			temp += val_ptr2[j] * vecB[col_ptr2[j]];
		}
		results[i] = temp;
	}
	/* 行列１の方が行数が多いとき */
	if(size1 > size2){
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
		for(slv_int i = min_size ; i < size1 ; i++){
			const slv_int the_size = end_pos1[i];
			DTypeA temp=0;
			for(slv_int j = start_pos1[i] ; j < the_size ; j++){
				temp += val_ptr1[j] * vecB[col_ptr1[j]];
			}
			results[i] = temp;
		}
		/* 行列２の方が行数が多いとき */
	}else if(size1 < size2){
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
		for(slv_int i = min_size ; i < size2 ; i++){
			const slv_int the_size2 = end_pos2[i];
			DTypeA temp=0;
			for(slv_int j = 0 ; j < the_size2 ; j++){
				temp += val_ptr2[j] * vecB[col_ptr2[j]];
			}
			results[i] = temp;
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	delete[] start_pos2;
	delete[] end_pos2;
	*ans = std::move(results);
	results = nullptr;
}


/*//=======================================================
// ● 行列AとBとCをかけて結果を返す
//=======================================================*/
template<typename Mat1, typename Mat2, typename Mat3, typename Mat0, typename DType1, typename DType2, typename DType3, typename DType0> 
void SparseMatOperators::MatProduct(Mat0& mat_ans, const Mat1& matA, const Mat2& matB, const Mat3& matC){
	auto tempSparse = matA.matrix * matB.matrix * matC.matrix;
	//tempSparse.pruned();
	Mat0 matABC_damy( std::move(tempSparse) );
	mat_ans = std::move(matABC_damy);
}





/*//=======================================================*/
/*//=======================================================*/
/*//=======================================================*/
/*//=======================================================*/
/*//=======================================================*/


/*//=======================================================
// ● 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、返す
//
//    (r1a,r2a)  ..... (r1a,r2b)
//     .          .
//     .          .
//    (r1b,r2a) ...... (r1b,r2b)
//=======================================================*/
template<typename Mat1, typename DType>
void SparseMatOperators::makeSubMat(Mat1& mat_ans, const Mat1& targetMat, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b){
	const slv_int the_size = targetMat.size;
	const slv_int new_size = range1b-range1a + 1;
	Mat1 tempB(new_size);

	slv_int* start_pos1 = new slv_int[the_size];
	slv_int* end_pos1 = new slv_int[the_size];
	targetMat.getCols(start_pos1, end_pos1);
	auto col_ptr1 = targetMat.getColPtr();
	auto val_ptr1 = targetMat.getValuePtr();
	for(slv_int i = 0 ; i < the_size ; i++){
		if(range1a <= i && i <= range1b){
			const slv_int col_size = end_pos1[i];
			for(slv_int j = start_pos1[i] ; j < col_size ; j++){
				slv_int col = col_ptr1[j];
				if(range2a <= col && col <= range2b){
					DType val = val_ptr1[j];
					tempB.add(i-range1a, col-range2a, val);
				}
			}
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	/* 確定させて引き渡す */
	tempB.fix();
	mat_ans = std::move(tempB);
}

/*//=======================================================
// ● 列rangeBを境に、行列を2つに分ける。行rangeAより下は分けた行列には含めない
//       K11←B→K12
//    (A11) ..|... (A1n)
//      .     |      .
//rangeA______|_______
// ↓無視     |      .
//    (Am1) ..|... (Amn)
//=======================================================*/
template<typename Mat1, typename DType>
void SparseMatOperators::MatDiv(const Mat1& targetMat, Mat1& matK11, Mat1& matK12, slv_int rangeA, slv_int rangeB){
	const slv_int the_size = targetMat.size;
	const slv_int gyo_limit = rangeA+1;
	const slv_int min_gyo = (gyo_limit >= the_size ? the_size : gyo_limit);
	slv_int* start_pos1 = new slv_int[the_size];
	slv_int* end_pos1 = new slv_int[the_size];
	targetMat.getCols(start_pos1, end_pos1);
	auto col_ptr1 = targetMat.getColPtr();
	auto val_ptr1 = targetMat.getValuePtr();

	Mat1 tempK11(min_gyo);
	Mat1 tempK12(min_gyo);
	for(slv_int i = 0 ; i < min_gyo ; i++){
		const slv_int col_size = end_pos1[i];
		for(slv_int j = start_pos1[i] ; j < col_size ; j++){
			slv_int col = col_ptr1[j];
			DType val = val_ptr1[j];
			if(col <= rangeB){
				tempK11.add(i, col, val);
			}else{
				tempK12.add(i, col, val);
			}
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	tempK11.fix();
	tempK12.fix();
	matK11 = std::move(tempK11);
	matK12 = std::move(tempK12);
}

/*//=======================================================
// ● 下三角行列を取得
//=======================================================*/
template<typename Mat1, typename DType>
void SparseMatOperators::getMatLower(Mat1& ans_mat, const Mat1& targetMat){
	const slv_int the_size = targetMat.size;
	Mat1 tempMat(the_size);

	slv_int* start_pos1 = new slv_int[the_size];
	slv_int* end_pos1 = new slv_int[the_size];
	targetMat.getCols(start_pos1, end_pos1);
	auto col_ptr1 = targetMat.getColPtr();
	auto val_ptr1 = targetMat.getValuePtr();
	/* この行列の下三角部分だけ持ってくる */
	std::vector<slv_int> temp_col_size;
	for (slv_int i = 0; i < the_size; i++) {
		const slv_int the_temp_col_num = end_pos1[i];
		for(slv_int j = start_pos1[i] ; j < the_temp_col_num ; j++){
			/* 対角を発見 */
			if( col_ptr1[j] == i){
				temp_col_size.push_back(j+1);
				break;
			/* 対角を超えてしまった */
			}else if(col_ptr1[j] > i){
				/* 最初から超えていたら～ゼロ */
				if (j == 0) {
					temp_col_size.push_back(-1);
				/* それ以外は1つ前に */
				}else {
					temp_col_size.push_back(j);

				}
				break;
			}
			/* 見つからなかったら最後まで */
			if(j == the_temp_col_num-1){
				temp_col_size.push_back(the_temp_col_num);
			}
		}
	}
	for (slv_int i = 0; i < the_size; i++) {
		const slv_int the_size2 = temp_col_size[i];//tempMat->column_size[i];
		for(slv_int j = start_pos1[i] ; j < the_size2 ; j++){
			tempMat.add(i, col_ptr1[j], val_ptr1[j]);
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	tempMat.fix();
	ans_mat = std::move(tempMat);
}

/*//=======================================================
// ● 上三角行列を取得
//=======================================================*/
template<typename Mat1, typename DType>
void SparseMatOperators::getMatUpper(Mat1& ans_mat, const Mat1& targetMat){
	const slv_int the_size = targetMat.size;
	Mat1 tempMat(the_size);

	slv_int* start_pos1 = new slv_int[the_size];
	slv_int* end_pos1 = new slv_int[the_size];
	targetMat.getCols(start_pos1, end_pos1);
	auto col_ptr1 = targetMat.getColPtr();
	auto val_ptr1 = targetMat.getValuePtr();

	/* この行列の上三角部分だけ持ってくる */
	std::vector<slv_int> temp_col_size;
	for (slv_int i = 0; i < the_size; i++) {
		const slv_int the_temp_col_num = end_pos1[i];
		if(the_temp_col_num <= 0){
			temp_col_size.push_back(-1);
			continue;
		}
		for(slv_int j = start_pos1[i] ; j < the_temp_col_num ; j++){
			/* 対角以上の位置を発見 */
			if( col_ptr1[j] >= i){
				temp_col_size.push_back(j);
				//start_pos = j;
				break;
			}
			/* 対角がなかったら-1 */
			if(j == the_temp_col_num - 1) {
				temp_col_size.push_back(-1);
			}
		}
	}
	for (slv_int i = 0; i < the_size; i++) {
		const slv_int the_temp_col_num = end_pos1[i];
		const slv_int the_size2 = temp_col_size[i];//tempMat->column_size[i];
		if(the_size2 > -1){
			for(slv_int j = the_size2; j < the_temp_col_num; j++) {
				tempMat.add(i, col_ptr1[j], val_ptr1[j]);
			}
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	tempMat.fix();
	ans_mat = std::move(tempMat);
}

/*//=======================================================
// ● 行列Aの逆行列を求めてInvに渡す(密用)
//=======================================================*/
template<typename Mat1, typename DType, typename MType>	
void SparseMatOperators::MatInv(Mat1& ans_mat, const Mat1& targetMat){
	const slv_int the_size = targetMat.size;
	/* 一度フル格納形式へ */
	MType tempMat0 = targetMat.matrix;
	MType tempMat2 = tempMat0.inverse();
	/* スパース形式に代入し、渡す */
	Mat1 tempMat(the_size);
	for(slv_int i = 0 ; i < the_size ; i++){
		for(slv_int j = 0 ; j < the_size ; j++){
			DType value = tempMat2(i, j);
			double abs_val = abs(value);
			if(abs_val > 1.0e-12){
				tempMat.add(i, j, value);
			}
		}
	}
	tempMat.fix();
	ans_mat = std::move(tempMat);
}


/*//=======================================================
// ● 不完全コレスキー分解
//=======================================================*/
template<typename Mat1, typename DType1, typename DTypeD>
void SparseMatOperators::IC_decomp(Mat1& mat_ans, const Mat1& mat1, DTypeD* diagD, const double accela){
	if(!mat1.is_fix || mat1.tempMat != nullptr ){
		std::cout << "Error ! target mat in ICdecomp must be clear!"<< std::endl;
		getchar();
		exit(1);
	}
	const slv_int the_size = mat1.size;
	for (slv_int i = 0; i < the_size; i++) {
		diagD[i] = 0;
	}

	/* この行列の下三角部分だけ持ってくる */
	Mat1 tempMat;
	SparseMatOperators::getMatLower<Mat1, DType1>(tempMat, mat1);

	slv_int* start_pos1 = new slv_int[the_size];
	slv_int* end_pos1 = new slv_int[the_size];
	tempMat.getCols(start_pos1, end_pos1);
	auto col_ptr1 = tempMat.getColPtr();
	auto val_ptr1 = tempMat.getValuePtr();

	/* LとDを求める */
	DType1 s;
	/* LとDを求める */
	for (slv_int i = 0; i < the_size; i++) {
		/* i行目にある非ゼロの数をゲット */
		const slv_int c_size = end_pos1[i];
		//DType1 *A_ptr = mat1.matrix[i];
		//slv_int *C_ptr = mat1.column[i];
		/* i行目の非ゼロの回数だけ列ループをまわす */
		for (slv_int jj = start_pos1[i]; jj < c_size; jj++) {
			/* 現在の列番号jを取得 */
			const slv_int j = col_ptr1[jj];
			if(j >= i) break;
			s = val_ptr1[jj];
			/* i行目の列の非ゼロ列をサーチ */
			const slv_int L_ksize = end_pos1[i];
			//DType1 *L_ptr = tempMat->matrix[i];
			//slv_int *K_ptr = tempMat->column[i];
			for(slv_int kk = start_pos1[i] ; kk < L_ksize ; kk++){
				if(col_ptr1[kk] <= j-1){
					const slv_int L_Jsize = end_pos1[j];;//tempMat->column_size[j];
					//DType1 *L_Jptr = tempMat->matrix[j];
					//slv_int *KJ_ptr = tempMat->column[j];
					for(slv_int ll = start_pos1[j] ; ll < L_Jsize ; ll++){
						if(col_ptr1[ll] == col_ptr1[kk]){
							s -= val_ptr1[kk] * val_ptr1[ll] * diagD[col_ptr1[ll]];
							break;
						}else if(col_ptr1[kk] < col_ptr1[ll]){
							break;
						}
					}
				}else{
					break;
				}
			}
			/* 値を[i][j]に代入 */
			tempMat.matrix.coeffRef(i, j) = s;
			//tempMat->column[i][jj] = j;
			//tempMat->matrix[i][jj] = s;
		}
		const slv_int last = end_pos1[i] - 1;
		//std::cout << "diag " << i << ", " << last << ", "<< col_ptr1[last] << std::endl;
		/* もし対角が無ければエラーで落とす */
		if(col_ptr1[last] != i){
			std::cout << "there is not diagonal element! "<< std::endl;
			exit(1);
		}
		s = val_ptr1[last] * accela;
		//slv_int * IC_ptr = tempMat->column[i];
		//DType1 *IL_ptr = tempMat->matrix[i];
		for (slv_int kk = start_pos1[i]; kk < last; kk++) {
			//double temp_x;
			if(col_ptr1[kk] <= i-1){
				s -= val_ptr1[kk] * val_ptr1[kk] * diagD[col_ptr1[kk]];
			}else{
				break;
			}
		}
		/* 対角に代入 */
		tempMat.matrix.coeffRef(i, col_ptr1[last]) = s;
		//tempMat->matrix[i][last] = s;
		try{
			diagD[i] = 1.0 / s;
		}catch(...){
			diagD[i] = 1.0 / (DBL_EPSILON);
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	/* コピーして終わる */
	mat_ans = std::move(tempMat);
}




/*//=======================================================
// ● A^T*A + epsIを作る（疑似逆行列用）
//=======================================================*/
template<typename Mat1, typename DType1>
void SparseMatOperators::AtA_eps(Mat1& mat_ans, const Mat1& mat1, double eps){
	const slv_int size1 = mat1.getSize();
	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	mat1.getCols(start_pos1, end_pos1);
	auto col_ptr1 = mat1.getColPtr();//.matrix.innerIndexPtr();
	auto val_ptr1 = mat1.getValuePtr();//matrix.valuePtr();

	const slv_int new_size = 1 + mat1.getMaxCol();
	/* まずA^TとAをかける */
	Mat1 mat_AtA(new_size);	

	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr1[aj];
			DType1 avalue = val_ptr1[aj];
			/* Aが列をもつ場所rowのに対応したCの行を探索 */
			const slv_int the_sizeb = end_pos1[i];
			for(slv_int bj = start_pos1[i] ; bj < the_sizeb ; bj++){
				slv_int rowB = col_ptr1[bj];
				DType1 bvalue = val_ptr1[bj];
				mat_AtA.add(row, rowB, avalue*bvalue);
			}
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;

	/* 最後に対角項を足す */
	for(slv_int i = 0 ; i <new_size ; i++){
		mat_AtA.add(i, i, eps);
	}
	mat_AtA.fix();
	mat_ans = std::move(mat_AtA);
}


/*============================================*/
/*============================================*/
/*============================================*/
/**/
/* nonzero既知固定用 */
/**/

/*//=======================================================
// ● (配列位置確定時)自身に、行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす
//=======================================================*/
template<typename Mat1, typename Mat2, typename Mat0, typename DType1, typename DType2, typename DType0> 
void SparseMatOperators::PlusMinusShiftFix(Mat0& matAB, const Mat1& matA, const Mat2& matB, double a1, double a2, slv_int pos1, slv_int pos2){
	if(!matAB.is_fix){
		return;
	}
	matAB.resetMat();

	const slv_int size1 = matA.size;
	const slv_int size2 = matB.size;
	//const slv_int max = (size1 > size2+pos1 ? size1 : size2+pos1);

	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr1[aj];
			DType0 avalue = val_ptr1[aj];
			avalue *= a1;
			matAB.matrix.coeffRef(i, row) = avalue;
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;

	slv_int* start_pos2 = new slv_int[size2];
	slv_int* end_pos2 = new slv_int[size2];
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr2[aj];
			DType0 avalue = val_ptr2[aj];
			avalue *= a2;
			matAB.matrix.coeffRef(i+pos1, row+pos2) += avalue;
		}
	}
	delete[] start_pos2;
	delete[] end_pos2;
}
/*//=======================================================
// ● (配列位置確定時)自身に、行列A±B±Cを加える(a1*A+a2*B+a3*C)
//=======================================================*/
template<typename Mat1, typename Mat2, typename Mat3, typename Mat0, typename DType1, typename DType2, typename DType3, typename DType0> 
void SparseMatOperators::PlusMinusShiftFix(Mat0& matABC, const Mat1& matA, const Mat2& matB, const Mat3& matC, double a1, double a2, double a3){
	if(!matABC.is_fix){
		return;
	}
	matABC.resetMat();

	const slv_int size1 = matA.size;
	const slv_int size2 = matB.size;
	const slv_int size3 = matC.size;
	slv_int max = size1;
	if(max < size2){
		max = size2;
}
	if(max < size3){
		max = size3;
	}

	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();//matrix.innerIndexPtr();
	auto val_ptr1 = matA.getValuePtr();//matrix.valuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr1[aj];
			DType0 avalue = val_ptr1[aj];
			avalue *= a1;
			matABC.matrix.coeffRef(i, row) = avalue;
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;

	slv_int* start_pos2 = new slv_int[size2];
	slv_int* end_pos2 = new slv_int[size2];
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.matrix.innerIndexPtr();
	auto val_ptr2 = matB.matrix.valuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr2[aj];
			DType0 avalue = val_ptr2[aj];
			avalue *= a2;
			matABC.matrix.coeffRef(i, row) += avalue;
		}
	}
	delete[] start_pos2;
	delete[] end_pos2;

	slv_int* start_pos3 = new slv_int[size3];
	slv_int* end_pos3 = new slv_int[size3];
	matC.getCols(start_pos3, end_pos3);
	auto col_ptr3 = matC.matrix.innerIndexPtr();
	auto val_ptr3 = matC.matrix.valuePtr();
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size3 ; i++){
		const slv_int the_size = end_pos3[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos3[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr3[aj];
			DType0 avalue = val_ptr3[aj];
			avalue *= a3;
			matABC.matrix.coeffRef(i, row) += avalue;
		}
	}
	delete[] start_pos3;
	delete[] end_pos3;
}



/*//=======================================================
// ● (配列位置確定時)行列AとBをかけて、自身ABに格納する
//=======================================================*/
template<typename Mat1, typename Mat2, typename Mat0, typename DType1, typename DType2, typename DType0> 
void SparseMatOperators::MatProductFix(Mat0& matAB, const Mat1& matA, const Mat2& matB){
	if(!matAB.is_fix){
		return;
	}
	matAB.resetMat();

	const slv_int size1 = matA.getSize();
	slv_int* start_pos1 = new slv_int[size1];
	slv_int* end_pos1 = new slv_int[size1];
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();//matrix.innerIndexPtr();
	auto val_ptr1 = matA.getValuePtr();//matrix.valuePtr();

	const slv_int size2 = matB.getSize();
	slv_int* start_pos2 = new slv_int[size2];
	slv_int* end_pos2 = new slv_int[size2];
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();//matrix.innerIndexPtr();
	auto val_ptr2 = matB.getValuePtr();//matrix.valuePtr();

#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr1[aj];
			DType0 avalue = val_ptr1[aj];
			/* Aが列をもつ場所rowのに対応したCの行を探索 */
			const slv_int the_sizeb = end_pos2[row];
			for(slv_int bj = start_pos2[row] ; bj < the_sizeb ; bj++){
				slv_int rowB = col_ptr2[bj];
				DType0 bvalue = val_ptr2[bj];
				matAB.matrix.coeffRef(i, rowB) += avalue*bvalue;
				//matAB.add_typeN(i, rowB, avalue*bvalue);
			}
		}
	}
	delete[] start_pos1;
	delete[] end_pos1;
	delete[] start_pos2;
	delete[] end_pos2;
}



/* end of namespace */
};


#endif


