/**
 * @file SprsOperator_MatOps.cpp
 * @brief  Four arithmetic matrix operations for SparseMatTMPL
 * 
 * このファイルは SparseMatTMPL の行列四則演算の実装本体です
 */

#include <SparseSolv/SprsOperator.hpp>
#include <SparseSolv/BasicDefinesEigenCore.hpp>


/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/* 疎行列演算用の専用名前空間 */
namespace SprsOperator {



/* 疎行列演算用内部空間間 */
/**
* @namespace SprsOperator::detai
* @brief Original namespace for internal use of SprsOperator
*/
namespace detail {


/*//=======================================================
// ● 行列ベクトル積・生ポインタ(内部用)
//=======================================================*/
/** 
 * @brief mat * pointer* operator (static, private), main body
*/
template<typename DTypeMat, typename DtypeV, typename DTypeResult> 
void mat_vec_product(DTypeResult*& ans, const DtypeV* vec, 
	const slv_int size, const std::vector<slv_int>& start_pos, const std::vector<slv_int>& end_pos, const slv_int* col_ptr, const DTypeMat* val_ptr, bool init){	

	/* 掛け算開始 */
	if(init){
		ans = new DTypeResult[size];
	}
	for(slv_int i = 0 ; i < size ; i++){
		DTypeResult temp = 0;
		const slv_int ss = end_pos[i];
		for(slv_int j = start_pos[i] ; j < ss ; j++){
			const slv_int pos = col_ptr[j];
			temp += vec[pos]*val_ptr[j];
		}
		ans[i] = temp;
	}
}


/*//=======================================================
// ● (mat1 + mat2)*vecBを計算・生ポインタ(内部用)
//=======================================================*/
/** 
 * @brief @brief vec = (matA+matB)*vecB, main body for pointer (static, private)
*/
template<typename DTypeMat, typename DtypeV, typename DTypeResult> 
void mat_vec_product2(DTypeResult*& ans, const DtypeV* const vec, 
	const slv_int size1, const std::vector<slv_int>& start_pos1, const std::vector<slv_int>& end_pos1, const slv_int* col_ptr1, const DTypeMat* val_ptr1,
	const slv_int size2, const std::vector<slv_int>& start_pos2, const std::vector<slv_int>& end_pos2, const slv_int* col_ptr2, const DTypeMat* val_ptr2,
	bool init){

	
	const slv_int min_size = (size1 > size2 ? size2 : size1);
	const slv_int max_size = (size1 > size2 ? size1 : size2);

	if(init){
		ans = new DTypeResult[max_size];
	}

#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
	for(slv_int i = 0 ; i < min_size ; i++){
		const slv_int the_size = end_pos1[i];
		DTypeResult temp=0;		
		for(slv_int j = start_pos1[i] ; j < the_size ; j++){
			temp += val_ptr1[j] * vec[col_ptr1[j]];
		}
		const slv_int the_size2 = end_pos2[i];
		for(slv_int j = start_pos2[i] ; j < the_size2 ; j++){
			temp += val_ptr2[j] * vec[col_ptr2[j]];
		}
		ans[i] = temp;
	}
	/* 行列１の方が行数が多いとき */
	if(size1 > size2){
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
		for(slv_int i = min_size ; i < size1 ; i++){
			const slv_int the_size = end_pos1[i];
			DTypeResult temp=0;
			for(slv_int j = start_pos1[i] ; j < the_size ; j++){
				temp += val_ptr1[j] * vec[col_ptr1[j]];
			}
			ans[i] = temp;
		}
		/* 行列２の方が行数が多いとき */
	}else if(size1 < size2){
#ifdef OMP_USING_MAT_SOL
#pragma omp parallel for
#endif
		for(slv_int i = min_size ; i < size2 ; i++){
			const slv_int the_size2 = end_pos2[i];
			DTypeResult temp=0;
			for(slv_int j = start_pos2[i] ; j < the_size2 ; j++){
				temp += val_ptr2[j] * vec[col_ptr2[j]];
			}
			ans[i] = temp;
		}
	}
}



/*//=======================================================
// ● 行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす
//=======================================================*/
/**
 * @brief a1*matA + a2*matB
*/
template<typename DMat>
void MatPlusMinusShift_sub(Eigen::SparseMatrix<DMat, Eigen::RowMajor>& matAB, double a1, double a2, slv_int pos1, slv_int pos2,
	const slv_int size1, const std::vector<slv_int>& start_pos1, const std::vector<slv_int>& end_pos1, const slv_int* col_ptr1, const DMat* val_ptr1,
	const slv_int size2, const std::vector<slv_int>& start_pos2, const std::vector<slv_int>& end_pos2, const slv_int* col_ptr2, const DMat* val_ptr2) {

	slv_int max1 = 0;
	slv_int max2 = 0;
	std::vector<Eigen::Triplet<DMat>> sparse_data;

	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr1[aj];
			DMat avalue = val_ptr1[aj];
			avalue *= a1;
			sparse_data.push_back(Eigen::Triplet<DMat>(i, row, avalue));
			if (max1 < i) {
				max1 = i;
			}
			if(max2 < row){
				max2 = row;
			}
		}
	}
	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr2[aj];
			DMat avalue = val_ptr2[aj];
			avalue *= a2;
			sparse_data.push_back(Eigen::Triplet<DMat>(i+pos1, row+pos2, avalue));
			if (max1 < i+pos1) {
				max1 = i+pos1;
			}
			if(max2 < row+pos2){
				max2 = row+pos2;
			}
		}
	}
	max1++;
	max2++;
	matAB = Eigen::SparseMatrix<DMat, Eigen::RowMajor>(max1, max2);
	matAB.setFromTriplets(sparse_data.begin(), sparse_data.end());
}


/*//=======================================================
// ● matAにa1*matBを加える。 開始位置はposだけずらす サブルーチン
//=======================================================*/
/** 
 * @brief add MatA += a1*MatB (matB row and col are shifted) subroutine
 * @param matA
 * @param matB
 * @param a1
 * @param pos1 shift param for row
 * @param pos2 shift param for col
*/
template<typename DType1> 
void MatSelfPlus_sub(SparseMatTMPL<DType1>& matA, const Eigen::SparseMatrix<DType1, Eigen::RowMajor>& matB, const DType1 a1, const slv_int pos1, const slv_int pos2,
	const slv_int size1, const std::vector<slv_int>& start_pos1, const std::vector<slv_int>& end_pos1, const slv_int* col_ptr1, const DType1* val_ptr1,
	const slv_int size2, const std::vector<slv_int>& start_pos2, const std::vector<slv_int>& end_pos2, const slv_int* col_ptr2, const DType1* val_ptr2) {

	const slv_int max_size = (size1 > size2+pos1 ? size1 : size2+pos1);

	SparseBuilderTMPL<DType1> builderAB(max_size);
	/* まず自身をコピー */
	for(int i=0 ; i < size1 ; i++){
		const slv_int the_size0 = end_pos1[i];
		for(slv_int aj = start_pos1[i] ; aj < the_size0 ; aj++){
			slv_int row = col_ptr1[aj];
			DType1 avalue = val_ptr1[aj];
			builderAB.add(i, row, avalue);
		}
	}

	/* 足す行列をコピー */
	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int row = col_ptr2[aj];
			DType1 avalue = val_ptr2[aj];
			avalue *= a1;
			builderAB.add(i+pos1, row+pos2, avalue);
		}
	}
	/* 確定 */
	matA.build(builderAB, false);
}



/*end of sub namespace*/
};


/*==================================================================*/
/*==================================================================*/
/* 実装本体*/
/*==================================================================*/
/*==================================================================*/
/*==================================================================*/



/*//=======================================================
// ● 行列の足し算
//=======================================================*/
/** 
 * @brief matrix addition
*/
SparseMat add(const SparseMat& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = matA.mat() + matB.mat();
	return SparseMat( std::move(out) );
}
SparseMatC add(const SparseMatC& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() + matB.mat().cast<dcomplex>();
	return SparseMatC( std::move(out) );
}
SparseMatC add(const SparseMat& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat().cast<dcomplex>() + matB.mat();
	return SparseMatC( std::move(out) );
}
SparseMatC add(const SparseMatC& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() + matB.mat();
	return SparseMatC( std::move(out) );
}


/*//=======================================================
// ● 行列の引き算
//=======================================================*/
/** 
 * @brief matrix subtraction
*/
SparseMat substruct(const SparseMat& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = matA.mat() - matB.mat();
	return SparseMat( std::move(out) );
}
SparseMatC substruct(const SparseMatC& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() - matB.mat().cast<dcomplex>();
	return SparseMatC( std::move(out) );
}
SparseMatC substruct(const SparseMat& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat().cast<dcomplex>() - matB.mat();
	return SparseMatC( std::move(out) );
}
SparseMatC substruct(const SparseMatC& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() - matB.mat();
	return SparseMatC( std::move(out) );
}



/*//=======================================================
// ● 行列の掛け算
//=======================================================*/
/** 
 * @brief matrix product
*/
SparseMat dot(const SparseMat& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = matA.mat() * matB.mat();
	return SparseMat( std::move(out) );
}
SparseMatC dot(const SparseMatC& matA, const SparseMat& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() * matB.mat().cast<dcomplex>();
	return SparseMatC( std::move(out) );
}
SparseMatC dot(const SparseMat& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat().cast<dcomplex>() * matB.mat();
	return SparseMatC( std::move(out) );
}
SparseMatC dot(const SparseMatC& matA, const SparseMatC& matB) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat() * matB.mat();
	return SparseMatC( std::move(out) );
}


/*//=======================================================
// ● 行列ベクトル演算（生ポインタver）
//=======================================================*/
/** 
 * @brief matrix-vector product for pointer (raw pointer version)
*/
double* dot(const SparseMat& matA, const double* vec) {
	double* ptr;
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const double* val_ptr = matA.getValuePtr();
	detail::mat_vec_product<double, double, double>(ptr, vec, size, start_pos, end_pos, col_ptr, val_ptr, true);
	return ptr;
}
dcomplex* dot(const SparseMat& matA, const dcomplex* vec) {
	dcomplex* ptr;
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const double* val_ptr = matA.getValuePtr();
	detail::mat_vec_product<double, dcomplex, dcomplex>(ptr, vec, size, start_pos, end_pos, col_ptr, val_ptr, true);
	return ptr;
}
dcomplex* dot(const SparseMatC& matA, const double* vec){
	dcomplex* ptr;
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const dcomplex* val_ptr = matA.getValuePtr();
	detail::mat_vec_product<dcomplex, double, dcomplex>(ptr, vec, size, start_pos, end_pos, col_ptr, val_ptr, true);
	return ptr;
}
dcomplex* dot(const SparseMatC& matA, const dcomplex* vec) {
	dcomplex* ptr;
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const dcomplex* val_ptr = matA.getValuePtr();
	detail::mat_vec_product<dcomplex, dcomplex, dcomplex>(ptr, vec, size, start_pos, end_pos, col_ptr, val_ptr, true);
	return ptr;
}


/*//=======================================================
// ● 行列ベクトル演算（std::vector ver）
//=======================================================*/
/** 
 * @brief matrix-vector product for vector (std::vector version)
*/
std::vector<double> dot(const SparseMat& matA, const std::vector<double>& vec) {
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const double* val_ptr = matA.getValuePtr();

	std::vector<double> ans;
	ans.resize(size);
	double* ptr = ans.data();
	detail::mat_vec_product<double, double, double>(ptr, vec.data(), size, start_pos, end_pos, col_ptr, val_ptr, false);
	return ans;
}
std::vector<dcomplex> dot(const SparseMat& matA, const std::vector<dcomplex>& vec) {	
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const double* val_ptr = matA.getValuePtr();

	std::vector<dcomplex> ans;
	ans.resize(size);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product<double, dcomplex, dcomplex>(ptr, vec.data(), size, start_pos, end_pos, col_ptr, val_ptr, false);
	return ans;
}
std::vector<dcomplex> dot(const SparseMatC& matA, const std::vector<double>& vec){
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const dcomplex* val_ptr = matA.getValuePtr();

	std::vector<dcomplex> ans;
	ans.resize(size);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product<dcomplex, double, dcomplex>(ptr, vec.data(), size, start_pos, end_pos, col_ptr, val_ptr, false);
	return ans;
}
std::vector<dcomplex> dot(const SparseMatC& matA, const std::vector<dcomplex>& vec) {
	slv_int size = matA.getSize();
	std::vector<slv_int> start_pos;
	std::vector<slv_int> end_pos;
	matA.getCols(start_pos, end_pos);
	const slv_int* col_ptr = matA.getColPtr();
	const dcomplex* val_ptr = matA.getValuePtr();

	std::vector<dcomplex> ans;
	ans.resize(size);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product<dcomplex, dcomplex, dcomplex>(ptr, vec.data(), size, start_pos, end_pos, col_ptr, val_ptr, false);
	return ans;
}


/*//=======================================================
// ● 行列ベクトル演算(Eigen) 
//=======================================================*/
/** 
 * @brief matrix-vector product based on Eigen defaults
*/
Eigen::VectorXd dot(const SparseMat& matA, const Eigen::VectorXd& vec){return matA.mat()*vec;};
Eigen::VectorXcd dot(const SparseMat& matA, const Eigen::VectorXcd& vec){return matA.mat().cast<dcomplex>()*vec;};
Eigen::VectorXcd dot(const SparseMatC& matA, const Eigen::VectorXd& vec){return matA.mat()*vec.cast<dcomplex>();};
Eigen::VectorXcd dot(const SparseMatC& matA, const Eigen::VectorXcd& vec){return matA.mat()*vec;};





/*--------------------------------------------------*/



/*//=======================================================
// ● 足し算オペレータ(matAにmatBを加える(a1*A+a2*B)。Bの位置を行方向にpos1、列方向にpos2だけずらす。)
//=======================================================*/
/** 
 * @brief add a1*MatA + a2*MatB (matB row and col are shifted)
 * @param matA
 * @param matB
 * @param a1
 * @param a2
 * @param pos1 shift param for row
 * @param pos2 shift param for col
*/
SparseMat ShiftPus(const SparseMat& matA, const SparseMat& matB, double a1, double a2, slv_int pos1, slv_int pos2) {
	const slv_int size1 = matA.getSize();
	const slv_int size2 = matB.getSize();
	const slv_int max = (size1 > size2+pos1 ? size1 : size2+pos1);

	/* 足した後の行列 */
	SparseBuilder builderAB(max);

	/* まずAを足す */
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr1[aj];
			double avalue = val_ptr1[aj];
			avalue *= a1;
			builderAB.add(i, col, avalue);
		}
	}

	/* Bをずらしながら足す */
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();

	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr2[aj];
			double avalue = val_ptr2[aj];
			avalue *= a2;
			builderAB.add(i+pos1, col+pos2, avalue);
		}
	}

	/* 確定させ、ムーブする */
	SparseMat matAB(builderAB, false);
	return matAB;
}

SparseMatC ShiftPus(const SparseMatC& matA, const SparseMatC& matB, double a1, double a2, slv_int pos1, slv_int pos2) {
	const slv_int size1 = matA.getSize();
	const slv_int size2 = matB.getSize();
	const slv_int max = (size1 > size2+pos1 ? size1 : size2+pos1);

	/* 足した後の行列 */
	SparseBuilderC builderAB(max);

	/* まずAを足す */
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	for(slv_int i = 0 ; i < size1 ; i++){
		const slv_int the_size = end_pos1[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos1[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr1[aj];
			dcomplex avalue = val_ptr1[aj];
			avalue *= a1;
			builderAB.add(i, col, avalue);
		}
	}

	/* Bをずらしながら足す */
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();

	for(slv_int i = 0 ; i < size2 ; i++){
		const slv_int the_size = end_pos2[i];
		/* i行目のAが持つ列を探索 */
		for(slv_int aj = start_pos2[i] ; aj < the_size ; aj++){
			slv_int col = col_ptr2[aj];
			dcomplex avalue = val_ptr2[aj];
			avalue *= a2;
			builderAB.add(i+pos1, col+pos2, avalue);
		}
	}

	/* 確定させ、ムーブする */
	SparseMatC matAB(builderAB, false);
	return matAB;
}




/*//=======================================================
// ● matAにa1*matBを加える。 開始位置はposだけずらす
//=======================================================*/
/** 
 * @brief add MatA += a1*MatB (matB row and col are shifted)
 * @param matA
 * @param matB
 * @param a1
 * @param pos1 shift param for row
 * @param pos2 shift param for col
*/
void MatSelfPlus(SparseMat& matA, const SparseMat& matB, const double a1, const slv_int pos1, const slv_int pos2){
	const slv_int size1 = matA.getSize();
	const slv_int size2 = matB.getSize();

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();

	detail::MatSelfPlus_sub<double>(matA, matB.mat(), a1, pos1, pos2,
		size1, start_pos1, end_pos1, col_ptr1, val_ptr1,
		size2, start_pos2, end_pos2, col_ptr2, val_ptr2);
}
void MatSelfPlus(SparseMatC& matA, const SparseMatC& matB, const dcomplex a1, const slv_int pos1, const slv_int pos2){
	const slv_int size1 = matA.getSize();
	const slv_int size2 = matB.getSize();

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	auto col_ptr2 = matB.getColPtr();
	auto val_ptr2 = matB.getValuePtr();

	detail::MatSelfPlus_sub<dcomplex>(matA, matB.mat(), a1, pos1, pos2,
		size1, start_pos1, end_pos1, col_ptr1, val_ptr1,
		size2, start_pos2, end_pos2, col_ptr2, val_ptr2);
}



/*//=======================================================
// ● (mat1 + mat2)*vecBを計算
//=======================================================*/
/** 
 * @brief vec = (matA+matB)*vecB
 * @param mat1
 * @param mat2
 * @param vecB
*/double* dot2(const SparseMat& matA, const SparseMat& matB, const double* vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();

	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	double* ptr;
	detail::mat_vec_product2<double, double, double>(ptr, vec, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, true);
	return ptr;
}
dcomplex* dot2(const SparseMatC& matA, const SparseMatC& MatB, const double* vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = MatB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	MatB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = MatB.getColPtr();
	const dcomplex* val_ptr2 = MatB.getValuePtr();
	
	dcomplex* ptr;
	detail::mat_vec_product2<dcomplex, double, dcomplex>(ptr, vec, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, true);
	return ptr;
}
dcomplex* dot2(const SparseMat& matA, const SparseMat& MatB, const dcomplex* vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();
	slv_int size2 = MatB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	MatB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = MatB.getColPtr();
	const double* val_ptr2 = MatB.getValuePtr();

	dcomplex* ptr;
	detail::mat_vec_product2<double, dcomplex, dcomplex>(ptr, vec, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, true);
	return ptr;
}
dcomplex* dot2(const SparseMatC& matA, const SparseMatC& MatB, const dcomplex* vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = MatB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	MatB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = MatB.getColPtr();
	const dcomplex* val_ptr2 = MatB.getValuePtr();

	dcomplex* ptr;
	detail::mat_vec_product2<dcomplex, dcomplex, dcomplex>(ptr, vec, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, true);
	return ptr;
}

/*//=======================================================
// ● (mat1 + mat2)*vecBを計算
//=======================================================*/
/** 
 * @brief vec = (matA+matB)*vecB for vector version
 * @param mat1
 * @param mat2
 * @param vecB
 */
std::vector<double> dot2(const SparseMat& matA, const SparseMat& matB, const std::vector<double>& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();

	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	std::vector<double> ans;
	ans.resize(max);
	double* ptr = ans.data();
	detail::mat_vec_product2<double, double, double>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}

std::vector<dcomplex> dot2(const SparseMatC& matA, const SparseMatC& matB, const std::vector<double>& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const dcomplex* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	std::vector<dcomplex> ans;
	ans.resize(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<dcomplex, double, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}

std::vector<dcomplex> dot2(const SparseMat& matA, const SparseMat& matB, const std::vector<dcomplex>& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	std::vector<dcomplex> ans;
	ans.resize(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<double, dcomplex, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}

std::vector<dcomplex> dot2(const SparseMatC& matA, const SparseMatC& matB, const std::vector<dcomplex>& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const dcomplex* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	std::vector<dcomplex> ans;
	ans.resize(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<dcomplex, dcomplex, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}


/*//=======================================================
// ● (mat1 + mat2)*vecBを計算 for Eigen
//=======================================================*/
/** 
 * @brief vec = (matA+matB)*vecB for Eigenr version
 * @param mat1
 * @param mat2
 * @param vecB
 */
Eigen::VectorXd dot2(const SparseMat& matA, const SparseMat& matB, const Eigen::VectorXd& vec) {
	
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();

	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	Eigen::VectorXd ans(max);
	double* ptr = ans.data();
	detail::mat_vec_product2<double, double, double>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}
Eigen::VectorXcd dot2(const SparseMatC& matA, const SparseMatC& matB, const Eigen::VectorXd& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const dcomplex* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	Eigen::VectorXcd ans(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<dcomplex, double, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}

Eigen::VectorXcd dot2(const SparseMat& matA, const SparseMat& matB, const Eigen::VectorXcd& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	Eigen::VectorXcd ans(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<double, dcomplex, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}
Eigen::VectorXcd dot2(const SparseMatC& matA, const SparseMatC& matB, const Eigen::VectorXcd& vec) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();
	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const dcomplex* val_ptr2 = matB.getValuePtr();

	const slv_int max = (size1 > size2 ? size1 : size2);
	Eigen::VectorXcd ans(max);
	dcomplex* ptr = ans.data();
	detail::mat_vec_product2<dcomplex, dcomplex, dcomplex>(ptr, vec.data(), size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2, false);
	return ans;
}


/*--------------------------------------------------*/


/*//=======================================================
// ● 行列AとBとCをかけて結果を返す
//=======================================================*/
/** 
* @brief matABC = (matA*matB*matC)
* @param matA
* @param matB
* @param matC
*/
SparseMat dot3(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = (matA.mat() * matB.mat()) * matC.mat();
	return SparseMat(std::move(out));
}
SparseMatC dot3(const SparseMat& matA, const SparseMatC& matB, const SparseMat& matC) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = (matA.mat().cast<dcomplex>() * matB.mat()) * matC.mat().cast<dcomplex>();
	return SparseMatC(std::move(out));
}
SparseMatC dot3(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = (matA.mat() * matB.mat()) * matC.mat();
	return SparseMatC(std::move(out));
}



/*//=======================================================
// ● 行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす 
//=======================================================*/
/** 
* @brief mataAC = (a1*A+a2*B)
*/
SparseMat MatPlusMinusShift(const SparseMat& matA, const SparseMat& matB, double a1, double a2, slv_int pos1, slv_int pos2) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const double* val_ptr1 = matA.getValuePtr();

	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const double* val_ptr2 = matB.getValuePtr();

	Eigen::SparseMatrix<double, Eigen::RowMajor> out;
	detail::MatPlusMinusShift_sub<double>(out, a1, a2, pos1, pos2, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2);
	return SparseMat(std::move(out));
}

SparseMatC MatPlusMinusShift(const SparseMatC& matA, const SparseMatC& matB, double a1, double a2, slv_int pos1, slv_int pos2) {
	slv_int size1 = matA.getSize();
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const slv_int* col_ptr1 = matA.getColPtr();
	const dcomplex* val_ptr1 = matA.getValuePtr();

	slv_int size2 = matB.getSize();
	std::vector<slv_int> start_pos2;
	std::vector<slv_int> end_pos2;
	matB.getCols(start_pos2, end_pos2);
	const slv_int* col_ptr2 = matB.getColPtr();
	const dcomplex* val_ptr2 = matB.getValuePtr();

	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out;
	detail::MatPlusMinusShift_sub<dcomplex>(out, a1, a2, pos1, pos2, size1, start_pos1, end_pos1, col_ptr1, val_ptr1, size2, start_pos2, end_pos2, col_ptr2, val_ptr2);
	return SparseMatC(std::move(out));
}



/*//=======================================================
// ● 行列A±B±Cを加える(a1*A+a2*B+a3*C)
//=======================================================*/
/** 
* @brief mataAC = a1*A+a2*B+a3*C
*/
SparseMat add3(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC, double a1, double a2, double a3) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = a1 * matA.mat() + a2 * matB.mat() + a3 * matC.mat();
	return SparseMat(std::move(out));
}
SparseMatC add3(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC, double a1, double a2, double a3) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = a1 * matA.mat() + a2 * matB.mat() + a3 * matC.mat();
	return SparseMatC(std::move(out));
}


/*  */

/*//=======================================================
// ● 行列Aの逆行列を求めてInvに渡す(密用)
//=======================================================*/
/** 
* @brief calc A^-1 using dense solver (do not use for large sparse matrices)
*/
SparseMat inv(const SparseMat& targetMat) {
	const slv_int the_size = targetMat.getSize();
	/* 一度フル格納形式へ */
	Eigen::MatrixXd tempMat0 = targetMat.mat();
	Eigen::MatrixXd tempMat2 = tempMat0.inverse();
	/* スパース形式に代入し、渡す */	
	SparseBuilder builder(the_size);
	for(slv_int i = 0 ; i < the_size ; i++){
		for(slv_int j = 0 ; j < the_size ; j++){
			double value = tempMat2(i, j);
			double abs_val = abs(value);
			if(abs_val > 1.0e-12){
				builder.add(i, j, value);
			}
		}
	}
	return( SparseMat(builder) );
}
SparseMatC inv(const SparseMatC& targetMat) {
	const slv_int the_size = targetMat.getSize();
	/* 一度フル格納形式へ */
	Eigen::MatrixXcd tempMat0 = targetMat.mat();
	Eigen::MatrixXcd tempMat2 = tempMat0.inverse();
	/* スパース形式に代入し、渡す */	
	SparseBuilderC builder(the_size);
	for(slv_int i = 0 ; i < the_size ; i++){
		for(slv_int j = 0 ; j < the_size ; j++){
			dcomplex value = tempMat2(i, j);
			double abs_val = std::abs(value);
			if(abs_val > 1.0e-12){
				builder.add(i, j, value);
			}
		}
	}
	return( SparseMatC(builder) );
}


/* end of namespace */
};
};
};
