/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file SparseMatTMPL.cpp
 * @brief Sparse Matrix Template Class
 * このファイルはテンプレートクラス SparseMatTMPL の実体を実装してtemplate公開します。
 */

#include <SparseSolv/SparseMat.hpp>
#include <string>
#include <fstream>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{


/*
//=======================================================
// ■ スパース行列テンプレート
//=======================================================*/



/*//=======================================================
// ● 行数を返す(サイズ範囲チェックありで)
//=======================================================*/
/** 
 * @brief get number of rows (with size range check)
*/
template<typename DType>
slv_int SparseMatTMPL<DType>::getSize() const {
    const Eigen::Index n = matrix.rows();
    if (n > std::numeric_limits<slv_int>::max()) {
        throw std::overflow_error("SparseMat size exceeds slv_int range");
    }
    return static_cast<slv_int>(n);
}

/*//=======================================================
// ● Eigenタプルから初期化
//=======================================================*/
/** 
 * @brief Constructor from triplet
*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(const std::vector<Eigen::Triplet<DType>>& sparse_data){
	const size_t data_size = sparse_data.size();
	/* 最大列をサーチ */
	slv_int max1 = 0;
	slv_int max2 = 0;
	for(size_t i = 0 ; i < data_size ; i++){
		slv_int gyo = sparse_data[i].row();
		slv_int retu = sparse_data[i].col();
		if(max1 < gyo) max1 = gyo;
		if(max2 < retu) max2 = retu;
	}
	max1++;
	max2++;
	/* 確定させる */
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(max1, max2);
	matrix.setFromTriplets(sparse_data.begin(), sparse_data.end());
}


/*//=======================================================
// ● 作成済み疎行列データから初期化
//=======================================================*/
/** 
 * @brief Constructor from non-zero data
*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<DType>& vals){
	/* EigenのTripleにデータを代入しつつ、最大行/列を探す */
	slv_int max_gyo = 0;
	slv_int max_retu = 0;
	std::vector<Eigen::Triplet<DType>> tripletList;
	tripletList.reserve(n_nonZero);
	for(slv_int i = 0; i < n_nonZero; i++) {
		const slv_int gyo = rows[i];
		const slv_int retu = cols[i];
		//Eigen::Triplet<DType> tripl = Eigen::Triplet<DType>(rows[i], cols[i], vals[i]);
		//tripletList.push_back( std::move(tripl) );
		tripletList.emplace_back(gyo, retu, vals[i]);
		if(max_gyo < gyo){
			max_gyo = gyo;
		}
		if(max_retu < retu){
			max_retu = retu;
		}
	}
	max_gyo++;
	max_retu++;
	/* 確定処理 */
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(max_gyo, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}

/*//=======================================================
// ● Builderから作成: to_square=長方形の場合に正方形にするか
//=======================================================*/
/**
 * @brief Constructor from builder
 * @param to_square true: make square matrix by adding zeros
*/
template<typename DType>
void SparseMatTMPL<DType>::build(SparseBuilderTMPL<DType>& builder, bool to_square) {
	const slv_int size = builder.getSize();
	/* builderからtripletを取得 */
	std::vector<Eigen::Triplet<DType>> tripletList;
	slv_int max_retu = builder.build(tripletList, to_square);
	/* Eigenに渡して疎行列を確定 */
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(size, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}


/*//=======================================================
// ● スパース内の最大の列位置を返す
//=======================================================*/
/** 
 * @brief get Max col pos
*/
template<typename DType>
slv_int SparseMatTMPL<DType>::getMaxCol() const{
	const slv_int size = this->getSize();
	/* ポインタをたどって探す */
	slv_int count=0;
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	const Eigen::Index total_size = matrix.nonZeros();
	Eigen::Index max = 0;
	for(Eigen::Index i = 0; i < size; i++) {
		const Eigen::Index num = (i==size-1 ? total_size : col_ptr[i+1]);
		for(Eigen::Index j = col_ptr[i]; j < num; j++) {
			if(max < row_ptr[count]) {
				max = row_ptr[count];
			}
			count++;
		}
	}
	return static_cast<slv_int>(max);
}


/*//=======================================================
// ● 各行の列のスタート・終了位置を返す
//=======================================================*/
/** 
 * @brief find each start and end col's positions at non-zero elements
*/
template<typename DType>
void SparseMatTMPL<DType>::getCols(slv_int* cols_data1, slv_int* cols_data2) const{
	const Eigen::Index size_n = this->getTrueSize() - 1;
	auto col_ptr = matrix.outerIndexPtr();
	for(Eigen::Index i = 0; i < size_n; i++) {
		cols_data1[i] = col_ptr[i];
		cols_data2[i] = col_ptr[i+1];
	}
	cols_data1[size_n] = col_ptr[size_n];
	cols_data2[size_n] = static_cast<slv_int>(matrix.nonZeros());
}
/** 
 * @brief find each start and end col's positions at non-zero elements (vector)
*/
template<typename DType>
void SparseMatTMPL<DType>::getCols(std::vector<slv_int>& cols_data1, std::vector<slv_int>& cols_data2) const{
	const Eigen::Index size = this->getTrueSize();
	cols_data1.resize(size);
	cols_data2.resize(size);
	this->getCols(cols_data1.data(), cols_data2.data());
}

/*//=======================================================
// ● 値だけゼロにする（位置は保持したまま）
//=======================================================*/
/** 
 * @brief reset matrix value to zero (col and row positions are unchanged)
*/
template<typename DType>
void SparseMatTMPL<DType>::resetMat(){
	/* ポインタをたどってゼロにする */
	auto val_ptr = matrix.valuePtr();
	const auto total_size = matrix.nonZeros();
	for(auto i = 0; i < total_size; i++) {
		val_ptr[i] = 0.0;
	}
}

/*//=======================================================
// ● i行目にtarget_r列があるかどうか（あったらそのindexを返す）
//=======================================================*/
/** 
 * @brief Does target_row include?
 * @param gyo target row
 * @param target_r target col
 * @param result_retu& stored result target col number
 * @return true/false if included, return true
*/
template<typename DType>
bool SparseMatTMPL<DType>::isInclude(slv_int gyo, slv_int target_r, slv_int& result_retu) const{
	//const slv_int size = this->getSize();
	Eigen::Index size = this->getTrueSize();
	Eigen::Index count=0;
	/* ポインタをたどって探す */
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	const Eigen::Index total_size = matrix.nonZeros();
	for(Eigen::Index i = 0; i < size; i++) {
		const Eigen::Index num = (i==size-1 ? total_size : col_ptr[i+1]);
		for(Eigen::Index j = col_ptr[i]; j < num; j++) {
			if(i == gyo && row_ptr[count] == target_r) {
				result_retu = static_cast<slv_int>(j-col_ptr[i]);
				return true;
			}
			count++;
		}
	}
	result_retu = 0;
	return false;
}


/*//=======================================================
// ● 指定した列の非ゼロの行位置と値をvectorに書き出す 
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<DType>& row_val) const{
	const Eigen::Index size = this->getTrueSize();
	row_pos.clear();
	row_val.clear();
	Eigen::Index count=0;
	/* ポインタをたどって探す */
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const Eigen::Index total_size = matrix.nonZeros();
	for(Eigen::Index i = 0; i < size; i++) {
		const Eigen::Index num = (i==size-1 ? total_size : col_ptr[i+1]);
		for(Eigen::Index j = col_ptr[i]; j < num; j++) {
			if(row_ptr[count] == target) {
				row_pos.push_back(static_cast<slv_int>(i));
				row_val.push_back(val_ptr[count]);
			}
			count++;
		}
	}
}
/*//=======================================================
// ● 指定した行の非ゼロの列位置と値をvectorに書き出す
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<DType>& col_val)const{
	const Eigen::Index size = this->getTrueSize();

	col_pos.clear();
	col_val.clear();
	/* ポインタをたどって探す */
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const Eigen::Index total_size = matrix.nonZeros();
	Eigen::Index i = target;
	auto count = col_ptr[i];
	const auto num = (i==size-1 ? total_size : col_ptr[i+1]);
	for(Eigen::Index j = col_ptr[i]; j < num; j++) {
		col_pos.push_back(static_cast<slv_int>(row_ptr[count]));
		col_val.push_back(val_ptr[count]);
		count++;
	}
}

/*//=======================================================
// ● ほぼゼロになっている非ゼロ位置の削除
//=======================================================*/
/** 
 * @brief delete row and column where flag is True
 * @param flag boolean flag array
*/
template<typename DType>
void SparseMatTMPL<DType>::prune(double thres){
	Eigen::SparseMatrix<DType> tempMat = matrix.pruned(thres);
	matrix = std::move(tempMat);
}



/*//=======================================================
// ● 書き出し
//=======================================================*/
/** 
 * @brief write matrix to file
 * @param str filename
*/
template<typename DType>
void SparseMatTMPL<DType>::printMat(const std::string& str) const{
	const slv_int size = this->getSize();

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	this->getCols(start_pos1, end_pos1);
	auto col_ptr1 = matrix.innerIndexPtr();;//matrix.outerIndexPtr();
	auto val_ptr1 = matrix.valuePtr();

	std::fstream fp(str, std::ios::out);
	fp << "Column size : " << size << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		fp << end_pos1[i]-start_pos1[i] <<  std::endl;
	}

	fp << "Column : "  << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		slv_int col_size = end_pos1[i]-start_pos1[i];
		if(col_size == 0){
			fp << "-1" << std::endl;
		}else{
			for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
				fp << col_ptr1[j] << " ";
			}
			fp << std::endl;
		}
	}
	fp <<"Mat value : "  << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		slv_int col_size = end_pos1[i]-start_pos1[i];
		if(col_size > 0){
			for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
				fp <<  std::scientific << val_ptr1[j] << " ";
			}
			fp << std::endl;
		}else {
			fp << std::endl;
		}
	}
	fp.close();

}
/*//=======================================================
// ● ファイルから行列構成
//=======================================================*/
/** 
 * @brief construct matrix from file
 * @param filename filename
*/
template<typename DType>
void SparseMatTMPL<DType>::readMat(const std::string& filename){
	std::fstream fp(filename, std::ios::in);
	slv_int read_size;
	std::string temp_str1, temp_str2, temp_str3;
	fp >> temp_str1 >> temp_str2  >> temp_str3 >> read_size;	
	std::cout << "read size " << read_size << std::endl;
	slv_int size = read_size;
	
	SparseBuilderTMPL<DType> builder(size);

	int temp1, temp2;
	slv_int* temp_sizes = new slv_int[size];
	slv_int** temp_cols = new slv_int*[size];
	for(slv_int i = 0 ; i < size ; i++){
		fp >> temp_sizes[i];
		temp_cols[i] = new slv_int[temp_sizes[i]];
		std::cout << temp_sizes[i] << std::endl;
	}
	/**/
	/**/
	fp >> temp_str1 >> temp_str2;
	for(slv_int i = 0 ; i < size ; i++){
		if(temp_sizes[i] == 0){
			fp >> temp1;
		}else{
			for(slv_int j = 0; j < temp_sizes[i]; j++){
				fp >> temp2;
				temp_cols[i][j] = temp2;
				std::cout << i << ", " << j << " >" << temp2 << std::endl;
			}
		}
	}
	/**/
	/**/
	fp >> temp_str1 >> temp_str2  >> temp_str3;
	std::cout << temp_str1 << ", " << temp_str2  << ", " <<  temp_str3 << std::endl;;
	for(slv_int i = 0 ; i < size ; i++){
		for(slv_int j = 0; j < temp_sizes[i]; j++){
			DType tempD;
			fp >> tempD;
			std::cout << i << ", " << j << " >" << tempD << std::endl;
			builder.add(i, temp_cols[i][j], tempD);
		}
	}
	fp.close();
	//
	
	delete[] temp_sizes;	
	for(slv_int i = 0 ; i < size ; i++){
		delete[] temp_cols[i];
	}
	delete[] temp_cols;

	std::vector<Eigen::Triplet<DType>> tripletList;
	slv_int max_retu = builder.build(tripletList, false);
	/* 確定させる */
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(size, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
}


/* ===== ここが本体===== */
template class SparseMatTMPL<double>;
template class SparseMatTMPL<dcomplex>;





/* end of namespace */
};
};
