/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file SprsTrans.cpp
 * @brief Matrix transform Operators for Sparse Matrix implementation
 * 
 * このファイルは SparseMatTMPL の行列構造を変換する処理を実装します。
 */

#include <SparseSolv/SprsTrans.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv {
/* 疎行列変換処理用の専用名前空間 */
namespace SprsTrans {

/* 疎行列変換処理用内部空間間 */
/**
* @namespace SprsTrans::detai
* @brief Original namespace for internal use of SprsTrans
*/
namespace detail {

/*//=======================================================
// ● 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作るサブルーチン
//=======================================================*/
/** 
 * @brief make sub matrix subrutin
 * @param r1a [r1a, r1b]×[r2a, r2b]
 * @param r1b [r1a, r1b]×[r2a, r2b]
 * @param r2a [r1a, r1b]×[r2a, r2b]
 * @param r2b [r1a, r1b]×[r2a, r2b]
*/
template<typename DType>
void makeSubMat_sub(SparseMatTMPL<DType>& resultMat, const SparseMatTMPL<DType>& matA, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b){
	const slv_int the_size = matA.getSize();
	const slv_int new_size = range1b-range1a + 1;
	
	/* 構築準備 */
	SparseBuilderTMPL<DType> builder(new_size);

	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const auto col_ptr1 = matA.getColPtr();
	const auto val_ptr1 = matA.getValuePtr();
	for(slv_int i = 0 ; i < the_size ; i++){
		if(range1a <= i && i <= range1b){
			const slv_int col_size = end_pos1[i];
			for(slv_int j = start_pos1[i] ; j < col_size ; j++){
				slv_int col = col_ptr1[j];
				if(range2a <= col && col <= range2b){
					DType val = val_ptr1[j];
					builder.add(i-range1a, col-range2a, val);
				}
			}
		}
	}
	/* 構築 */
	resultMat.build(builder, false);
}

/*//=======================================================
// ● 列rangeBを境に、行列を2つに分ける。行rangeAより下は分けた行列には含めない サブルーチン
//       K11←B→K12
//    (A11) ..|... (A1n)
//      .     |      .
//rangeA______|_______
// ↓無視     |      .
//    (Am1) ..|... (Amn)
//=======================================================*/
/** 
 * @brief make sub matrix subturine for dividing matrix into 2 by col rangeB (K11←B→K12, rows below rangeA are not included in K11 and K12)
 * @param rangeA row number deleted from this
 * @param rangeB col number deleted from this
*/
template<typename DType>
void MatDiv_sub(SparseMatTMPL<DType>& matK11, SparseMatTMPL<DType>& matK12, const SparseMatTMPL<DType>& matA, slv_int rangeA, slv_int rangeB){
	const slv_int the_size = matA.getSize();
	const slv_int gyo_limit = rangeA+1;
	const slv_int min_gyo = (gyo_limit >= the_size ? the_size : gyo_limit);
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	const auto col_ptr1 = matA.getColPtr();
	const auto val_ptr1 = matA.getValuePtr();

	/* 構築準備 */
	SparseBuilderTMPL<DType> builder11(min_gyo);
	SparseBuilderTMPL<DType> builder12(min_gyo);

	for(slv_int i = 0 ; i < min_gyo ; i++){
		const slv_int col_size = end_pos1[i];
		for(slv_int j = start_pos1[i] ; j < col_size ; j++){
			slv_int col = col_ptr1[j];
			DType val = val_ptr1[j];
			if(col <= rangeB){
				builder11.add(i, col, val);
			}else{
				builder12.add(i, col, val);
			}
		}
	}
	/* 構築 */
	matK11.build(builder11, false);
	matK12.build(builder12, false);
}


/*//=======================================================
// ● フラグがTrueの位置の列・行のデータを削除（FEM用） - サブルーチン
//=======================================================*/
/** 
 * @brief delete row and column where flag is True (subroutine)
 * @param flag boolean flag array
*/
template<typename DType>
void delFlagPosition_sub(SparseMatTMPL<DType>& matA, const bool* flag){
	slv_int size = matA.getSize();

	/* 行列位置を取得 */
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	/* 削除分のID再計算 */
	std::map<slv_int, slv_int> id_map;
	slv_int counter=0;
	for(slv_int i = 0; i < size; i++){
		if(flag[i]){
			continue;
		}
		id_map[i] = counter;
		counter++;
	}
	slv_int max_retu = 0;
	std::vector<Eigen::Triplet<DType>> tripletList;
	tripletList.reserve(matA.mat().nonZeros());
	/* 再度探索 */
	for(slv_int i = 0 ; i < size ; i++){
		/* フラグがonなら無視 */
		if(flag[i]){
			continue;
		}
		for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
			slv_int col = col_ptr1[j];
			/* フラグがonなら無視 */
			if(flag[col]){
				continue;
			}
			/* フラグがoffなら加える */
			DType val = val_ptr1[j];
			//tripletList.push_back(Eigen::Triplet<DType>(id_map[i], id_map[col], val));
			tripletList.emplace_back(id_map[i], id_map[col], val);
			if(max_retu < id_map[col]){
				max_retu = id_map[col];
			}
		}
	}
	/* 再セット */
	max_retu++;
	Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(counter, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	/* 確定 */
	matA.setMatrix(std::move(matrix));
}



/*//=======================================================
// ● フラグがTrueの位置の列・行のデータをゼロにして対角だけ１に（FEM用）
//=======================================================*/
/** 
 * @brief set one on diagonal and zero elsewhere where flag is True
 * @param flag boolean flag array
*/
template<typename DType>
void DiagFlagPosition_sub(SparseMatTMPL<DType>& matA, const bool* flag){
	slv_int size = matA.getSize();

	/* 行列位置を取得 */
	std::vector<slv_int> start_pos1;
	std::vector<slv_int> end_pos1;
	matA.getCols(start_pos1, end_pos1);
	auto col_ptr1 = matA.getColPtr();
	auto val_ptr1 = matA.getValuePtr();

	slv_int max_retu = 0;
	std::vector<Eigen::Triplet<DType>> tripletList;
	tripletList.reserve(matA.mat().nonZeros());
	/* 再度探索 */
	for(slv_int i = 0 ; i < size ; i++){
		/* フラグがonなら、対角だけ探して無視 */
		if(flag[i]){
			for(slv_int j = start_pos1[i]; j < end_pos1[i]; j++){
				slv_int col = col_ptr1[j];
				if(col == i){
					//tripletList.push_back(Eigen::Triplet<DType>(i, i, 1.0));
					tripletList.emplace_back(i, i, 1.0);
					if(max_retu < i){
						max_retu = i;
					}
					break;
				}
			}
			continue;
		}
		for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
			slv_int col = col_ptr1[j];
			/* フラグがonで */
			if(flag[col]){
				continue;
			}
			/* フラグがoffなら加える */
			DType val = val_ptr1[j];
			//tripletList.push_back(Eigen::Triplet<DType>(i, col, val));
			tripletList.emplace_back(i, col, val);
			if(max_retu < col){
				max_retu = col;
			}
		}
	}
	/* 再セット */
	max_retu++;
	Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(size, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());
	/*  */
	matA.setMatrix(std::move(matrix));	
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
// ● 行列の値を四捨五入）
//=======================================================*/
/** 
 * @brief matrix rounding
*/
void round(SparseMat& matA) {
	double* val_ptr = matA.getValuePtr();
	const Eigen::Index total_size = matA.mat().nonZeros();
	/* 非ゼロ位置をすべて見て、それぞれ四捨五入 */
	for(Eigen::Index i = 0; i < total_size; ++i) {
		val_ptr[i] = std::round(val_ptr[i]);
	}
	matA.mat().prune(1.0e-12); /* 四捨五入してゼロになった位置を削除 */
}

/*//=======================================================
// ● 行列Aにa1をかける
//=======================================================*/
/** 
 * @brief matA = a*matA
*/
void SelfDot(SparseMat& matA, double a1) {
	double* val_ptr = matA.getValuePtr();
	const Eigen::Index total_size = matA.mat().nonZeros();
	/* 非ゼロ位置をすべて見て、それぞれにa1をかける */
	for(Eigen::Index i = 0; i < total_size; ++i) {
		val_ptr[i] = a1 * val_ptr[i];
	}
}
void SelfDot(SparseMatC& matA, double a1) {
	dcomplex* val_ptr = matA.getValuePtr();
	const Eigen::Index total_size = matA.mat().nonZeros();
	/* 非ゼロ位置をすべて見て、それぞれにa1をかける */
	for(Eigen::Index i = 0; i < total_size; ++i) {
		val_ptr[i] = a1 * val_ptr[i];
	}
}

/*//=======================================================
// ● ほぼゼロになっている非ゼロ位置の削除
//=======================================================*/
/** 
 * @brief prune matA's non-zero entries that are almost zero
*/
void prune(SparseMat& matA, double thres) {
	matA.mat().prune(thres);
	matA.mat().makeCompressed();
}


/*//=======================================================
// ● 転置
//=======================================================*/
/** 
 * @brief Matrix transpose
*/
SparseMat trans(const SparseMat& matA) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = matA.mat().transpose();
	return( SparseMat( std::move(out) ) );
}
SparseMatC trans(const SparseMatC& matA) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = matA.mat().transpose();
	return( SparseMatC( std::move(out) ) );
}


/*//=======================================================
// ● 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作る
//=======================================================*/
/** 
 * @brief make sub matrix
 * @param r1a [r1a, r1b]×[r2a, r2b]
 * @param r1b [r1a, r1b]×[r2a, r2b]
 * @param r2a [r1a, r1b]×[r2a, r2b]
 * @param r2b [r1a, r1b]×[r2a, r2b]
*/
SparseMat makeSubMat(const SparseMat& matA, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b) {
	SparseMat out;
	detail::makeSubMat_sub<double>(out, matA, range1a, range1b, range2a, range2b);
	return out;
}
SparseMatC makeSubMat(const SparseMatC& matA, slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b) {
	SparseMatC out;
	detail::makeSubMat_sub<dcomplex>(out, matA, range1a, range1b, range2a, range2b);
	return out;
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
/** 
 * @brief make sub matrix subturine for dividing matrix into 2 by col rangeB (K11←B→K12, rows below rangeA are not included in K11 and K12)
 * @param rangeA row number deleted from this
 * @param rangeB col number deleted from this
*/
void MatDiv(SparseMat& matK11, SparseMat& matK12, const SparseMat& matA, slv_int rangeA, slv_int rangeB) {
	detail::MatDiv_sub<double>(matK11, matK12, matA, rangeA, rangeB);
}
void MatDiv(SparseMatC& matK11, SparseMatC& matK12, const SparseMatC& matA, slv_int rangeA, slv_int rangeB) {
	detail::MatDiv_sub<dcomplex>(matK11, matK12, matA, rangeA, rangeB);
}


/*//=======================================================
// ● 下三角行列の取得
//=======================================================*/
/** 
 * @brief lower tri matrix
*/
SparseMat MatLower(const SparseMat& matA) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = (matA.mat()).triangularView<Eigen::Lower>();
	return( SparseMat(std::move(out)) );
}
SparseMatC MatLower(const SparseMatC& matA) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = (matA.mat()).triangularView<Eigen::Lower>();
	return( SparseMatC(std::move(out)) );
}

/*//=======================================================
// ● 上三角行列の取得
//=======================================================*/
/** 
 * @brief upper tri matrix
*/
SparseMat MatUpper(const SparseMat& matA) {
	Eigen::SparseMatrix<double, Eigen::RowMajor> out = (matA.mat()).triangularView<Eigen::Upper>();
	return( SparseMat(std::move(out)) );
}
SparseMatC MatUpper(const SparseMatC& matA) {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> out = (matA.mat()).triangularView<Eigen::Upper>();
	return( SparseMatC(std::move(out)) );
}



/*//=======================================================
// ● フラグがTrueの位置の列・行のデータを削除（FEM用）
//=======================================================*/
/** 
 * @brief delete row and column where flag is True
 * @param flag boolean flag array
*/
void DelFlagPosition(SparseMat& matA, const bool* flag) {
	detail::delFlagPosition_sub<double>(matA, flag);

}
void DelFlagPosition(SparseMatC& matA, const bool* flag) {
	detail::delFlagPosition_sub<dcomplex>(matA, flag);
}

/*//=======================================================
// ● フラグがTrueの位置の列・行のデータをゼロにして対角だけ１に（FEM用）
//=======================================================*/
/** 
 * @brief set one on diagonal and zero elsewhere where flag is True
 * @param flag boolean flag array
*/
void DiagFlagPosition(SparseMat& matA, const bool* flag){
	detail::DiagFlagPosition_sub<double>(matA, flag);
}
void DiagFlagPosition(SparseMatC& matA, const bool* flag){
	detail::DiagFlagPosition_sub<dcomplex>(matA, flag);
}







/* end of namespace */
};
};
};

