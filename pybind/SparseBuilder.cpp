/*
* <SparseSolv>
* Copyright (c) 2026 Takahiro Sato
*
* This source code is licensed under the MPL2 License.
* See the LICENSE file in the root directory for details.
*/

/**
 * @file SparseBuilder.cpp
 * @brief Sparse Matrix Builder Template Class
 * このファイルはテンプレートクラス SparseBuilder を提供します。
 * このクラスの実体を実装してtemplate公開します。
 */

#include <SparseSolv/SparseBuilder.hpp>

/* オリジナル名前空間(JPMARsライブラリ) */
namespace JPMRspace{
/* 疎行列ソルバ用オリジナル名前空間 */
namespace SparseSolv{

/*
//=======================================================
// ■ スパース行列構築補助クラステンプレート（自作形式）
//=======================================================*/



/*//=======================================================
// ● 一時行列を作成
//=======================================================*/
/** 
 * @brief initialize temp. data
 * @param size0 row size
*/
template<typename DType>
void SparseBuilderTMPL<DType>::initialize(slv_int ss){
	size = ss;
	this->initialize();
}


/*//=======================================================
// ● 一時行列を作成
//=======================================================*/
/** 
 * @brief initialize temp. data (main)
*/
template<typename DType>
void SparseBuilderTMPL<DType>::initialize(){
	tempMat.assign(size, std::map<slv_int, DType>());
}

/*//=======================================================
// ● 構築済み部分の値を０に再セット
//=======================================================*/
/** 
 * @brief reset values to zeros
*/
template<typename DType>
void SparseBuilderTMPL<DType>::resetMat(){
	for(slv_int i = 0 ; i < size ; i++){
		/* kv は pair<const slv_int, DType>*/
		for(auto& kv : tempMat[i]){   
			/* 値をゼロ相当で初期化*/
			kv.second = DType{};      
		}
	}
}

/*//=======================================================
// ● 一次配列の状況をprint
//=======================================================*/
/** 
 * @brief print Temp.Mat
*/
template<typename DType>
void SparseBuilderTMPL<DType>::printTempMat()const{
	for(slv_int i = 0 ; i < size ; i++){
		std::cout << i << "-th size = " << tempMat[i].size() << " ----- ";
		for(const auto& kv : tempMat[i]){   
			/* 値をゼロ相当で初期化*/
			std::cout << kv.first << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << "----" << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		for(const auto& kv : tempMat[i]){   
			/* 値をゼロ相当で初期化*/
			std::cout << kv.second << ", ";
		}
		std::cout << std::endl;
	}
}

/*//=======================================================
// ● 一次配列からtripletを作り、一次配列を削除
//=======================================================*/
/** 
 * @brief build Eigen::triplet from tempMat
 * @param to_square true: make square matrix by adding zeros
 * @return max col number
*/
template<typename DType>
slv_int SparseBuilderTMPL<DType>::build(std::vector<Eigen::Triplet<DType>>& tripletList, bool to_square){
	slv_int max_retu = 0;
	tripletList.clear();
	/* 予めサイズを確保しておく（適当な値） */
	tripletList.reserve(size*5); 
	/*  tempMatの内容をtripletListに変換していく */
	for(slv_int i = 0; i < size; i++) {
		const auto ss = tempMat[i].size();
		if(ss == 0) {
			continue;
		}
		for(auto& itr : tempMat[i]) {
			const slv_int pos = itr.first;
			const DType val = itr.second;
			//Eigen::Triplet<DType> triplet = Eigen::Triplet<DType>(i, pos, val);
			//tripletList.push_back( std::move(triplet) );
			tripletList.emplace_back(i, pos, val);
			if(pos > max_retu){
				max_retu = pos;
			}
		}
		std::map<slv_int, DType> empty; empty.clear();
		tempMat[i].swap(empty);
	}
	max_retu++;
	/* toSquare=true なら、行数の方が大きいなら正方サイズにする */
	/* toSquare=false なら、長方形のまま */
	if(to_square){		
		if(max_retu < size){
			tripletList.push_back(Eigen::Triplet<DType>(size-1, size-1, 0.0));
			max_retu = size;
		}	
	}
	return max_retu;
}



/* ===== ここが本体===== */
template class SparseBuilderTMPL<double>;
template class SparseBuilderTMPL<dcomplex>;





/* end of namespace */
};
};
