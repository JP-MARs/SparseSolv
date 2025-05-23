﻿
#include "SparseMat.hpp"
#include "SparseMatOperators.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{

/*//=======================================================
// ● コンストラクタ
//=======================================================*/
/* Eigen疎行列指定 */
SparseMat::SparseMat(const Eigen::SparseMatrix<double, Eigen::RowMajor>& Mat0){
	matrix.initialize(Mat0);
}
/* サイズ指定 */
SparseMat::SparseMat(slv_int x){
	matrix.tempInitialize(x);
}
/* 作成済み疎行列データから初期化 */
SparseMat::SparseMat(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<double>& vals){
	matrix.initialize(n_nonZero, rows, cols, vals);
}


/*//=======================================================
// ● コピーコンストラクタたち
//=======================================================*/
/* 同クラスコピー */
SparseMat::SparseMat(const SparseMat& mat){
	matrix = mat.matrix;
}
/* 同クラスムーブ */
SparseMat::SparseMat(SparseMat&& mat) noexcept{
	matrix = std::move(mat.matrix);
}
/* 行列クラスコピー */
SparseMat::SparseMat(const SparseMatBaseD& mat) {
	matrix = mat;
}
/* 行列クラスムーブ */
SparseMat::SparseMat(SparseMatBaseD&& mat) noexcept {
	matrix =std::move(mat);
}

/*//=======================================================
// ● 代入オペレータ
//=======================================================*/
SparseMat& SparseMat::operator=(const SparseMat& mat){
	matrix = mat.matrix;
	return *this;
}
/*//=======================================================
// ● ムーブオペレータ
//=======================================================*/
SparseMat& SparseMat::operator=(SparseMat&& mat) noexcept{
	matrix = std::move(mat.matrix);
	return *this;
}

/*
＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
 ■ オペレータ群
＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝
*/

/*//=======================================================
// ● 行列ベクトル積
//=======================================================*/
double* SparseMat::operator*(const double* vec) const{
	double* ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseD, double, double, double>(&ans, matrix, vec);
	return ans;
}
dcomplex* SparseMat::operator*(const dcomplex* vec) const{
	dcomplex* ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseD, double, dcomplex, dcomplex>(&ans, matrix, vec);
	return ans;
}
Eigen::VectorXd SparseMat::operator*(const Eigen::VectorXd& vec) const{
	Eigen::VectorXd ans = matrix.matrix * vec;
	return ans;
}
Eigen::VectorXcd SparseMat::operator*(const Eigen::VectorXcd& vec) const{
	Eigen::VectorXcd ans = matrix.matrix * vec;
	return ans;
}
std::vector<double> SparseMat::operator*(const std::vector<double>& vec) const{
	std::vector<double> ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseD, double>(ans, matrix, vec);
	return ans;
}
std::vector<dcomplex> SparseMat::operator*(const std::vector<dcomplex>& vec) const{
	std::vector<dcomplex> ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseD>(ans, matrix, vec);
	return ans;
}

/*//=======================================================
// ● 行列スカラ積
//=======================================================*/
SparseMat SparseMat::operator*(const double x) const {
	SparseMatBaseD mat_ori(matrix);
	mat_ori *= x;
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMat::operator*(const dcomplex x) const {
	SparseMatBaseC mat_oriC;
	SparseMatOperators::copyDtoC(mat_oriC, matrix);
	mat_oriC *= x;
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_oriC));
	return mat_ans;
}

/*//=======================================================
// ● 行列-行列積
//=======================================================*/
SparseMat SparseMat::operator*(const SparseMat& mat) const {
	Eigen::SparseMatrix<double, Eigen::RowMajor> tempMat = ( (matrix.matrix) * (mat.matrix.matrix) ).pruned();
	SparseMatBaseD mat_ori(std::move(tempMat));
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMat::operator*(const SparseMatC& mat) const {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> tempMat = ((matrix.matrix) * (mat.matrix.matrix)).pruned();
	SparseMatBaseC mat_ori(std::move(tempMat));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 行列-行列和
//=======================================================*/
SparseMat SparseMat::operator+(const SparseMat& mat) const {
	Eigen::SparseMatrix<double, Eigen::RowMajor> tempMat = (matrix.matrix) + (mat.matrix.matrix);
	SparseMatBaseD mat_ori(std::move(tempMat));
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMat::operator+(const SparseMatC& mat) const {
	SparseMatBaseC mat_ori;
	SparseMatOperators::plus_operators<SparseMatBaseC, SparseMatBaseD, SparseMatBaseC, dcomplex, double>(mat_ori, mat.matrix, matrix, 1.0, 1.0, 0, 0);
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}



/*================================*/
/*================================*/
/*================================*/

/*//=======================================================
// ● 転置
//=======================================================*/
SparseMat SparseMat::trans() const {
	SparseMatBaseD mat_ori(matrix.matrix.transpose());
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}
/*//=======================================================
// ● 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、Bに渡す
//=======================================================*/
SparseMat SparseMat::makeSubMat(slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b) {
	SparseMatBaseD mat_ori;
	SparseMatOperators::makeSubMat<SparseMatBaseD, double>(mat_ori, matrix, range1a, range1b, range2a, range2b);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;

}
/*//=======================================================
// ● 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する
//=======================================================*/
void SparseMat::MatDiv(SparseMat& matK11, SparseMat& matK12, slv_int rangeA, slv_int rangeB) {
	SparseMatBaseD mat_ori1;
	SparseMatBaseD mat_ori2;
	SparseMatOperators::MatDiv<SparseMatBaseD, double>(matrix, mat_ori1, mat_ori2, rangeA, rangeB);
	/* 結果を本体行列にムーブし、終わる */
	matK11.matrix = std::move(mat_ori1);
	matK12.matrix = std::move(mat_ori2);
}

/*//=======================================================
// ● 下三角行列の取得
//=======================================================*/
SparseMat SparseMat::getMatLower() const {
	SparseMatBaseD mat_ori;
	SparseMatOperators::getMatLower<SparseMatBaseD, double>(mat_ori, matrix);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}
/*//=======================================================
// ● 上三角行列の取得
//=======================================================*/
SparseMat SparseMat::getMatUpper() const {
	SparseMatBaseD mat_ori;
	SparseMatOperators::getMatUpper<SparseMatBaseD, double>(mat_ori, matrix);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 逆行列（密アルゴリズム・大行列でやるとメモリが吹き飛ぶ！）
//=======================================================*/
SparseMat SparseMat::inv() const{
	SparseMatBaseD mat_ori;
	SparseMatOperators::MatInv<SparseMatBaseD, double, Eigen::MatrixXd>(mat_ori, matrix);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 行列値を丸めて整数にする
//=======================================================*/
void SparseMat::round(){
	/* ポインタをたどって探す */
	auto val_ptr = matrix.matrix.valuePtr();
	const slv_int total_size = matrix.matrix.nonZeros();
	for(slv_int i = 0; i < total_size; i++) {
		double val = val_ptr[i];
		if(val >= 0.0){
			val += 0.5;
		}else{
			val -= 0.5;
		}
		val_ptr[i] = (int)(val);
	}
}

/*//=======================================================
// ● 疑似逆行列用の行列（A^T*A+epsI）を作成
//=======================================================*/
SparseMat SparseMat::makePrsdInv(double eps) const{
	SparseMatBaseD mat_ori;
	SparseMatOperators::AtA_eps<SparseMatBaseD, double>(mat_ori, matrix, eps);
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}


/*//=======================================================
// ● 不完全コレスキー分解
//=======================================================*/
SparseMat SparseMat::IC_decomp(double* diagD, const double accela) const{
	SparseMatBaseD mat_ori;
	SparseMatOperators::IC_decomp<SparseMatBaseD, double, double>(mat_ori, matrix, diagD, accela);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 対角スケーリング
//=======================================================*/
SparseMat SparseMat::diagScaling(double* trans_vec, const double* ori_vec) const{
	SparseMatBaseD mat_ori;
	SparseMatOperators::diagScaling(mat_ori, matrix, trans_vec, ori_vec);
	/* 結果を本体行列にムーブし、終わる */
	SparseMat mat_ans(std::move(mat_ori));
	return mat_ans;
}


/* end of namespace */
};
