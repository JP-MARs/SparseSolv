
#include "SparseMatC.hpp"
#include "SparseMatOperators.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*//=======================================================
// ● コンストラクタ
//=======================================================*/
SparseMatC::SparseMatC(){
	matrix = new SparseMatBaseC;
}
/* サイズ指定 */
SparseMatC::SparseMatC(slv_int x){
	matrix = new SparseMatBaseC(x);
}
/* 作成済み疎行列データから初期化 */
SparseMatC::SparseMatC(slv_int n_nonZero, const std::vector<slv_int>& rows, const std::vector<slv_int>& cols, const std::vector<dcomplex>& vals){
	matrix = new SparseMatBaseC(n_nonZero, rows, cols, vals);
}

/*//=======================================================
// ● デストラクタ
//=======================================================*/
SparseMatC::~SparseMatC(){
	delete matrix;
}

/*//=======================================================
// ● コピーコンストラクタたち
//=======================================================*/
/* 同クラスコピー */
SparseMatC::SparseMatC(const SparseMatC& mat){
	matrix = new SparseMatBaseC(*mat.matrix);
}
/* 同クラスムーブ */
SparseMatC::SparseMatC(SparseMatC&& mat) noexcept{
	matrix = new SparseMatBaseC(std::move(*mat.matrix));
	mat.matrix = nullptr;
}
/* 行列クラスコピー */
SparseMatC::SparseMatC(const SparseMatBaseC& mat) {
	matrix = new SparseMatBaseC(mat);
}
/* 行列クラスムーブ */
SparseMatC::SparseMatC(SparseMatBaseC&& mat) noexcept {
	matrix = new SparseMatBaseC(std::move(mat));
}

/*//=======================================================
// ● 代入オペレータ
//=======================================================*/
SparseMatC& SparseMatC::operator=(const SparseMatC& mat){
	if(matrix != nullptr){
		delete matrix;
	}
	*matrix = *mat.matrix;
	return *this;
}
/*//=======================================================
// ● ムーブオペレータ
//=======================================================*/
SparseMatC& SparseMatC::operator=(SparseMatC&& mat) noexcept{
	if(matrix != nullptr){
		delete matrix;
	}
	*matrix = std::move(*mat.matrix);
	mat.matrix = nullptr;
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
dcomplex* SparseMatC::operator*(const double* vec) const{
	dcomplex* ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseC, dcomplex, double, dcomplex>(&ans, *(this->matrix), vec);
	return ans;
	//dcomplex* ans = nullptr;
	//return ans;
}
dcomplex* SparseMatC::operator*(const dcomplex* vec) const{
	dcomplex* ans;
	SparseMatOperators::mat_vec_product<SparseMatBaseC, dcomplex, dcomplex, dcomplex>(&ans, *(this->matrix), vec);
	return ans;
}
Eigen::VectorXcd SparseMatC::operator*(const Eigen::VectorXd& vec) const{
	Eigen::VectorXcd ans = matrix->matrix * vec;
	return ans;
}
Eigen::VectorXcd SparseMatC::operator*(const Eigen::VectorXcd& vec) const{
	Eigen::VectorXcd ans = matrix->matrix * vec;
	return ans;
}


/*//=======================================================
// ● 行列スカラ積
//=======================================================*/
SparseMatC SparseMatC::operator*(const double x) const {
	SparseMatBaseC mat_ori(*(this->matrix));
	mat_ori *= x;
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMatC::operator*(const dcomplex x) const {
	SparseMatBaseC mat_oriC(*(this->matrix));
	mat_oriC *= x;
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_oriC));
	return mat_ans;
}

/*//=======================================================
// ● 行列-行列積
//=======================================================*/
SparseMatC SparseMatC::operator*(const SparseMat& mat) const {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> tempMat = ((matrix->matrix) * (mat.matrix->matrix)).pruned();
	SparseMatBaseC mat_ori(std::move(tempMat));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMatC::operator*(const SparseMatC& mat) const {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> tempMat = ((matrix->matrix) * (mat.matrix->matrix)).pruned();
	SparseMatBaseC mat_ori(std::move(tempMat));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 行列-行列和
//=======================================================*/
SparseMatC SparseMatC::operator+(const SparseMat& mat) const {
	SparseMatBaseC mat_ori;
	SparseMatOperators::plus_operators<SparseMatBaseC, SparseMatBaseD, SparseMatBaseC, dcomplex, double>(mat_ori, *matrix, *(mat.matrix), 1.0, 1.0, 0, 0);
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
SparseMatC SparseMatC::operator+(const SparseMatC& mat) const {
	Eigen::SparseMatrix<dcomplex, Eigen::RowMajor> tempMat = ((matrix->matrix) + (mat.matrix->matrix));
	SparseMatBaseC mat_ori(std::move(tempMat));
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
SparseMatC SparseMatC::trans() const {
	SparseMatBaseC mat_ori(matrix->matrix.transpose());
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
/*//=======================================================
// ● 範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、Bに渡す
//=======================================================*/
SparseMatC SparseMatC::makeSubMat(slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b) {
	SparseMatBaseC mat_ori;
	SparseMatOperators::makeSubMat<SparseMatBaseC, dcomplex>(mat_ori, *(this->matrix), range1a, range1b, range2a, range2b);
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
/*//=======================================================
// ● 列rangeBを境に、行列を2つに分ける。行rangeAより下は削除する
//=======================================================*/
void SparseMatC::MatDiv(SparseMatC& matK11, SparseMatC& matK12, slv_int rangeA, slv_int rangeB) {
	SparseMatBaseC mat_ori1;
	SparseMatBaseC mat_ori2;
	SparseMatOperators::MatDiv<SparseMatBaseC, dcomplex>(*(this->matrix), mat_ori1, mat_ori2, rangeA, rangeB);
	/* 結果を本体行列にムーブし、終わる */
	*(matK11.matrix) = std::move(mat_ori1);
	*(matK12.matrix) = std::move(mat_ori2);
}

/*//=======================================================
// ● 下三角行列の取得
//=======================================================*/
SparseMatC SparseMatC::getMatLower() const {
	SparseMatBaseC mat_ori;
	SparseMatOperators::getMatLower<SparseMatBaseC, dcomplex>(mat_ori, *(this->matrix));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}
/*//=======================================================
// ● 上三角行列の取得
//=======================================================*/
SparseMatC SparseMatC::getMatUpper() const {
	SparseMatBaseC mat_ori;
	SparseMatOperators::getMatUpper<SparseMatBaseC, dcomplex>(mat_ori, *(this->matrix));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 逆行列（密アルゴリズム・大行列でやるとメモリが吹き飛ぶ！）
//=======================================================*/
SparseMatC SparseMatC::inv() const{
	SparseMatBaseC mat_ori;
	SparseMatOperators::MatInv<SparseMatBaseC, dcomplex, Eigen::MatrixXcd>(mat_ori, *(this->matrix));
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;

}

/*//=======================================================
// ● 疑似逆行列用の行列（A^T*A+epsI）を作成
//=======================================================*/
SparseMatC SparseMatC::makePrsdInv(double eps) const{
	SparseMatBaseC mat_ori;
	SparseMatOperators::AtA_eps<SparseMatBaseC, dcomplex>(mat_ori, *(this->matrix), eps);
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}

/*//=======================================================
// ● 不完全コレスキー分解
//=======================================================*/
SparseMatC SparseMatC::IC_decomp(dcomplex* diagD, const double accela) const{
	SparseMatBaseC mat_ori;
	SparseMatOperators::IC_decomp<SparseMatBaseC, dcomplex>(mat_ori, *(this->matrix), diagD, accela);
	/* 結果を本体行列にムーブし、終わる */
	SparseMatC mat_ans(std::move(mat_ori));
	return mat_ans;
}



/* end of namespace */
};
