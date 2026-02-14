<div style="text-align: justify;">

# Introduction
## 基本クラス
本ライブラリを利用する際に重要なクラスは以下の通りです。
* SparseMat: 実数疎行列のクラス
* SparseMatC: 複素実数疎行列のクラス
* MatSolvers: 線形ソルバクラス

## 疎行列の構築方法
&ensp;疎行列は、まず行数を指定して実体を作成します。この段階ではc++のmap[]を用いて、一時保存状態になります。
疎行列の値は、addメソッドを用いて代入していきます。FEMなどでは各位置の値が一度に決まることはないと思いますので、addを使って順次足していってください。
仕様上、最初に指定した以上の行に代入するとエラーになるので注意してください。<br>
&ensp;疎行列の構築が終わったら、fixメソッドを使うことで疎行列の配置を確定させます。具体的には、Eigenの疎行列を作成します。
行列の確定後（fixを呼び出した後）は、四則演算をEigenの機能を用いることで高速に行えます。
（なお、確定後もaddを使って行列の更新は可能です。）

### 疎行列の構築サンプル
&ensp;疎行列の構築サンプルです。
```
double matA[5][5] ={
	{ 61,    0,   3.1,  -0.1,   1},
	{  0,   42.5, 2.0,   0.5,   4},
	{ 3.1,  2.0, 10.5,     0,   0},
	{-0.1, 0.5,     0,  80.4, 1.5},
	{   1,   4,     0,   1.5,  1000.0}
};
size = 5
//疎行列を定義
SRLfem::SparseMat matAs(size);
/* 行列A,B,Cの値をセット */
for(int i = 0 ; i < size ; i++){
	for(int j = 0 ; j < size ; j++){
		if( fabs(matA[i][j]) > 1.0e-12 ){
			matAs.add(i, j, matA[i][j]);
		}
	}
}
//疎行列位置を確定
matAs1.fix();
```

なお、fixを呼ぶと、現在の要素位置に従ってEigen疎行列を作成します。
つまり、現在、長方形になっている場合はそのまま長方形で確定させます。
数値計算上、正方形にしたい場合は、fixに引数trueを追加してください。
```
matAs1.fix(true);
```
こうすることで、正方形になるようにゼロを代入してEigen疎行列を構築します。


### 疎行列の四則演算
&ensp;fixで確定後は、疎行列の演算を実行可能です。
（なお、fix前の行列に演算しようと思うと未定義動作となりバグります。）
```
SRLfem::SparseMat matAs(5);
SRLfem::SparseMat matBs(5);
......

SRLfem::SparseMat matCs1 = matAs1 + matBs;
SRLfem::SparseMat matCs2 = matAs1 * matBs;

///...etc


//行列・ベクトル積も対応可能

double vec[5] ={1,2,3,4,5};
double* vec2 = matAs*vec;

vector<double> vec2{1,2,3,4,5};
vector<double> ans;
ans = matAs*vec2;

Eigen::VectorXd vec3
vec3 << 1,2,3,4,5;
Eigen::VectorXd ans2 = matAs*vec3;

```

&ensp;ほかにも、転置、部分抽出、上三角、下三角取得、など、数値計算に必要そうな処理は大体実装済みです。
SparseMat.hppの、メソッド定義をご覧ください。



## 線形ソルバの利用
&ensp;線形方程式Ax=bを解くソルバを複数実装しています。
特に電磁気のFEMで有効とされる加速係数付きICCG法を実装している、数少ないOSSと思います。
本OSS独自実装のソルバは下記のとおりです。

* 加速係数付きICCG法（実数・複素）
* 加速係数付きMRTR法（実数・複素）
* Eisenstatの方法を導入した前処理付きMRTR法（実数・複素）(SGS-MRTR)

&ensp;なお、現状は、**大規模並列化には対応していません**。Eigenの並列計算機能を用いた小規模並列のみです。数コア程度なら並列化で多少早くなりますが、4,6コア以上を使ってもほとんど効果がないと思います。
デフォルトではEigenの並列化がonになっていると思いますので、Eigen機能を使って並列数を調整ください。
不要な場合はEigenの並列化機能をOFFすればソルバ内の並列化もOFFになります。
```
Eigen::setNbThreads(4);
int n = Eigen::nbThreads( );
cout << "Thred Eigen " << n << endl;
```

### ソルバの実行サンプル
&ensp;独自実装のソルバはMatSolversクラスを使って実行できます。
収束判定を設定するために内部状態を持っているので、インスタンス化して利用してください。
以下は、ICCG法の実行例です。
```
SRLfem::MatSolvers solvers;

//@param1: 行数
//@param2: 収束判定
//@param3: 最大反復数
//@param4: 加速係数（1以上だとその値で固定。負の値にすると、その絶対値を超えない範囲で自動計算）
//@param5: 疎行列A
//@param6: 右辺ベクトルb
//@param7: 解ベクトルx
//@param8: ソルバ開始前の解ベクトルゼロ初期化有無。trueでxをゼロで初期化。falseだと何もしない(省略時はfalse)。
bool bl = solvers.solveICCG(size, 1.0e-8, 10000, 1.05, matAs, vec, results, true);
```

&ensp;MatSolverクラスは、収束判定などにいくつかの内部変数を持っています。
代表的なものは以下。
* **setSaveLog(true)**: 残差のログを保存します。デフォルトはOFF。
* **setSaveBest(true)**: 最良解を保存するか。デフォルトはOFF。ONのとき、発散時には保存している最良解を解ベクトルxに代入します。
* **setConvNormalizeType(int x)**: 収束判定の正規化タイプ(0:右辺、1:初期残差、2:外部指定)
* **setDiagScale(true)**: 対角スケーリングの有無をセット。デフォルトはOFFです。

残差を記録している場合は、以下のように取得できます。
```
SRLfem::MatSolvers solvers;
vector<double> log1;
solvers.getResidualLog(log1);
for(auto x : log1){
    cout << x << endl;
}
```

ICCG以外のソルバの実行例です。
```
SRLfem::MatSolvers solvers;

// [IC-MRTR]
//@param1: 行数
//@param2: 収束判定
//@param3: 最大反復数
//@param4: 加速係数（1以上だとその値で固定。負の値にすると、その絶対値を超えない範囲で自動計算）
//@param5: 疎行列A
//@param6: 右辺ベクトルb
//@param7: 解ベクトルx
//@param8: ソルバ開始前の解ベクトルゼロ初期化有無。trueでxをゼロで初期化。falseだと何もしない(省略時はfalse)。
bool bl = solvers.solveICMRTR(size, 1.0e-8, 10000, 1.05, matAs, vec, results, true);

// [SGS-MRTR]: 加速係数不要。内部で対角スケーリング必須なので、フラグがOFFでも対角スケーリングは実行されます
//@param1: 行数
//@param2: 収束判定
//@param3: 最大反復数
//@param5: 疎行列A
//@param6: 右辺ベクトルb
//@param7: 解ベクトルx
//@param8: ソルバ開始前の解ベクトルゼロ初期化有無。trueでxをゼロで初期化。falseだと何もしない(省略時はfalse)。
bool bl = solvers.solveSGSMRTR(size, 1.0e-8, 10000, matAs, vec, results, true);

```

</div>