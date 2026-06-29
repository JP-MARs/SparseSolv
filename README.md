# SparseSolv
SparseMatrix and Matrix Solvers including 
- shifted-ICCG
- shifted-IC+MRTR
- Eisenstat's Symmetric Gauss-Seidel-MRTR.

### Contributors
 - Takahiro Sato (Muroran institute of technology, JAPAN)
 - Shingo Hiruma (Hokkaido University, JAPAN)
 - Kengo Sugahara (Kindai University, JAPAN)
 - Hideaki Nagamine (Gifu University)
 - Tomonori Tsuburaya (Fukuoka University, JAPAN)
 

# 本OSSについて
日本の磁界解析分野で独自に発展してきた「加速係数付きICCG法」を含む、C++の疎行列ライブラリです。Eigenの疎行列クラスのラッパーとしての疎行列クラスと、疎行列クラスに対する線形ソルバで構成されています。疎行列クラスに対する基本的な四則演算や行列の変形操作、ソルバ用の前処理操作も複数用意しています。

電磁界の不定な方程式を解く際に有効な加速係数付きICCG法をはじめとして、日本の磁界系研究者によって検討されてきた特有のソルバを複数実装しています。電磁界のFEMにおいて、ゲージ固定などで不定性を取り除くことなく、ICCG法で不定な方程式をそのまま解くことが可能です。実装されているソルバは下記です。
- 加速係数付きICCG法(実数対称/複素対称)
- 加速係数付きMRTR法(実数対称/複素対称)
- Eisenstatの方法を用いたMRTR法(実数対称/複素対称)
- Eigenのデフォルトソルバのラッパー

なお、Eigenのデフォルト並列機能を用いた並列化は可能ですが、大規模行列用ソルバとしての高度な並列化(MPI対応や多並列時のスケーラビリティ)には対応していません。

## Pythonラッパー
本ソルバはPybind11を利用したPythonバインディングも提供しています。Pythonで利用したい場合はそちらを利用可能です。

## 利用例
examplesフォルダに、本ライブラリの利用例をアップしています。

### Pybind_example.py
Pythonでの基本的な利用法のスクリプトです。適当な疎行列を作成し、ソルバで適当な線形方程式を解くサンプルです。

### VoxelFEM
3次元静磁場FEMを加速係数付きICCG法で解く簡単なサンプルです。
ボクセル状の単純なメッシュを読み込み、静磁場のFEM方程式を構築します。その際、ゲージ固定で不定性を除去せず、元の不定な方程式のまま有限要素方程式を構築します。そして、その方程式を加速係数付きICCG法でそのまま解きます。（ゲージ固定をするより、不定な方程式をそのまま解いたほうが高速に解けることが証明されています。）

ボクセルモデルのメッシュは下記にあります（独自フォーマットですが）

https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_850000ele.zip<br>
https://u.muroran-it.ac.jp/it-elec-lab/open_data/voxel_data/MeshData_2000000ele.zip<br>

zipファイルを解凍した後、`*.vxldata`の拡張子のファイル群を`examples/VoxelFEM/MeshData/`においてください。

## そのほか
詳しい内容や、学術的な背景の説明は、解説用ホームページをご覧ください（本ページ右上にリンクがあります）。

