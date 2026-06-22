<div style="text-align: justify;">

# Technical background
## [はじめに]　本OSS（加速係数ICCG）の技術背景
&ensp;本ライブラリで提供している加速係数付きICCG法が発達した技術的背景を紹介します。
特に、**日本国内の電磁界解析コミュニティで独自に発展・定着した実装ノウハウ**（Shifted ICCGや $LDL^{T}$ 分解の採用など）や、辺要素FEM特有の収束性に関する議論に焦点を当てています。


## １．理論的背景 (Theoretical Background)

&ensp;電磁界解析（特にFEM）において、離散化により得られる連立一次方程式 $A x=b$ は、大規模かつ疎（Sparse）な構造を持ちます。
この行列 $A$ の性質により、適用すべきクリロフ部分空間法が異なります。

### 1.A. CG法とICCG法
**適用対象**: 実対称正定値行列 (Real Symmetric Positive Definite)

関数 $\phi(x)$ の最小化問題と等価です。

$$
\phi(x) = \frac{1}{2} x^T A x - x^T b
$$

探索方向 $p_k$ が $A$-直交 (共役) であるように選ばれます。

$$
p_i^T A p_j = 0 \quad (i \neq j)
$$

&ensp;収束性を改善するために、不完全コレスキー分解 (Incomplete Cholesky Decomposition) を前処理行列 $M$ として用いるのが **ICCG法** です。

$$
M = LDL^{T} \approx A
$$

### 1.B. COCG法 (Conjugate Orthogonal CG)
**適用対象**: 複素対称行列 (Complex Symmetric)

$$
A = A^T \quad (\text{ただし通常 } A \neq A^H)
$$
渦電流解析などの周波数領域FEMでは、行列は複素対称になります。これに対して、BiCG法よりも効率的なCOCG法 (van der Vorst, 1990) が用いられます。

#### 直交性の定義の違い
通常の複素内積（エルミート内積）ではなく、双線形形式 (bilinear form) を用います。

$$
(x, y)_* = x^T y \quad (\text{共役をとらない})
$$

これにより、CG法と同様の3項漸化式で更新が可能となり、メモリ効率と計算速度が向上します。

### 1.C. 加速係数付きIC分解(Incomplete Cholesky factorization with shitfed factor)
&ensp;上述の通り、CG法を使う場合は、通常、行列の前処理を行うことで収束性が改善されます。
特に「不完全コレスキー分解」を前処理として利用するCG法をICCG法と呼びます。一般的な行列 $A$ ではICCG法が有効ですが、ここで電磁気のFEMに特有の現象が影響します。<br>
&ensp;電磁気FEMの場合、 $\text{rot grad} \phi = 0$ の関係により、係数行列 $A$ が不定（Singular）になり、通常のコレスキー分解やICCG法がそのままでは適用できなくなります。（ピボットがゼロになる）
不定性をなくすため、木構造によるゲージを課し不定性を無くすことができるのですが、極端に収束が遅くなることが知られていました。そこで日本で考案され、日本では広く利用されているのが「加速係数付きのICCG法」です。<br>
&ensp;加速係数付きのIC分解では、前処理行列を作成する際、元の行列 $A$ そのものではなく、対角成分を $\alpha$ 倍した行列 $A'$ に対して不完全コレスキー分解を実行します。
つまり、元の行列 $A$ を、その対角成分 $D$ とそれ以外 $K$ に分けて以下のように表記します。

$$
A = D+K
$$

そして、IC分解を、 $A$ ではなく、対角成分を $\alpha$ 倍した行列 $A'$ に対して適用します。

$$
A' = \alpha D+K
$$

この $\alpha$ は加速係数と日本では呼ばれています。このようにすれば、不定な行列であってもIC分解を適用することができ、ICCG法を適用することができます。そして、この方法を用いたほうが、ゲージを課し不定性を無くすよりも非常に高速に解けることが証明されています(後述)。そのため、日本の電磁界解析では、加速係数付きICCG法が広く利用されています。<br>
&ensp;加速係数の値は問題依存であり、適切な値を設定することで収束性が改善されます。ただし固定値でも通常のICCG法より十分に高速化されます。固定値を使う場合は、1.05や1.2を使うことが多いようです。

> **重要：**
> 行列 $A$ の修正は「前処理行列 $M$」を作るためだけに使用します。それ以外、例えば、残差 $r = b - Ax$ の計算には、**必ず元の行列 $A$** を使用します。（そうしないと、解くべき問題自体が変わってしまう。）

---

## 2. 電磁界解析における行列解法の歴史的変遷 (Historical Evolution)

&ensp;日本の計算電磁気学（Computational Electromagnetics）分野において、現在の標準的な解法（Shifted ICCG）が確立されるまでには、1980年代から90年代にかけてのハードウェアの進化と、定式化に関する深い議論の歴史がありました。

### 2.A. バンド行列とベクトル計算機の時代 (~1980年代中盤)
&ensp;初期の3次元FEM解析への挑戦は、メモリリソースと計算速度の問題でした。
*   当初は「バンド行列ソルバ」や「スカイライン法」などの直接法が用いられていましたが、未知数が数1,000を超えると当時のメインフレームでも計算不能になりました。
*   **Meijerink & van der Vorst [6]** (1977) により、不完全コレスキー分解を用いた前処理付きCG法（ICCG法）が提案されました。日本国内では、特に日立・富士通・NECなどの国産スーパーコンピュータ（ベクトル計算機）上でのICCG法の高速化（超平面法によるベクトル化など）が、企業・大学で熱心に研究されました。
*   **岡山大学・中田高義研究室の貢献**: 中田高義教授率いる岡山大学の研究グループ（藤原耕二氏を含む）は、ICCG法の収束特性改善に関して日本をリードする研究を行いました **[3]**。特に：
    *   リナンバリング（reverse Cuthill-McKee法など）による行列帯幅縮小と収束改善
    *   **シフト係数（加速係数）** の導入による収束安定化
    *   磁性材料の標準磁気特性測定法の確立
    
など、実用化に不可欠な基盤技術を体系的に研究し、「良い前処理を行えば反復法は速い」という認識の定着に大きく貢献しました。中田研出身の研究者は現在も国内外で電磁界解析の第一線で活躍しています。

### 2.B. 辺要素の登場と「特異性の壁」 (1980年代後半〜)
T. Weiland (FIT) や A. Bossavit (Edge Elements) らにより、スプリアス解を含まないベクトル解析として「辺要素」が理論化されました。しかし、実装段階で大きな問題が生じました。
*   **特異行列 (Singular Matrix)**: $\text{rot grad} \phi = 0$ の関係により、係数行列 $A$ が不定（Singular）になり、通常のコレスキー分解やICCG法はそのまま適用できないという問題が生じました。（ピボットがゼロになる）。
*   **当時の世界（日本以外）の状況**:
    *   欧米では、「数学的に正しい行列にする」ことを重視し、Tree-Cotree分解などのグラフ理論を用いてゲージ不定性を除去（Gauge Fixing）する方法が検討されました。しかし、これは前処理の実装を複雑にし、スパース性を損なう可能性がありました。
    *   ほかに、ペナルティ項を加えて無理やり正定値化する方法も試みられました。しかし、ペナルティ係数の調整という新たな問題が生じました。

### 2.C. 日本独自のブレイクスルー：Shifted ICCG (1990年代初頭)
一方、日本では、企業の現場からの実用的なアプローチが生まれました。
*   **亀有昭久氏（元・三菱原子力工業、現SSIL）による辺要素FEMへの適用 [2]**:
    *   Kershaw **[7]** が1978年にプラズマシミュレーションの負ピボット問題を回避するために提案し、Manteuffel **[5]** が理論を体系化した対角シフト付き不完全分解（Shifted IC）を、辺要素FEMの特異行列問題に適用し、その有効性を実証しました（前述の加速係数付きIC分解）。
    *   数学的な厳密さよりも「解けること」を優先した実用的アプローチにより、ゲージ固定のような複雑な前処理を必要とせず、 $LDL^{T}$ 分解（ルートフリー / square-root-free）を用いることで計算も高速な **加速係数付きICCG法** 法が確立されました。
*   **藤原耕二氏（岡山大、現同志社大）らによるシフト係数の最適化 [3]**:
    *   加速係数の最適値を求める方法が研究され、実対称行列に対しては自動決定法が提案されています **[8]**。
*   **複素対称行列への拡張（ICCOCG）**:
    *   日本国内では、Van der Vorst の1990年の論文 **[4]** が発表される以前から、複素対称行列に対する ICCG法の拡張（ICCOCG）が実用化されていました。渦電流解析などの周波数応答解析において、複素数に対しても 加速係数付きICCG法 を適用し「特異行列のまま解く」スタイルが確立されていました。
    *   オランダの H. A. van der Vorst 教授は、電磁界解析コミュニティ（COMPUMAG等）とも深い関わりを持ち、1990年に複素対称行列向けの **COCG法 [4]** を発表しました。この論文により、双線形形式に基づく直交性の定義と3項漸化式による効率的な更新という**理論的裏付け**が明確になり、日本で先行して実用化されていた手法の数学的正当性が確立されました。

### 2.D. 理論的裏付けとデファクトスタンダード化 (2000年代〜)
*   長らく「経験的に動く」とされてきた「加速係数付きICCG法」ですが、**五十嵐一教授 (北海道大学) [1]** により、数学的な理論が研究されました。
    *   「ゲージ固定を行うと、むしろ行列のスペクトル特性（固有値分布）が悪化し、収束が遅くなること」を理論的に示しました。
    *   逆に、「特異行列のままであっても、適切な初期値（勾配成分ゼロ）とKrylov部分空間法を組み合わせれば、実質的に良条件の問題として解ける」ことを証明しました。
*   これにより、加速係数付きICCG法は「理論的裏付けのある最適解」となり、日本発の技術として、現在も多くの商用・アカデミックコードのコアエンジンとして利用されています。

---

## 実用的実装の検討 (Practical Implementation)

**岡本教授（現 法政大学）ら [9]** (2014) は、国内10研究機関でICCG法の実装をベンチマーキングし、高速化のための実用的な実装法をまとめました。

### ベンチマークモデル
*   **ボックスシールドモデル**: 節点数558,256、要素数539,105、未知数（DoF）1,594,320、係数行列の非零要素数26,618,443
*   六面体一次辺要素、磁気ベクトルポテンシャル定式化（線形静磁界）
*   CRS（Compressed Row Storage）形式で疎行列を格納

### 対角スケーリングの数学的等価性
対角スケーリング（点ヤコビ前処理）の有無がIC前処理後の収束特性に影響しないことが数学的に証明されました。

スケーリング行列 $S$ を対角成分 $\sqrt{a_{ii}}$ で定義すると：

$$
\bar{A} = S^{-1}AS^{-T}, \quad \bar{x} = S^{T}x, \quad \bar{b} = S^{-1}b
$$

IC分解 $A \approx CC^{T}$, $\bar{A} \approx \bar{C}\bar{C}^{T}$ に対して、 $\bar{C} = S^{-1}C$ の関係より：

$$
(C^{-1}AC^{-T})(C^{T}x) = C^{-1}b
$$

となり、対角スケーリングの有無に関係なく前処理後の方程式は同値となります。

### 実装上の最適化ポイント
1.  **IC分解の形式**: $LDL^{T}$ 分解（ルートフリー）において、 $L$ の対角成分を1とする形式（IC type 2）を採用すると、前進代入で40%以上の高速化が可能
2.  **前進・後退代入**: 式 $(LD)D^{-1}(LD)^{T}u = r$ の形式に変形し、 $LD$ を格納することで対角成分のメモリアクセスを削減
3.  **後退代入のループ方向**: 逆方向にループすることで連続メモリアクセスによる高速化
4.  **行列ベクトル積**: 下三角部分のループ内で上三角部分の寄与も同時に計算し、メモリアクセス回数を半減
5.  **キャッシュミス回避**: 間接参照によるキャッシュミスを最小化する実装が重要

### 加速係数
*   加速係数 $\gamma$ は問題依存であり、最適値は行列の性質によって異なります
*   岡本らのベンチマーキング **[9]** では比較のため $\gamma = 1.06$ を統一して使用
*   **自動決定法**: 北尾・高橋・藤原ら **[8]** は、残差と汎関数に基づくシフト係数の自動決定法を提案しており、実対称行列に対して適用可能です

---

## 参考文献 (References)

1.  **H. Igarashi, T. Honma**, "On the property of the ICCG method applied to the finite element equation for quasi-static fields," *IEEE Transactions on Magnetics*, vol. 38, no. 2, pp. 565-568, 2002.
2.  **A. Kameari** (SSIL), "Calculation of transient 3D eddy current using edge-elements," *IEEE Transactions on Magnetics*, vol. 26, no. 2, pp. 466-469, 1990.
3.  **K. Fujiwara, T. Nakata**, "Acceleration of convergence characteristic of the ICCG method," *IEEE Transactions on Magnetics*, vol. 29, no. 2, pp. 1958-1961, 1993.
4.  **H. A. van der Vorst** and J. B. M. Melissen, "A Petrov-Galerkin type method for solving Ax=b, where A is symmetric complex," *IEEE Transactions on Magnetics*, vol. 26, no. 1, pp. 706-708, 1990.
5.  **T. A. Manteuffel**, "An incomplete factorization technique for positive definite linear systems," *Mathematics of Computation*, vol. 34, no. 150, pp. 473-497, 1980. [DOI: 10.1090/S0025-5718-1980-0559197-0](https://doi.org/10.1090/S0025-5718-1980-0559197-0)
6.  **J. A. Meijerink, H. A. van der Vorst**, "An iterative solution method for linear systems of which the coefficient matrix is a symmetric M-matrix," *Mathematics of Computation*, vol. 31, no. 137, pp. 148-162, 1977. [DOI: 10.1090/S0025-5718-1977-0438681-4](https://doi.org/10.1090/S0025-5718-1977-0438681-4)
7.  **D. S. Kershaw**, "The incomplete Cholesky-conjugate gradient method for the iterative solution of systems of linear equations," *Journal of Computational Physics*, vol. 26, no. 1, pp. 43-65, 1978. [DOI: 10.1016/0021-9991(78)90098-0](https://doi.org/10.1016/0021-9991(78)90098-0)
8.  **J. Kitao, Y. Takahashi, K. Fujiwara, T. Mifune, T. Iwashita**, "Automatic determination of acceleration factor in shifted ICCG method for 3-D electromagnetic field analyses," *IEEE Transactions on Magnetics*, vol. 49, no. 5, pp. 1741-1744, 2013.
9.  **Y. Okamoto, Y. Takahashi, K. Fujiwara, A. Ahagon, T. Mifune, T. Iwashita**, "辺要素有限要素解析から得られる実対称線形方程式求解におけるICCG法の実用的実装の検討," *電気学会論文誌B*, vol. 134, no. 9, pp. 767-776, 2014. [DOI: 10.1541/ieejpes.134.767](https://doi.org/10.1541/ieejpes.134.767)

</div>