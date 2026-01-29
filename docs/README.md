# 行列解法・ICCG法関連文献

このフォルダには、有限要素法(FEM)などの電磁界解析で用いられる疎行列解法、主にICCG法やCOCG法に関する文献が格納されています。
特に、**日本国内の電磁界解析コミュニティで独自に発展・定着した実装ノウハウ**（Shifted ICCGや $LDL^{T}$ 分解の採用など）や、辺要素FEM特有の収束性に関する議論に焦点を当てています。

## 理論的背景 (Theoretical Background)

電磁界解析（特にFEM）において、離散化により得られる連立一次方程式 $A x=b$ は、大規模かつ疎（Sparse）な構造を持ちます。
この行列 $A$ の性質により、適用すべきクリロフ部分空間法が異なります。

### 1. CG法とICCG法
**適用対象**: 実対称正定値行列 (Real Symmetric Positive Definite)

関数 $\phi(x)$ の最小化問題と等価です。

```math
\phi(x) = \frac{1}{2} x^T A x - x^T b
```

探索方向 $p_k$ が $A$-直交 (共役) であるように選ばれます。

```math
p_i^T A p_j = 0 \quad (i \neq j)
```

収束性を改善するために、不完全コレスキー分解 (Incomplete Cholesky Decomposition) を前処理行列 $M$ として用いるのが **ICCG法** です。

```math
M = LDL^{T} \approx A
```

### 2. COCG法 (Conjugate Orthogonal CG)
**適用対象**: 複素対称行列 (Complex Symmetric)

```math
A = A^T \quad (\text{ただし通常 } A \neq A^H)
```
渦電流解析などの周波数領域FEMでは、行列は複素対称になります。これに対して、BiCG法よりも効率的なCOCG法 (van der Vorst, 1990) が用いられます。

#### 直交性の定義の違い
通常の複素内積（エルミート内積）ではなく、双線形形式 (bilinear form) を用います。

```math
(x, y)_* = x^T y \quad (\text{共役をとらない})
```
これにより、CG法と同様の3項漸化式で更新が可能となり、メモリ効率と計算速度が向上します。

---

## 電磁界解析における行列解法の歴史的変遷 (Historical Evolution)

日本の計算電磁気学（Computational Electromagnetics, CEM）分野において、現在の標準的な解法（Shifted ICCG）が確立されるまでには、1980年代から90年代にかけてのハードウェアの進化と、定式化のパラダイムシフトに伴う深い議論の歴史がありました。

### 1. バンド行列とベクトル計算機の時代 (~1980年代中盤)
初期の3次元FEM解析への挑戦は、メモリリソースと計算速度との戦いでした。
*   当初は「バンド行列ソルバ」や「スカイライン法」などの直接法が用いられていましたが、未知数が数1,000を超えると当時のメインフレームでも計算不能に陥りました。
*   **Meijerink & van der Vorst [6]** (1977) による不完全コレスキー分解を用いた前処理付きCG法（ICCG法）の提案は、この状況を一変させました。日本国内では、特に日立・富士通・NECなどの国産スーパーコンピュータ（ベクトル計算機）上でのICCG法の高速化（超平面法によるベクトル化など）が、企業・大学で熱心に研究されました。
*   **岡山大学・中田高義研究室の貢献**: 中田高義教授率いる岡山大学の研究グループ（藤原耕二氏を含む）は、ICCG法の収束特性改善に関して日本をリードする研究を行いました **[3]**。特に：
    *   リナンバリング（reverse Cuthill-McKee法など）による行列帯幅縮小と収束改善
    *   **シフト係数（加速係数）** の導入による収束安定化
    *   磁性材料の標準磁気特性測定法の確立
    など、実用化に不可欠な基盤技術を体系的に研究し、「良い前処理を行えば反復法は速い」という認識の定着に大きく貢献しました。中田研出身の研究者は現在も国内外で電磁界解析の第一線で活躍しています。

### 2. 辺要素の登場と「特異性の壁」 (1980年代後半〜)
T. Weiland (FIT) や A. Bossavit (Edge Elements) らにより、スプリアス解を含まないベクトル解析として「辺要素」が理論化されました。しかし、これを実装しようとした研究者たちは、大きな壁に直面しました。
*   **特異行列 (Singular Matrix)**: $\text{rot grad} \phi = 0$ の関係により、係数行列 $A$ が不定（Singular）になり、通常のコレスキー分解やICCG法がそのままでは適用できなかったのです（ピボットがゼロになる）。
*   **世界的な迷走**:
    *   欧米の研究者の多くは、「数学的に正しい行列にする」ことを重視し、Tree-Cotree分解などのグラフ理論を用いてゲージ不定性を除去（Gauge Fixing）する方法を模索しました。しかし、これは前処理の実装を極めて複雑にし、スパース性を損なうものでした。
    *   あるいは、ペナルティ項を加えて無理やり正定値化する方法も試みられましたが、ペナルティ係数の調整という新たな問題を招きました。

### 3. 日本独自のブレイクスルー：Shifted ICCG (1990年代初頭)
この停滞を打破したのが、日本国内の現場からの実用的なアプローチでした。
*   **亀有昭久氏（元・三菱原子力工業、現SSIL）による辺要素FEMへの適用 [2]**:
    *   Kershaw **[7]** が1978年にプラズマシミュレーションの負ピボット問題を回避するために提案し、Manteuffel **[5]** が理論を体系化した対角シフト付き不完全分解（Shifted IC）を、辺要素FEMの特異行列問題に適用し、その有効性を実証しました。
    *   数学的な厳密さよりも「解けること」を優先した実用的アプローチにより、ゲージ固定のような複雑な前処理を必要とせず、 $LDL^{T}$ 分解（ルートフリー / square-root-free）を用いることで計算も高速な **Shifted ICCG** 法が確立されました。
*   **藤原耕二氏（岡山大、現同志社大）らによるシフト係数の最適化 [3]**:
    *   シフト係数（加速係数）の最適値を求める方法が研究され、実対称行列に対しては自動決定法が提案されています **[8]**。
*   **Van der Vorst との交流**:
    *   オランダの H. A. van der Vorst 教授は、電磁界解析コミュニティ（COMPUMAG等）とも深い関わりを持ち、1990年に複素対称行列向けの **COCG法 [4]** を提案しました。
    *   日本の研究者は即座にこれを取り入れ、複素数に対しても Shifted IC + COCG を適用することで、渦電流解析などの周波数応答解析においても「特異行列のまま解く」スタイルを確立しました。

### 4. 理論的裏付けとデファクトスタンダード化 (2000年代〜)
*   長らく「経験的に動く」とされてきたこの手法に対し、**五十嵐 (北大) [1]** は数学的なメスを入れました。
    *   「ゲージ固定を行うと、むしろ行列のスペクトル特性（固有値分布）が悪化し、収束が遅くなること」を理論的に示しました。
    *   逆に、「特異行列のままであっても、適切な初期値（勾配成分ゼロ）とKrylov部分空間法を組み合わせれば、実質的に良条件の問題として解ける」ことを証明しました。
*   これにより、Shifted ICCGは単なる「便宜的な手法」から「理論的裏付けのある最適解」へと昇華され、日本発の技術として、現在も多くの商用・アカデミックコードのコアエンジンとして稼働しています。

---

## 実用的実装の検討 (Practical Implementation)

**岡本ら [9]** (2014) は、国内10研究機関でICCG法の実装をベンチマーキングし、高速化のための実用的な実装法を体系化しました。

### ベンチマークモデル
*   **ボックスシールドモデル**: 節点数558,256、要素数539,105、未知数（DoF）1,594,320、係数行列の非零要素数26,618,443
*   六面体一次辺要素、磁気ベクトルポテンシャル定式化（線形静磁界）
*   CRS（Compressed Row Storage）形式で疎行列を格納

### 対角スケーリングの数学的等価性
対角スケーリング（点ヤコビ前処理）の有無がIC前処理後の収束特性に影響しないことが数学的に証明されました。

スケーリング行列 $S$ を対角成分 $\sqrt{a_{ii}}$ で定義すると：

```math
\bar{A} = S^{-1}AS^{-T}, \quad \bar{x} = S^{T}x, \quad \bar{b} = S^{-1}b
```

IC分解 $A \approx CC^{T}$, $\bar{A} \approx \bar{C}\bar{C}^{T}$ に対して、 $\bar{C} = S^{-1}C$ の関係より：

```math
(C^{-1}AC^{-T})(C^{T}x) = C^{-1}b
```

となり、対角スケーリングの有無に関係なく前処理後の方程式は同値となります。

### 実装上の最適化ポイント
1.  **IC分解の形式**: $LDL^{T}$ 分解（ルートフリー）において、 $L$ の対角成分を1とする形式（IC type 2）を採用すると、前進代入で40%以上の高速化が可能
2.  **前進・後退代入**: 式 $(LD)D^{-1}(LD)^{T}u = r$ の形式に変形し、 $LD$ を格納することで対角成分のメモリアクセスを削減
3.  **後退代入のループ方向**: 逆方向にループすることで連続メモリアクセスによる高速化
4.  **行列ベクトル積**: 下三角部分のループ内で上三角部分の寄与も同時に計算し、メモリアクセス回数を半減
5.  **キャッシュミス回避**: 間接参照によるキャッシュミスを最小化する実装が重要

### シフト係数（加速係数）
*   シフト係数 $\gamma$ は問題依存であり、最適値は行列の性質によって異なります
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

---
*Generated by geminiさん on 2026-01-29*
