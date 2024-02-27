/*
★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆

　定数などの定義

★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆
*/
#ifndef DEFINE_SET_VOX
#define DEFINE_SET_VOX

/* コイルの設定 */
constexpr double COIL_MIN = (0.03);
constexpr double COIL_MAX = (0.04);
constexpr double COIL_Z = (0.02);

constexpr double IRON_Z = (0.04);

/* 磁性体の透磁率（線形計算用）*/
/* set relativepermeability in iron region */
constexpr double PI = (3.14159265358979);

/*真空透磁率*/
constexpr double MYU0 = (4.0*PI*1.0e-7);

constexpr double Air_perm = MYU0;
constexpr double Iron_perm = 1000.0*MYU0;

constexpr int MAT_AIR  = 1;
constexpr int MAT_IRON  = 2;
constexpr int MAT_COIL  = 3;

/*電流[A]*/
constexpr double CURRENT = (1.0);

/*コイルターン数[turns]*/
constexpr double NUMTURNS = (1000.0);

///*電流ターン[A*turns]*/
constexpr double TOTAL_CURRENT = (CURRENT*NUMTURNS);

/*連立方程式の収束判定*/
constexpr double EPS_SOLVER  = (1.0e-6);

/*反復回数の上限*/
constexpr int MAX_ITR_ICCG = 10000;

/*IC分解の加速係数*/
constexpr double GA_IC = 1.1;

#endif
