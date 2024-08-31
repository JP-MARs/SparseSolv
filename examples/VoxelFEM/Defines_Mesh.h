/*
★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆

　メッシュ分割数などの定数の定義

★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆★☆
*/

#ifndef DEFINE_MESH_VOX
#define DEFINE_MESH_VOX

/* メッシュを切る範囲X（０から～） */
constexpr double X_WIDTH = 0.1;
/* メッシュを切る範囲Y（０から～） */
constexpr double Y_WIDTH = 0.1;
/* メッシュを切る範囲Z（０から～） */
constexpr double Z_WIDTH = 0.1;

/* メッシュのx長さ（細かい方） */
constexpr double X_L1 = 0.001;
/* メッシュのy長さ（細かい方） */
constexpr double Y_L1 = 0.001;
/* メッシュのz長さ（細かい方） */
constexpr double Z_L1 = 0.001;

/* メッシュのx長さ（粗い方） */
constexpr double X_L2 = 0.01;
/* メッシュのy長さ（粗い方） */
constexpr double Y_L2 = 0.01;
/* メッシュのz長さ（粗い方） */
constexpr double Z_L2 = 0.01;

/* メッシュの粗さ変更点のx座標 */
constexpr double CHANGE_X = 0.01;
/* メッシュのy長さ（粗い方） */
constexpr double CHANGE_Y = 0.01;
/* メッシュのz長さ（粗い方） */
constexpr double CHANGE_Z = 0.01;

/* 細かいメッシュの数X */
constexpr int MESH_X1 = ((int)((CHANGE_X / X_L1)));
/* 細かいメッシュの数Y */
constexpr int MESH_Y1 = ((int)((CHANGE_Y / Y_L1)));
/* 細かいメッシュの数Z */
constexpr int MESH_Z1 = ((int)((CHANGE_Z / Z_L1)));

/* 粗いメッシュの数X */
constexpr int MESH_X2 = ((int)( (X_WIDTH-CHANGE_X) / X_L2 ));
/* 粗いメッシュの数Y */
constexpr int MESH_Y2 = ((int)( (Y_WIDTH-CHANGE_Y) / Y_L2 ));
/* 粗いメッシュの数Z */
constexpr int MESH_Z2 = ((int)( (Z_WIDTH-CHANGE_Z) / Z_L2 ));

/* メッシュの数合計X */
constexpr int TOTAL_MESH_X = (MESH_X1+MESH_X2);
/* メッシュの数合計Y */
constexpr int TOTAL_MESH_Y = (MESH_Y1+MESH_Y2);
/* メッシュの数合計Z */
constexpr int TOTAL_MESH_Z = (MESH_Z1+MESH_Z2);

#endif
