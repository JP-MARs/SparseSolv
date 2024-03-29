﻿
#include "Node.h"

/*
//=======================================================
// ■ Node
//=======================================================
// FEMの節点定義クラス
//=======================================================*/

/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Node::Node(int I, double xx, double yy, double zz){
	set(I, xx, yy, zz);
}
/*//=======================================================
  // ● 節点座標セッタ
  //=======================================================*/
void Node::set(int I, double xx, double yy, double zz){
	ID = I;
	x = xx; y = yy; z = zz;
}