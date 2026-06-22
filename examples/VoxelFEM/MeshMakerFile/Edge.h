/*
//=======================================================
// ■ Edge
//=======================================================
// FEMの辺定義クラス
//=======================================================*/
class Edge {
private:
	int ID;												/* この辺のID */
	int stat_node_id;									/* 辺の始点NodeのID */
	int end_node_id;									/* 辺の終点NodeのID */
	bool zero_edge;										/* 既知ゼロの辺か */
	bool T_edge;										/* 電流ベクトルポテンシャルの辺か */
	double *T_potential;								/* 電流ベクトルポテンシャルの値(ポインタ参照・実態は上位クラスが保持) */
public:
	Edge();												/* コンストラクタ */
	Edge(int sID, int eID);								/* コンストラクタ2 */
	void set(int sID, int eID);							/* 座標セット */
	int getSID(){return stat_node_id;};
	int getEID(){return end_node_id;};
	int getID(){return ID;};
	void setID(int id){ID=id;};
	void setZeroEdge(bool bl){zero_edge=bl;};
	bool isZeroEdge(){return zero_edge;};
	void setT_Edge(bool bl){T_edge=bl;};
	bool isT_Edge(){return T_edge;};
	double get_Tvalue(){return *T_potential;};
	void set_Tvalue(double *x){T_potential = x;};
};
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Edge::Edge(){
	zero_edge = false;
	T_edge = false;
}
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Edge::Edge(int sID, int eID){
	zero_edge = false;
	T_edge = false;
	set(sID, eID);
}
/*//=======================================================
  // ● 始点・終点セッタ
  //=======================================================*/
void Edge::set(int sID, int eID){
	stat_node_id = sID; end_node_id = eID;
}