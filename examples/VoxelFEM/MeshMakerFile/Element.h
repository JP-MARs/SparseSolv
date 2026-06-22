/*
//=======================================================
// ■ Element
//=======================================================
// FEMの要素定義クラス
//=======================================================*/
class Element {
private:
	int ID;																	/* この要素のID */
	int x;																	/* 原点から何番目の要素か(ｘ方向) */
	int y;																	/* 原点から何番目の要素か(ｙ方向) */
	int z;																	/* 原点から何番目の要素か(ｚ方向) */
	int *nodes;																/* この要素の節点たち */
	int *edges;																/* この要素の辺たち */
	bool isT;																/* この要素の辺でＴの値持ちがあるか判別 */
	static Node_Mesh *NodeData;													/* 節点データ保持配列 */
	static Edge *EdgeData;													/* 辺データ保持配列 */
public:
	Element();																/* コンストラクタ */
	~Element();
	void set(int xx, int yy, int zz, int *n, int *e, int id);				/* セット */
	int getID(){return ID;};
	int *getNodes(){return nodes;};
	int *getEdges(){return edges;};
	int getNodeID(int i){return nodes[i];};
	int getEdgeID(int i){return edges[i];};
	int getX(){return x;};
	int getY(){return y;};
	int getZ(){return z;};
	void setIsT(bool bl){isT = bl;};
	bool is_T(){return isT;};
	static void setNodeData(Node_Mesh *n){NodeData = n;};
	static void setEdgeData(Edge *n){EdgeData = n;};
	void getGrav(double *x);													/* 要素の中心位置返し */
	double getMinX(){return  NodeData[nodes[0]].getPosX();};					/* 要素の最小X位置返し */
	double getMinY(){return  NodeData[nodes[0]].getPosY();};					/* 要素の最小Y位置返し */
	double getMinZ(){return  NodeData[nodes[0]].getPosZ();};					/* 要素の最小Z位置返し */
	double getMaxX(){return  NodeData[nodes[6]].getPosX();};					/* 要素の最大X位置返し */
	double getMaxY(){return  NodeData[nodes[6]].getPosY();};					/* 要素の最大Y位置返し */
	double getMaxZ(){return  NodeData[nodes[6]].getPosZ();};					/* 要素の最大Z位置返し */
};
Node_Mesh *Element::NodeData=NULL;
Edge *Element::EdgeData=NULL;
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Element::Element(){
	nodes = new int[8];
	edges = new int[12];
}
Element::~Element(){
	delete[] nodes;
	delete[] edges;
}
/*//=======================================================
  // ● データのセッタ
  //=======================================================*/
void Element::set(int xx, int yy, int zz, int *n, int *e, int id){
	ID = id;
	x = xx;
	y = yy; 
	z = zz;
	for(int i = 0 ; i < 8 ; i++){
		nodes[i] = n[i];
	}
	for(int i = 0 ; i < 12 ; i++){
		edges[i] = e[i];
	}
}
/*//=======================================================
  // ● 要素の中心位置返し
  //=======================================================*/
void Element::getGrav(double *x){
	double g_x=0, g_y=0, g_z=0;
	for(int i = 0 ; i < 8 ; i++){
		g_x += NodeData[nodes[i]].getPosX();
		g_y += NodeData[nodes[i]].getPosY();
		g_z += NodeData[nodes[i]].getPosZ();
	}
	g_x /= 8.0;
	g_y /= 8.0;
	g_z /= 8.0;
	x[0] = g_x; x[1] = g_y; x[2] = g_z;
}