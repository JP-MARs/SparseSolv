/*
//=======================================================
// ■ SNode
//=======================================================
// サーチ木の節点を定義するクラス
//=======================================================*/
class SNode;
class SNode {
private:
	int id;												/* このノードの値 */
	int x;
	int y;
	int z;
	bool active;
	SNode *right;										/* 右の子 */
	SNode *center;										/* 中心の子 */
	SNode *left;										/* 左の子 */
	int edge_1;											/* x側の辺番号 */
	int edge_2;											/* y側の辺番号 */
	int edge_3;											/* z側の辺番号 */
public:
	SNode();											/* コンストラクタ */
	~SNode(){;};
	int getID(){return id;};							/* 節の内容ゲッタ */
	void setID(int x){id = x;};							/* 節の内容セッタ */
	SNode *getRight(){return right;};					/* 右の子のポインタ返し */
	SNode *getCenter(){return center;};					/* 中心の子のポインタ返し */
	SNode *getLeft(){return left;};						/* 左の子のポインタ返し */
	SNode **getRightPtr(){return &right;};				/* 右の子のポインタポインタ返し */
	SNode **getCenterPtr(){return &center;};			/* 中心の子のポインタポインタ返し */
	SNode **getLeftPtr(){return &left;};				/* 左の子のポインタポインタ返し */
	void setRight(SNode *n){right = n;};				/* 右の子のポインタセット */
	void setCenter(SNode *n){center = n;};				/* 中心の子のポインタセット */
	void setLeft(SNode *n){left = n;};					/* 左の子のポインタセット */
	void setX(int xx){x = xx;};
	void setY(int yy){y = yy;};
	void setZ(int zz){z = zz;};
	int getX(){return x;};
	int getY(){return y;};
	int getZ(){return z;};
	bool isActive(){return active;};
	void Activate(){active=true;};
	void setEdge1(int x){edge_1=x;};
	void setEdge2(int x){edge_2=x;};
	void setEdge3(int x){edge_3=x;};
	int getEdge1(){return edge_1;};
	int getEdge2(){return edge_2;};
	int getEdge3(){return edge_3;};
};
SNode::SNode(){
	right=NULL;
	center=NULL;
	left=NULL;
	active=false;
	edge_1 = -9;
	edge_2 = -9;
	edge_3 = -9;
}

/*
//=======================================================
// ■　ThreeTree
//=======================================================
//   Voxel節点サーチ木の基本クラス
//=======================================================*/
class ThreeTree{
protected:
	int size_x;												/* メッシュのx向きの要素数 */
	int size_y;												/* メッシュのy向きの要素数 */
	int size_z;												/* メッシュのz向きの要素数 */
	SNode *root;											/* ルート節 */
	void deleteChildren(SNode **node);						/* 子供たち大破壊 */
public:
	ThreeTree(){;}											/* コンストラクタ */
	ThreeTree(int x, int y, int z);							/* コンストラクタ */
	~ThreeTree();											/* デストラクタ */
	int getSNodeID(int x, int y, int z);					/* 節、ｘ、ｙ、ｚのID取得 */
	bool make_new_node(int x, int y, int z);				/* 新たな節作成 */
	bool search(int x, int y, int z);						/* サーチ */
	bool search2(int x, int y, int z);						/* サーチその２　節がなければ探したうえで作成 */
};
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
ThreeTree::ThreeTree(int x, int y, int z){
	size_x = x;
	size_y = y;
	size_z = z;
	try{
		root = new SNode();
	}catch(bad_alloc){
		cout << "Memory Error!!\r\n";
		exit(1);
	}
}
/*//=======================================================
  // ● デストラクタ
  //=======================================================*/
ThreeTree::~ThreeTree(){
	deleteChildren(&root);
}
/*//=======================================================
  // ● 子供たち大破壊
  //=======================================================*/
void ThreeTree::deleteChildren(SNode **node){
	if(*((*node)->getLeftPtr()) != NULL){
		deleteChildren((*node)->getLeftPtr());
	}
	if(*((*node)->getCenterPtr()) != NULL){
		deleteChildren((*node)->getCenterPtr());
	}
	if(*((*node)->getRightPtr()) != NULL){
		deleteChildren((*node)->getRightPtr());
	}
	delete (*node);
}
/*//=======================================================
  // ● 新たな節作成
  //=======================================================*/
bool ThreeTree::make_new_node(int x, int y, int z){
	SNode **node = &root;
	for(int i = 0 ; i < x ; i++){
		node = (*node)->getLeftPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	for(int j = 0 ; j < y ; j++){
		node = (*node)->getCenterPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	for(int k = 0 ; k < z ; k++){
		node = (*node)->getRightPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	if( (*node)->isActive() ){
		return false;
	}else{
		(*node)->Activate();
		return true;
	}
}
/*//=======================================================
  // ● サーチ
  //=======================================================*/
bool ThreeTree::search(int x, int y, int z){
	SNode **node = &root;
	for(int i = 0 ; i < x ; i++){
		node = (*node)->getLeftPtr();
		if(*node == NULL){
			return false;
		}
	}
	for(int j = 0 ; j < y ; j++){
		node = (*node)->getCenterPtr();
		if(*node == NULL){
			return false;
		}
	}
	for(int k = 0 ; k < z ; k++){
		node = (*node)->getRightPtr();
		if(*node == NULL){
			return false;
		}
	}
	return true;
}
/*//=======================================================
  // ● サーチその２　節がなければ探したうえで作成
  //=======================================================*/
bool ThreeTree::search2(int x, int y, int z){
	bool bl = make_new_node(x, y, z);
	return( !bl );
}
/*//=======================================================
  // ● 節のｘ、ｙ、ｚのID取得
  //=======================================================*/
int ThreeTree::getSNodeID(int x, int y, int z){
	SNode **node = &root;
	for(int i = 0 ; i < x ; i++){
		node = (*node)->getLeftPtr();
	}
	for(int j = 0 ; j < y ; j++){
		node = (*node)->getCenterPtr();
	}
	for(int k = 0 ; k < z ; k++){
		node = (*node)->getRightPtr();
	}
	int id = (*node)->getID();
	return id;
}
/*
//=======================================================
// ■　SearchTree
//=======================================================
// 
//=======================================================*/
class SearchTree : public ThreeTree{
private:
	int count;
	SNode **ptr;
	void makeSearchTree();									/* サーチ木作成 */
	void deleteChildren(SNode *node);						/* 子供たち大破壊 */
	bool make_new_nodes(int x, int y, int z, int id);		/* 新たな節作成 */
public:
	SearchTree(int x, int y, int z);						/* コンストラクタ(与えるサイズは、要素数) */
	~SearchTree();
	SNode *getPtr(int i ){return *(ptr+i);};
};
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
SearchTree::SearchTree(int x, int y, int z) : ThreeTree(x,y,z){
	const int sss = (x+1)*(y+1)*(z+1);
	ptr = new SNode*[sss];
	for(int i = 0 ; i < sss ; i++){
		ptr[i] = NULL;
	}
	*ptr = root;
	count = 1;
	makeSearchTree();
}
SearchTree::~SearchTree(){
	delete[] ptr;
}
/*//=======================================================
  // ● サーチ木作成
  //=======================================================*/
void SearchTree::makeSearchTree(){
	root->setID(0);
	root->setX(0);
	root->setY(0);
	root->setZ(0);
	queue<int> Qx;
	queue<int> Qy;
	queue<int> Qz;
	Qx.push(1); 	Qy.push(0); 	Qz.push(0);
	Qx.push(0); 	Qy.push(1); 	Qz.push(0);
	Qx.push(0); 	Qy.push(0); 	Qz.push(1);
	ThreeTree temp_tree(size_x, size_y, size_z);
	temp_tree.search2(0,0,0);
	temp_tree.search2(1,0,0);
	temp_tree.search2(0,1,0);
	temp_tree.search2(0,0,1);
	int xx,yy,zz;
	int temp;
	int id = 1;
	while(!Qx.empty()){
		/* キューから１つ取り出す */
		xx = Qx.front();
		Qx.pop();
		yy = Qy.front();
		Qy.pop();
		zz = Qz.front();
		Qz.pop();
		/* 取りだした節にIDセット */
		bool bl = make_new_nodes(xx, yy, zz, id);		
		if(bl) id++;
		/* キューに次の探索代入 */
		if(xx != size_x){
			temp =xx + 1;
			if(!temp_tree.search2(temp, yy, zz)){
				Qx.push(temp);
				Qy.push(yy);
				Qz.push(zz);
			}
		}
		if(yy != size_y){
			temp = yy + 1;
			if(!temp_tree.search2(xx, temp, zz)){
				Qx.push(xx);
				Qy.push(temp);
				Qz.push(zz);
			}
		}
		if(zz != size_z){
			temp = zz + 1;
			if(!temp_tree.search2(xx, yy, temp)){
				Qx.push(xx);
				Qy.push(yy);
				Qz.push(temp);
			}
		}
	}
}
/*//=======================================================
  // ● 新たな節作成
  //=======================================================*/
bool SearchTree::make_new_nodes(int x, int y, int z, int id){
	SNode **node = &root;
	for(int i = 0 ; i < x ; i++){
		node = (*node)->getLeftPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	for(int j = 0 ; j < y ; j++){
		node = (*node)->getCenterPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	for(int k = 0 ; k < z ; k++){
		node = (*node)->getRightPtr();
		try{
			if(*node == NULL){
				*node = new SNode();
			}
		}catch(bad_alloc){
			cout << "Memory Error!!\r\n";
			exit(1);
		}
	}
	if((*node)->isActive()){
		return false;
	}else{
		(*node)->setID(id);
		(*node)->Activate();
		(*node)->setX(x);
		(*node)->setY(y);
		(*node)->setZ(z);		
		*(ptr+count) = (*node);
		count++;
		return true;
	}
}

