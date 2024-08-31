/*
//=======================================================
// ■ MeshMaker
//=======================================================
// FEMのボクセルメッシュ自動生成プログラム
//=======================================================*/
class MeshMaker {
private:
	int x_num;												/* ｘ方向の要素数 */
	int y_num;												/* ｙ方向の要素数 */
	int z_num;												/* ｚ方向の要素数 */
	int node_num;											/* 節点の総数 */
	int edge_num;											/* 辺の総数 */
	int element_num;										/* 要素の総数 */
	double x_length1;										/* 各要素のｘ長さ */
	double y_length1;										/* 各要素のｙ長さ */
	double z_length1;										/* 各要素のｚ長さ */
	double x_length2;										/* 各要素のｘ長さ、その２ */
	double y_length2;										/* 各要素のｙ長さ、その２ */
	double z_length2;										/* 各要素のｚ長さ、その２ */
	int change_x;											/* 要素のサイズ変更になる位置Ｘ */
	int change_y;											/* 要素のサイズ変更になる位置Ｙ */
	int change_z;											/* 要素のサイズ変更になる位置Ｚ */
	Node_Mesh *NodeData;									/* 節点データ保持配列 */
	Edge *EdgeData;											/* 辺データ保持配列 */
	Element *elements;										/* 要素データ保持配列 */
	int A0_size;											/* ゼロ既知になる辺のサイズ */
	int T_size;												/* 電流ポテンシャルを持つ辺のサイズ */
	int *A0_ID;												/* ゼロ既知になる辺ID保持配列 */
	int *T_ID;												/* 電流ポテンシャルを持つ辺ID保持配列 */

	void setNodeNum(SearchTree *S_tree);					/* 節点番号作成 */
	void setEdgeNum(SearchTree *S_tree);					/* 辺番号作成 */
	void setElementNum(SearchTree *S_tree);					/* 要素番号作成 */
	/* サブルーチン群 */
	/* ------------------------------------------------------------------------------------- */
	void setElementSubRtin(int x, int y, int z, SearchTree *S_tree);
	void setEdgeSubRutin(int i, SearchTree *S_tree, int dir, int& count);
	/* ------------------------------------------------------------------------------------- */
public:
	MeshMaker(int xn, int yn, int zn, double *xyz_L1, double *xyz_L2, int *change_xyz);		/* コンストラクタ */
	~MeshMaker();
	int getNodeNum(){return node_num;};
	int getEdgeNum(){return edge_num;};
	int getElementNum(){return element_num;};
	void FileWrite();																		/* ファイル書き出し */
	void WriteVTK();
};
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
MeshMaker::MeshMaker(int xn, int yn, int zn, double *xyz_L1, double *xyz_L2, int *change_xyz){
	/* 各種サイズ設定 */
	x_length1 = xyz_L1[0];
	y_length1 = xyz_L1[1];
	z_length1 = xyz_L1[2];
	x_length2 = xyz_L2[0];
	y_length2 = xyz_L2[1];
	z_length2 = xyz_L2[2];
	change_x = change_xyz[0];
	change_y = change_xyz[1];
	change_z = change_xyz[2];
	x_num = xn;
	y_num = yn;
	z_num = zn;
	node_num = (x_num+1)*(y_num+1)*(z_num+1);
	edge_num = ( x_num*(y_num+1) + (x_num+1)*y_num )*(z_num+1) + z_num*(x_num+1)*(y_num+1);
	element_num = x_num*y_num*z_num;
	A0_size = 0;
	T_size = 0;
	/* 節点作成 */
	SearchTree S_tree(x_num, y_num, z_num);
	setNodeNum(&S_tree);
	/* 辺作成 */
	setEdgeNum(&S_tree);
	/* 要素作成 */
	setElementNum(&S_tree);
	Element::setNodeData(NodeData);
	Element::setEdgeData(EdgeData);
	/* ファイル書き出し */
	FileWrite();
	WriteVTK();
}
MeshMaker::~MeshMaker(){
	delete[] NodeData;
	delete[] EdgeData;
	delete[] elements;
	delete[] A0_ID;
	delete[] T_ID;
}
/*//=======================================================
  // ● ファイル書き出し
  //=======================================================*/
void MeshMaker::FileWrite(){
	FILE *fp_n = fopen("../MeshData/Node.vxldata", "w");
	FILE *fp_ed = fopen("../MeshData/Edge.vxldata", "w");
	FILE *fp_el = fopen("../MeshData/Element.vxldata", "w");
	int ID;
	int x, y, z;
	/* ノードデータ書き出し */
	Node_Mesh *tempN = NodeData;
	fprintf(fp_n, "%d %d %d %d\n", node_num, x_num, y_num, z_num);
	for(int i = 0 ; i < node_num ; i++){
		ID = tempN->getID();
		//x = tempN->getX(); y = tempN->getY();  z = tempN->getZ(); 
		double xx = tempN->getPosX();
		double yy = tempN->getPosY();
		double zz = tempN->getPosZ();
		//fprintf(fp_n, "%d %d %d %d\n", ID, x, y, z);
		fprintf(fp_n, "%d %lf %lf %lf\n", ID, xx, yy, zz);
		tempN++;
	}
	/* 辺データ書き出し */
	Edge *tempEd = EdgeData;
	int sID, eID, A0, T;
	fprintf(fp_ed, "%d\n", edge_num);
	for(int i = 0 ; i < edge_num ; i++){
		ID = tempEd->getID(); sID = tempEd->getSID(); eID = tempEd->getEID();
		A0 = tempEd->isZeroEdge(); T = tempEd->isT_Edge();
		fprintf(fp_ed, "%d %d %d\n", ID, sID, eID);
		tempEd++;
	}
	/* 要素データ書き出し */
	Element *tempEl = elements;
	int *Ns, *Es, xx, yy, zz;
	fprintf(fp_el, "%d\n", element_num);
	for(int i = 0 ; i < element_num ; i++){
		ID = tempEl->getID(); 
		xx = tempEl->getX(); yy = tempEl->getY(); zz = tempEl->getZ();
		Ns = tempEl->getNodes(); Es = tempEl->getEdges();
		fprintf(fp_el, "%d %d %d %d ", ID, xx, yy, zz);
		for(int j = 0 ; j < 8 ; j++){
			fprintf(fp_el, "%d ", Ns[j]);
		}
		for(int j = 0 ; j < 12 ; j++){
			fprintf(fp_el, "%d ", Es[j]);
		}
		fprintf(fp_el, "\n");
		tempEl++;
	}
	fclose(fp_n);
	fclose(fp_ed);
	fclose(fp_el);
}
/*//=======================================================
  // ● 節点の番号作成～！
  //=======================================================*/
void MeshMaker::setNodeNum(SearchTree *S_tree){
	/* 節点データ作成 */
	NodeData = new Node_Mesh[node_num];
	double vx,vy,vz;
	Node_Mesh *tempPtr = NodeData;
	for(int i = 0 ; i < node_num ; i++){
		SNode *temp_n = S_tree->getPtr(i);
		int temp_x = temp_n->getX();
		int temp_y = temp_n->getY();
		int temp_z = temp_n->getZ();
		tempPtr->set(temp_x, temp_y, temp_z);
		tempPtr->setID(i);
		vx = ( temp_x > change_x ? x_length1*change_x + (temp_x-change_x)*x_length2 : temp_x*x_length1 );
		vy = ( temp_y > change_y ? y_length1*change_y + (temp_y-change_y)*y_length2 : temp_y*y_length1 );
		vz = ( temp_z > change_z ? z_length1*change_z + (temp_z-change_z)*z_length2 : temp_z*z_length1 );
		tempPtr->set_pos(vx, vy, vz);
		tempPtr++;
	}
}
/*//=======================================================
  // ● 辺の番号作成～！
  //=======================================================*/
void MeshMaker::setEdgeNum(SearchTree *S_tree){
	/* 辺作成 */
	EdgeData = new Edge[edge_num];
	int count = 0;
	for(int i = 0 ; i < node_num ; i++){
		setEdgeSubRutin(i, S_tree, 0, count);
		setEdgeSubRutin(i, S_tree, 1, count);
		setEdgeSubRutin(i, S_tree, 2, count);
	}
	A0_ID = new int[A0_size];
	T_ID = new int[T_size];
	Edge *tempE = EdgeData;
	/* 既知のポテンシャルと、Tを持つ辺IDを記憶 */
	int tt=0, tt2=0;
	for(int i = 0 ; i < edge_num ; i++){
		if( tempE->isZeroEdge() ){
			A0_ID[tt] = tempE->getID();
			tt++;
		}
		if( tempE->isT_Edge() ){
			T_ID[tt2] = tempE->getID();
			tt2++;
		}
		tempE++;
	}
}
/*//=======================================================
  // ● 辺番号設定用サブルーチン
  //=======================================================*/
void MeshMaker::setEdgeSubRutin(int i, SearchTree *S_tree, int dir, int& count){
	int id1,id2;
	int xx1,xx2,yy1,yy2,zz1,zz2;
	/* 始点代入 */
	Node_Mesh *tempN = NodeData+i;
	id1 = tempN->getID();
	xx1 = tempN->getX(); yy1 = tempN->getY(); zz1 = tempN->getZ();
	bool flag = false;
	xx2 = xx1; yy2 = yy1; zz2 = zz1;
	/* 終点代入 */
	if(dir == 0){
		if(xx1 != x_num){
			xx2++;
			flag = true;
		}
	}else if(dir == 1){
		if(yy1 != y_num){
			yy2++;
			flag = true;
		}
	}else{
		if(zz1 != z_num){
			zz2++;
			flag = true;
		}
	}
	/* 辺データ代入 */
	SNode *temp_node;
	if(flag){
		/* 辺にIDを与え、始点・終点をセット */
		Edge *tempE = EdgeData+count;
		id2 = S_tree->getSNodeID(xx2, yy2, zz2);
		tempE->set(id1, id2);
		tempE->setID(count);
		temp_node = S_tree->getPtr(id1);
		/* その節点から出ている辺3本のID定義 */
		if(dir == 0 && xx1 != x_num){
			temp_node->setEdge1(count);
		}else if(dir == 1 && yy1 != y_num){
			temp_node->setEdge2(count);
		}else if(dir == 2 && zz1 != z_num){
			temp_node->setEdge3(count);
		}
		count++;
	}
}
/*//=======================================================
  // ● 要素番号作成
  //=======================================================*/
void MeshMaker::setElementNum(SearchTree *S_tree){
	elements = new Element[element_num];
	for(int k = 0 ; k < z_num ; k++){
		for(int j = 0 ; j < y_num ; j++){
			for(int i = 0 ; i < x_num ; i++){
				setElementSubRtin(i, j, k, S_tree);
			}
		}
	}
}
/*//=======================================================
  // ● 要素番号作成のサブルーチン
  //=======================================================*/
void MeshMaker::setElementSubRtin(int x, int y, int z, SearchTree *S_tree){
	const int E_ID = x + y*x_num + z * (x_num*y_num);
	/* 辺と節点のIDを全部取得！ */
	int id0 = S_tree->getSNodeID(x, y, z);
	SNode *node = S_tree->getPtr(id0);
	int edge1 = node->getEdge1();
	int id1 = (EdgeData+edge1)->getEID();
	int edge2 = node->getEdge2();
	int id2 = (EdgeData+edge2)->getEID();
	int edge3 = node->getEdge3();
	int id3 = (EdgeData+edge3)->getEID();
	node = S_tree->getPtr(id1);
	int edge4 = node->getEdge2();
	int id4 = (EdgeData+edge4)->getEID();
	int edge5 = node->getEdge3();
	int id5 = (EdgeData+edge5)->getEID();
	node = S_tree->getPtr(id2);
	int edge6 = node->getEdge1();
	int edge7 = node->getEdge3();
	int id6 = (EdgeData+edge7)->getEID();
	node = S_tree->getPtr(id3);
	int edge8 = node->getEdge1();
	int edge9 = node->getEdge2();
	node = S_tree->getPtr(id4);
	int edge10 = node->getEdge3();
	int id7 = (EdgeData+edge10)->getEID();
	node = S_tree->getPtr(id5);
	int edge11 = node->getEdge2();
	node = S_tree->getPtr(id6);
	int edge12 = node->getEdge1();
	/* 配列にする */
	int ns[] = {id0, id1, id4, id2,  id3, id5, id7, id6};
	int es[] = {edge12, edge8, edge1, edge6,  edge11, edge4, edge2, edge9,  edge10, edge7, edge3, edge5};
	/* 要素に代入 */
	(elements+E_ID)->set(x, y, z, ns, es, E_ID);
}

/*//=======================================================
  // ● VTK作成
  //=======================================================*/
void MeshMaker::WriteVTK(){
	fstream fp("./mesh.vtk", std::ios::out);
    /* ヘッダ書き出し */
    fp << "# vtk DataFile Version 2.0" << endl;
    fp << "hexahedron" << endl;
    fp << "ASCII" << endl;
    fp << "DATASET UNSTRUCTURED_GRID" << endl;

    /* 節点書き出し */
    fp << "POINTS " << node_num << " float" << endl;
    for(int mm = 0 ; mm < node_num ;mm++){
        double xx = NodeData[mm].getPosX();
        double yy = NodeData[mm].getPosY();
        double zz = NodeData[mm].getPosZ();
        fp << xx << " " << yy << " " << zz << endl;
    }
    /* 要素データの書き出し(要素数＋書き出す情報の数) */
    fp << "CELLS " << element_num << " " << (element_num*9) << endl;
    for(int i = 0 ; i < element_num ; i++){
        int inode[8];
		inode[0] = elements[i].getNodeID(0);
		inode[1] = elements[i].getNodeID(1);
		inode[2] = elements[i].getNodeID(2);
		inode[3] = elements[i].getNodeID(3);
		inode[4] = elements[i].getNodeID(4);
		inode[5] = elements[i].getNodeID(5);
		inode[6] = elements[i].getNodeID(6);
		inode[7] = elements[i].getNodeID(7);
        fp << "8 " << inode[0] << " " << inode[1] << " " << inode[2] <<  " " << inode[3] <<   " " << inode[4] <<   " " << inode[5] <<   " " << inode[6] <<   " " << inode[7] << endl;
    }
    /* 要素タイプの指定 */
    fp << "CELL_TYPES " << element_num << endl;
    for(int i = 0 ; i < element_num ; i++){
		fp << 12 << endl;
    }
    /* 材料番号の書き出し */
    fp << "CELL_DATA " << element_num << endl;
    fp << "SCALARS Material float" << endl;
    fp << "LOOKUP_TABLE default" << endl;
    for(int i = 0 ; i < element_num ; i++){
        fp << 1 << endl;
    }

    fp.close();


}