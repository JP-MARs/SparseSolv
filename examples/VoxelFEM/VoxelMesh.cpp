
#include "VoxelMesh.h"


/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
VoxelMesh::VoxelMesh(double *xyz_L1, double *xyz_L2, int *change_xyz){
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
	/* Nodeデータ読み取り */
	setNodeData();
	/* Edgeデータ読み取り */
	setEdgeData();
	/* Elementデータ読み取り */
	setElementData();
	Element::setNodeData(NodeData);
	Element::setEdgeData(EdgeData);
	/* 辺のサブIDセット */
}

VoxelMesh::~VoxelMesh(){
	delete[] NodeData;
	delete[] EdgeData;
	delete[] elements;
}
/*//=======================================================
  // ● Nodeデータ読み取り
  //=======================================================*/
void VoxelMesh::setNodeData(){
	FILE *fp = fopen("./MeshData/Node.vxldata", "r");
	/* 辺の数受け取り */
	fscanf(fp, "%d %d %d %d\n", &node_num, &x_num, &y_num, &z_num);
	NodeData = new Node[node_num];
	int id, x, y, z;
	double vx,vy,vz;
	Node *temp = NodeData;
	for(int i = 0 ; i < node_num ; i++){
		fscanf(fp, "%d %lf %lf %lf\n", &id, &vx, &vy, &vz);
		temp->set(id, vx, vy, vz);
		temp++;
	}
	fclose(fp);
}
/*//=======================================================
  // ● Edgeデータ読み取り
  //=======================================================*/
void VoxelMesh::setEdgeData(){
	FILE *fp = fopen("./MeshData/Edge.vxldata", "r");
	/* 辺の数受け取り */
	fscanf(fp, "%d\n", &edge_num);
	EdgeData = new Edge[edge_num];
	int id, s_id, e_id;
	Edge *temp = EdgeData;
	for(int i = 0 ; i < edge_num ; i++){
		fscanf(fp, "%d %d %d\n", &id, &s_id, &e_id);
		temp->setID(id);
		temp->set(s_id, e_id);
		temp++;
	}
	fclose(fp);
}
/*//=======================================================
  // ● Elementデータ読み取り
  //=======================================================*/
void VoxelMesh::setElementData(){
	FILE *fp = fopen("./MeshData/Element.vxldata", "r");
	/* 辺の数受け取り */
	fscanf(fp, "%d\n", &element_num);
	elements = new Element[element_num];
	int id, x, y, z;
	int n[8], e[12];
	bool bl;
	Element *temp = elements;
	for(int i = 0 ; i < element_num ; i++){
		fscanf(fp, "%d %d %d %d", &id, &x, &y, &z);
		for(int j = 0 ; j < 8 ; j++){
			fscanf(fp, "%d ", n+j);
		}
		bl = false;
		for(int j = 0 ; j < 12 ; j++){
			fscanf(fp, "%d ", e+j);
			bl |= EdgeData[e[j]].isT_Edge();
		}
		temp->set(x, y, z, n, e, id);
		temp++;
	}
	fclose(fp);
}

/*//=======================================================
  // ● SubIDセット
  //=======================================================*/
void VoxelMesh::setSubIDs(){
	/* 通常 */
	int count_x=0;
	Edge *temp_E = EdgeData;
	vector<int> zero_edges;
	for(int i = 0 ; i < edge_num ; i++){
		/* ゼロ辺なら次へ */
		bool bl  = temp_E->isZeroEdge();
		if(bl){
			zero_edges.push_back(i);
			temp_E++;
			continue;
		}
		/* 既知辺だけ、行列用のサブIDをセット */
		temp_E->setSubID(count_x);
		count_x++;
		temp_E++;
	}
	for(int i = 0 ; i < zero_edges.size() ; i++){
		EdgeData[zero_edges[i]].setSubID(count_x);
		count_x++;
	}
}

/*//=======================================================
  // ● 結果書き出し
  //=======================================================*/
void VoxelMesh::writeVTK(double* resultA){
	fstream fp("./post_result.vtk", std::ios::out);
    /* ヘッダ書き出し */
    fp << "# vtk DataFile Version 2.0" << endl;
    fp << "hexahedron" << endl;
    fp << "ASCII" << endl;
    fp << "DATASET UNSTRUCTURED_GRID" << endl;

    /* 節点書き出し */
    fp << "POINTS " << node_num << " float" << endl;
    for(int mm = 0 ; mm < node_num ;mm++){
        double xx = NodeData[mm].getX();
        double yy = NodeData[mm].getY();
        double zz = NodeData[mm].getZ();
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
		int mat = elements[i].getMatType();
        fp << mat << endl;
    }

    /*------------------------------------------------*/
    /* 磁束密度書き出し */
    /*-----------------------------*/
	if(resultA != nullptr){
		fp << "VECTORS FluxDensityVec float" << endl;
		for(int i = 0 ; i < element_num ; i++){
			double pos[3] = {0,0,0};
			double vecB[3];
			elements[i].calcB(vecB, resultA, pos);
			fp << vecB[0] << " " << vecB[1] << " " << vecB[2] << endl;
		}
	}

    /*------------------------------------------------*/
    /* 電流密度書き出し */
    /*-----------------------------*/
	fp << "VECTORS CurrentDensityVec float" << endl;
    for(int i = 0 ; i < element_num ; i++){
		double pos[3] = {0,0,0};
		double vecJ[3];
		elements[i].calcJ(vecJ, pos);
		fp << vecJ[0] << " " << vecJ[1] << " " << vecJ[2] << endl;
    }
	
    fp.close();	
}
