/*
//=======================================================
// ■ Node_Mesh
//=======================================================
// FEMの節点定義クラス（メッシュ生成用）
//=======================================================*/
class Node_Mesh {
private:
	int ID;												/* このノードのID */
	int x;												/* X座標 */
	int y;												/* Y座標 */
	int z;												/* Z座標 */
	double pos_x;										/* ホントのX座標 */
	double pos_y;										/* ホントのY座標 */
	double pos_z;										/* ホントのZ座標 */
public:
	Node_Mesh(){;};													/* コンストラクタ */
	Node_Mesh(int xx, int yy, int zz);								/* コンストラクタ2 */
	void set(int xx, int yy, int zz);								/* 座標セット */
	void set_pos(double xx, double yy, double zz);					/* ホントの座標セット */
	int getX(){return x;};
	int getY(){return y;};
	int getZ(){return z;};
	double getPosX(){return pos_x;};
	double getPosY(){return pos_y;};
	double getPosZ(){return pos_z;};
	int getID(){return ID;};
	void setID(int id){ID=id;};
};
/*//=======================================================
  // ● コンストラクタ
  //=======================================================*/
Node_Mesh::Node_Mesh(int xx, int yy, int zz){
	set(xx, yy, zz);
}
/*//=======================================================
  // ● 節点座標セッタ
  //=======================================================*/
void Node_Mesh::set(int xx, int yy, int zz){
	x = xx; y = yy; z = zz;
}
/*//=======================================================
  // ● ホントの座標セット
  //=======================================================*/
void Node_Mesh::set_pos(double xx, double yy, double zz){
	pos_x = xx; pos_y = yy; pos_z = zz;
}

