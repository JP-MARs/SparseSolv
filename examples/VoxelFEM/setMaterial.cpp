
#include "VoxelMesh.h"


/*//=======================================================
  // ● 各要素の材料情報作成
  //=======================================================*/
void VoxelMesh::setMaterials(){
	double air_perm = Air_perm ;
	double iron_perm = Iron_perm ;

	Element *temp = elements;
	double g_xyz[3];
	for(int i = 0 ; i < element_num ; i++){
		temp->getGrav(g_xyz);
		bool in_coil1 = (g_xyz[0] < COIL_MAX && g_xyz[1] < COIL_MAX);
		bool in_coil2 = (g_xyz[0] < COIL_MIN && g_xyz[1] < COIL_MIN);
		//Iron
		if(g_xyz[0] < COIL_MIN && g_xyz[1] < COIL_MIN && g_xyz[2] < IRON_Z){
			temp->setMaterial(MAT_IRON, Iron_perm);
		// Coil
		}else if(in_coil1 && !in_coil2 && g_xyz[2] < COIL_Z){
			temp->setMaterial(MAT_COIL, Air_perm);
		// Other
		}else{
			temp->setMaterial(MAT_AIR, Air_perm);
		}
		temp++;
	}
}

/*//=======================================================
  // ● 磁気抵抗率ゲット
  //=======================================================*/
void VoxelMesh::calc_v_tensor(double *v_tensor, Element *Ele_ptr){
	/* まず透磁率を取得 */
	double mat_permi = Ele_ptr->getMatPermeability();
	double v = 1.0 / mat_permi;

	v_tensor[0] = v; v_tensor[1] = 0; v_tensor[2] = 0;
	v_tensor[3] = 0; v_tensor[4] = v; v_tensor[5] = 0;
	v_tensor[6] = 0; v_tensor[7] = 0; v_tensor[8] = v;
	return;
}


