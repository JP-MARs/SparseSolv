#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <queue>

using namespace std;

#include "../Defines_Mesh.h"
#include "../Defines_Other.h"
#include "SearchTree.h"
#include "Node_Mesh.h"
#include "Edge.h"
#include "Element.h"
#include "Boundary.h"
#include "T_Vector.h"
#include "MeshMaker.h"


/*//=======================================================
  // ●　メイン関数
  //=======================================================*/
int main(int argc, char *argv[]){
/*
	int x_size = 15;
	int y_size = 15;
	int z_size = 15;

	double length1[] = {0.1, 0.1, 0.1};
	double length2[] = {0.5, 0.5, 0.5};

	int change[] = {5, 5, 5};
*/

	int x_size = TOTAL_MESH_X;
	int y_size = TOTAL_MESH_Y;
	int z_size = TOTAL_MESH_Z;

	double length1[] = {X_L1, Y_L1, Z_L1};
	double length2[] = {X_L2, Y_L2, Z_L2};

	int change[] = {MESH_X1, MESH_Y1, MESH_Z1};

	MeshMaker mesh(x_size, y_size, z_size, length1, length2, change);

	return 0;

}
