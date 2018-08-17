#include <stdio.h>
#include <stdlib.h>

//This code assumes the mesh is in AWP format (fast y, x, z)
//This will perform bilinear interpolation (but not tri)

int main(int argc, char** argv) {
	if (argc<8) {
		printf("Usage: %s <mesh file> <nx> <ny> <nz> <x-index> <y-index> <z-index>\n", argv[0]);
		exit(1);
	}

	char* mesh_filename = argv[1];
	int nx = atoi(argv[2]);
	int ny = atoi(argv[3]);
	int nz = atoi(argv[4]);
	float x_index = atof(argv[5]);
	float y_index = atof(argv[6]);
	int z_index = atoi(argv[7]);
	float data[2][2][3];

	int x_grid_pts[2];
	int y_grid_pts[2];

	x_grid_pts[0] = (int)x_index;
	x_grid_pts[1] = x_grid_pts[0]+1;
	y_grid_pts[0] = (int)y_index;
	y_grid_pts[1] = y_grid_pts[0]+1;


	FILE* fp_in = fopen(mesh_filename, "rb");
	//Fast y, x, z
	int i, j;
	for (i=0; i<2; i++) {
		for (j=0; j<2; j++) {
			size_t offset = 3*sizeof(float)*(((size_t)z_index)*nx*ny + x_grid_pts[i]*ny + y_grid_pts[j]);
			fseek(fp_in, offset, SEEK_SET);
			fread(data[i][j], sizeof(float), 3, fp_in);
		}
	}
	fclose(fp_in);
	
	//Perform interpolation
	float x_frac = x_index - x_grid_pts[0];
	float y_frac = y_index - y_grid_pts[0];
	//P1: interp btw (x0,y0) and (x0,y1), P2: (x1,y0) and (x1,y1), P3: P1 and P2
	float p1[3], p2[3], p3[3];
	
	for (i=0; i<3; i++) {
		p1[i] = (1-y_frac)*data[0][0][i] + y_frac*data[0][1][i];
		p2[i] = (1-y_frac)*data[1][0][i] + y_frac*data[1][1][i];
		p3[i] = (1-x_frac)*p1[i] + x_frac*p2[i];
	}

	printf("%f %f %f\n", p3[0], p3[1], p3[2]);

	return 0;

}
