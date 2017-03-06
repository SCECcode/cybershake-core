#include <stdio.h>
#include <stdlib.h>

#include "string.h"

inline int min(int a, int b) {
	if (a<b) 
		return a;
	return b;
}

inline int max(int a, int b) {
	if (a>b) 
		return a;
	return b;
}

void average_point(float* slice, int x, int y, int nx, int ny, int smooth_dist, float (*vals)[3]) {
	int i, j, k;
	int x_start, x_stop, y_start, y_stop;
	//X-dir
	x_start = max(x-smooth_dist, 0);
	x_stop = min(x+smooth_dist, nx-1);
        y_start = max(y-smooth_dist, 0);
        y_stop = min(y+smooth_dist, ny-1);
	float tot[3] = {0, 0, 0};
	int num_vals = 0;
	//Calclate average using diamond
	for (i=x_start; i<=x_stop; i++) {
		for (j=max(y_start+abs(i-x), 0); j<=min(y_stop-abs(i-x), ny-1); j++) {
			int index = 3*(i*ny + j);
			for (k=0; k<3; k++) {
				tot[k] += slice[index+k];
			}
			num_vals++;	
		}
	}
	for (i=0; i<3; i++) {
		(*vals)[i] = tot[i]/((float)num_vals);
	}
}

void process_mesh(char* mesh_in, char* pts_file, int nx, int ny, int nz, int smoothing_dist, char* mesh_out) {
	FILE* fp_in = fopen(mesh_in, "rb");
	FILE* fp_out = fopen(mesh_out, "wb");
	FILE* pts_in = fopen(pts_file, "r");

	float* horizontal_slice_in = malloc(sizeof(float)*3*nx*ny);
	float* horizontal_slice_out = malloc(sizeof(float)*3*nx*ny);

	//Get points
	int pt_x, pt_y;
	float lon, lat;
	float vals[3];
	int num_pts = 0;
	//int points[1000000][2];
	int** points;
	int i;
	int MAX_PTS = 4000000;
	points = malloc(sizeof(int*)*MAX_PTS);
	for (i=0; i<MAX_PTS; i++) {
		points[i] = malloc(sizeof(int)*2);
	}
	while (fscanf(pts_in, "%d %d %f %f\n", &pt_x, &pt_y, &lon, &lat)!=EOF) {
		points[num_pts][0] = pt_x;
		points[num_pts][1] = pt_y;
		num_pts++;
	}
	fclose(pts_in);

	int point_index = 0;
	int j, k;
	//Iterate over horizontal slices
	//Fast y, x
	for (i=0; i<nz; i++) {
		printf("Horizontal slice %d of %d.\n", i, nz);
		point_index = 0;
		fread(horizontal_slice_in, sizeof(float), 3*nx*ny, fp_in);
		for (j=0; j<nx; j++) {
			//Copy X-strip over, will overwrite with new values
			memcpy(horizontal_slice_out+3*j*ny, horizontal_slice_in+3*j*ny, 3*sizeof(float)*ny);
			//any points with this X-value?
			while (point_index<num_pts && points[point_index][0]==j) {
                                //Do averaging for point
                                average_point(horizontal_slice_in, points[point_index][0], points[point_index][1], nx, ny, smoothing_dist, &vals);
                                //Update point
                                for (k=0; k<3; k++) {
                                        horizontal_slice_out[3*(j*ny + points[point_index][1]) + k] = vals[k];
                                }
                                //Get next point
				point_index++;
			}
		}
		//Write slice to file
		fwrite(horizontal_slice_out, sizeof(float), 3*nx*ny, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);
	for (i=0; i<MAX_PTS; i++) {
		free(points[i]);
	}
	free(points);
}


int main(int argc, char** argv) {
	if (argc<8) {
		printf("Usage: %s <AWP-format mesh in> <list of smoothing pts> <RWG nx> <RWG ny> <nz> <smoothing range in pts> <velocity mesh out>", argv[0]);
		exit(1);
	}

	char* mesh_in = argv[1];
	char* pts_file = argv[2];
	int nx = atoi(argv[3]);
	int ny = atoi(argv[4]);
	int nz = atoi(argv[5]);
	int smoothing_dist = atoi(argv[6]);
	char* mesh_out = argv[7];

	process_mesh(mesh_in, pts_file, nx, ny, nz, smoothing_dist, mesh_out);

	return 0;
}
