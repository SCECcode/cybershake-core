#include <stdio.h>
#include <stdlib.h>

#include "string.h"
#include "mpi.h"
#include <time.h>

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
	double tot[3] = {0.0, 0.0, 0.0};
	int num_vals = 0;
	//Calclate average using diamond
	for (i=x_start; i<=x_stop; i++) {
		int j_start = max(y_start+abs(i-x), 0);
		int j_end = min(y_stop-abs(i-x), ny-1);
		for (j=j_start; j<=j_end; j++) {
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

void process_mesh(int my_id, char* mesh_in, char* pts_file, int nx, int ny, int nz, int starting_slice, int ending_slice, int smoothing_dist, char* mesh_out) {
	printf("%d) processing slices %d to %d.\n", my_id, starting_slice, ending_slice);
	//FILE* fp_in = fopen(mesh_in, "rb");
	//FILE* fp_out = fopen(mesh_out, "wb");
	MPI_File fp_in, fp_out;
	int rc = MPI_File_open(MPI_COMM_WORLD, mesh_in, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_in);
	if (rc!=0) {
		fprintf(stderr, "Error opening file %s for reading.\n", mesh_in);
		MPI_Finalize();
		exit(2);
	}
	rc = MPI_File_open(MPI_COMM_WORLD, mesh_out, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp_out);
	if (rc!=0) {
	        fprintf(stderr, "Error opening file %s for writing.\n", mesh_out);
                MPI_Finalize();
                exit(2);
        }

	//For now, let everyone read pts_in; if an issue try broadcasting instead
	FILE* pts_in = fopen(pts_file, "r");

	float* horizontal_slice_in = malloc(sizeof(float)*3*nx*ny);
	float* horizontal_slice_out = malloc(sizeof(float)*3*nx*ny);

	//Get points
	printf("%d) Getting points.\n", my_id);
	fflush(stdout);
	int pt_x, pt_y;
	float lon, lat;
	float vals[3];
	int num_pts = 0;
	//int points[1000000][2];
	int** points;
	int i;
	int max_points = 100000;
	points = malloc(sizeof(int*)*max_points);
	for (i=0; i<max_points; i++) {
		points[i] = malloc(sizeof(int)*2);
	}
	while (fscanf(pts_in, "%d %d %f %f\n", &pt_x, &pt_y, &lon, &lat)!=EOF) {
		points[num_pts][0] = pt_x;
		points[num_pts][1] = pt_y;
		num_pts++;
		if (num_pts==max_points) {
			//expand array
			printf("%d) expanding array\n", my_id);
			max_points += 100000;
			points = realloc(points, sizeof(int*)*max_points);
			for (i=num_pts; i<max_points; i++) {
				points[i] = malloc(sizeof(int)*2);
			}
		}
	}
	fclose(pts_in);
	printf("%d) done.\n", my_id);

	int point_index = 0;
	int j, k;
	//Iterate over horizontal slices
	//Fast y, x
	//for (i=0; i<nz; i++) {

	for (i=starting_slice; i<ending_slice; i++) {
		time_t now = time(NULL);
		char* time_str = ctime(&now);
		printf("%d @ %s) Horizontal slice %d of %d.\n", my_id, time_str, (i-starting_slice), (ending_slice-starting_slice));
		point_index = 0;
		//fread(horizontal_slice_in, sizeof(float), 3*nx*ny, fp_in);
		MPI_Offset offset = sizeof(float)*3*nx*ny*i;
		MPI_Status status;
		MPI_File_read_at(fp_in, offset, horizontal_slice_in, 3*nx*ny, MPI_FLOAT, &status);
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
		//fwrite(horizontal_slice_out, sizeof(float), 3*nx*ny, fp_out);
		now = time(NULL);
                time_str = ctime(&now);
		printf("%d @ %s) Writing slice %d to file.\n", my_id, time_str, i-starting_slice);
		MPI_File_write_at(fp_out, offset, horizontal_slice_out, 3*nx*ny, MPI_FLOAT, &status);
	}
	//fflush(fp_out);
	//fclose(fp_out);
	//fclose(fp_in);
	MPI_File_close(&fp_in);
	MPI_File_close(&fp_out);
	for (i=0; i<max_points; i++) {
		free(points[i]);
	}
	free(points);
}


int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

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

	//Figure out what horizontal slices you're responsible for
	int my_id, num_procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if (num_procs>nz) {
		if (my_id==0) {
			fprintf(stderr, "Can't use more processors than z-slices.  nz=%d but num_procs=%d, aborting.\n", nz, num_procs);
		}
		MPI_Finalize();
		exit(1);
	}

	float slices_per_proc = ((float)nz/((float)num_procs));
	int starting_slice = (int)(my_id*slices_per_proc);
	int ending_slice = (int)((my_id+1)*slices_per_proc);
	if (ending_slice>nz) {
		ending_slice = nz;
	}
	printf("%d) Responsible for slices %d - %d\n", my_id, starting_slice, ending_slice);
	fflush(stdout);

	process_mesh(my_id, mesh_in, pts_file, nx, ny, nz, starting_slice, ending_slice, smoothing_dist, mesh_out);

	MPI_Finalize();
	return 0;
}
