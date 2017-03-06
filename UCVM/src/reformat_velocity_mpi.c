#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "mpi.h"

void *check_malloc(size_t len);

int my_id, num_procs;

int X = 0;
int Y = 1;
int Z = 2;

int read_RWG_files(int n[], char* prefix, float*** data, int my_bounds[][2]) {
	MPI_Datatype input_filetype;
	MPI_Status status;
	MPI_Offset my_disp;
	MPI_File fp_in;
	int size, err;
	int i, j;
	char filename[256];

	long num_points = (long)n[0]*n[1]*n[2];

	//determine my starting and ending points
	//Slab decomp along y
	my_bounds[X][0] = 0;
        my_bounds[X][1] = n[X];
        my_bounds[Z][0] = 0;
        my_bounds[Z][1] = n[Z];
        my_bounds[Y][0] = my_id*n[Y]/num_procs;
        my_bounds[Y][1] = (my_id+1)*n[Y]/num_procs;

	if (my_id==num_procs-1) {
		//Pick up extras
		my_bounds[Y][1] = n[Y];
	}

	printf("%d)  x:  %d to %d, y: %d to %d, z: %d to %d\n", my_id, my_bounds[X][0], my_bounds[X][1], my_bounds[Y][0], my_bounds[Y][1], my_bounds[Z][0], my_bounds[Z][1]);
	
	int my_num_pts = (my_bounds[Y][1]-my_bounds[Y][0])*(my_bounds[Z][1]-my_bounds[Z][0])*(my_bounds[X][1]-my_bounds[X][0]);

	//MPI_Type_contiguous(my_num_pts, MPI_FLOAT, &input_filetype);
	//MPI_Type_commit(&input_filetype);

	*data = check_malloc(sizeof(float*)*3);

	char suffix[3] = {'p', 's', 'd'};

	for (i=0; i<3; i++) {
		sprintf(filename, "%s.%c", prefix, suffix[i]);
		//(*data)[i] = check_malloc(sizeof(float)*(my_bounds[X][1]-my_bounds[X][0])*(my_bounds[Y][1]-my_bounds[Y][0])*(my_bounds[Z][1]-my_bounds[Z][0]));
		(*data)[i] = check_malloc(sizeof(float)*my_num_pts);
		my_disp = sizeof(float)*(my_bounds[Y][0]*n[X]*n[Z] + my_bounds[Z][0]*n[X] + my_bounds[X][0]);
		printf("%d) reading %d points from displacement %ld.\n", my_id, my_num_pts, my_disp);

		err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_in);
		//err = MPI_File_set_view(fp_in, my_disp, MPI_FLOAT, input_filetype, "native", MPI_INFO_NULL);
		err = MPI_File_set_view(fp_in, my_disp, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
		if (err!=MPI_SUCCESS) {
	                fprintf(stderr, "MPI set view error on process %d, file %s.\n", my_id, filename);
	                char string[256];
	                int err_len;
	                MPI_Error_string(err, string, &err_len);
	                fprintf(stderr, "Error message: %s\n", string);
	                MPI_Finalize();
	                exit(1);
        	}
		if ((long)my_num_pts * (long)sizeof(float) > INT_MAX) {
			//Read in 1 GB chunks
			int ONE_GB = 1024*1024*1024;
			int ONE_GB_PTS = ONE_GB/sizeof(float);
			int num_pts_to_read = ONE_GB_PTS;
			int num_pts_read = 0;
			int data_read;
			while (num_pts_read<my_num_pts) {
				err = MPI_File_read_at_all(fp_in, 0, (*data)[i]+num_pts_read, num_pts_to_read, MPI_FLOAT, &status);
	                        MPI_Get_count(&status, MPI_INT, &data_read);
	                        printf("%d read.\n", data_read);
	                        if (err!=MPI_SUCCESS) {
	                                fprintf(stderr, "MPI error on read_at_all, process %d, file %s.\n", my_id, filename);
	                                char string[256];
	                                int err_len;
	                                MPI_Error_string(err, string, &err_len);
	                                fprintf(stderr, "Error message: %s\n", string);
	                                fprintf(stderr, "Count argument is %d\n", num_pts_to_read);
	                                MPI_Finalize();
	                                exit(1);
	                        }
				num_pts_read += num_pts_to_read;
				if (my_num_pts-num_pts_read<ONE_GB_PTS) {
					num_pts_to_read = my_num_pts-num_pts_read;
				} else {
					num_pts_to_read = ONE_GB_PTS;
				}
			}
		} else {
			//err = MPI_File_read_at_all(fp_in, 0, (*data)[i], 1, input_filetype, &status);
	                err = MPI_File_read_at_all(fp_in, 0, (*data)[i], my_num_pts, MPI_FLOAT, &status);
			int data_read;
			MPI_Get_count(&status, MPI_INT, &data_read);
			printf("%d read.\n", data_read);
	                if (err!=MPI_SUCCESS) {
	                        fprintf(stderr, "MPI error on read_at_all, process %d, file %s.\n", my_id, filename);
	                        char string[256];
	                        int err_len;
	                        MPI_Error_string(err, string, &err_len);
	                        fprintf(stderr, "Error message: %s\n", string);
				fprintf(stderr, "Count argument is %d\n", my_num_pts);
        	                MPI_Finalize();
        	                exit(1);
        	        }
		}
		err = MPI_File_close(&fp_in);
	}
	return 0;
}


int write_AWP_files(int n[], char* output_file, float** data, int my_bounds[][2]) {
	MPI_Datatype output_filetype;
	MPI_Status status;
	MPI_File fp_out;
	int i, j, k, m;
	//maximum buffer size
	float* buffer = check_malloc(sizeof(float)*3*(my_bounds[Y][1]-my_bounds[Y][0]));

	//Open output file
	printf("%d) opening file %s\n", my_id, output_file);
	int err = MPI_File_open(MPI_COMM_WORLD, output_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
	if (err!=MPI_SUCCESS) {
		fprintf(stderr, "MPI open error on process %d, file %s.\n", my_id, output_file);
		char string[256];
                int err_len;
                MPI_Error_string(err, string, &err_len);
		fprintf(stderr, "Error message: %s\n", string);
                MPI_Finalize();
                exit(1);
	}

	err = MPI_File_set_view(fp_out, 0, MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
        if (err!=MPI_SUCCESS) {
                fprintf(stderr, "MPI set view error on process %d, file %s.\n", my_id, output_file);
                char string[256];
                int err_len;
                MPI_Error_string(err, string, &err_len);
                fprintf(stderr, "Error message: %s\n", string);
                MPI_Finalize();
                exit(1);
        }

	//a contiguous chunk is ny/py elements

	//Rearrange your data in fast y, x, z, from fast x, z, y
	for (k=0; k<n[Z]; k++) {
		if (my_bounds[Z][0]<=k && my_bounds[Z][1]>k) {
			printf("%d) %d of %d z-slices.\n", my_id, k-my_bounds[Z][0], my_bounds[Z][1]-my_bounds[Z][0]);
			for (i=0; i<n[X]; i++) {
				if (my_bounds[X][0]<=i && my_bounds[X][1]>i) {
					//If I have any of these X or Z values, interleave my p, s, d data and write it
					for (j=my_bounds[Y][0]; j<my_bounds[Y][1]; j++) {
						int local_offset = (i-my_bounds[X][0])+(k-my_bounds[Z][0])*(my_bounds[X][1]-my_bounds[X][0])+(j-my_bounds[Y][0])*(my_bounds[X][1]-my_bounds[X][0])*(my_bounds[Z][1]-my_bounds[Z][0]);
						//copy Vp, Vs, rho
						buffer[3*(j-my_bounds[Y][0])] = 1000.0*data[0][local_offset];
						buffer[3*(j-my_bounds[Y][0])+1] = 1000.0*data[1][local_offset];
						buffer[3*(j-my_bounds[Y][0])+2] = 1000.0*data[2][local_offset];
						if (data[0][local_offset]*data[1][local_offset]*data[2][local_offset]==0) {
							printf("Error:  Found a 0 at local_offset=%d.\n", local_offset);
							printf("i=%d, j=%d, k=%d.\n", i, j, k);
							printf("data[0]=%f, data[1]=%f, data[2]=%f.\n", data[0][local_offset], data[1][local_offset], data[2][local_offset]);
							exit(1);
						}
					}
					//Write data
					long write_offset = 3*(my_bounds[Y][0] + ((long)i)*n[Y] + ((long)k)*n[X]*n[Y]);
					printf("i=%d, k=%d, write_offset=%ld\n", i, k, write_offset);
					MPI_File_write_at(fp_out, write_offset, buffer, 3*(my_bounds[Y][1]-my_bounds[Y][0]), MPI_FLOAT, &status);
				}
				
			}
		}
	}

	err = MPI_File_close(&fp_out);
        if (err!=MPI_SUCCESS) {
                fprintf(stderr, "MPI file close error on process %d, file %s.\n", my_id, output_file);
                char string[256];
                int err_len;
                MPI_Error_string(err, string, &err_len);
                fprintf(stderr, "Error message: %s\n", string);
                MPI_Finalize();
                exit(1);
        }

}

//Convert from fast x, z, y in 3 files to fast y, x, z in 1 file

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);

	if (argc<6) {
		printf("Usage: %s <nx> <ny> <nz> <velocity file prefix> <awp velocity out>", argv[0]);
                exit(0);
        }

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int i;
	int n[3];

	n[0] = atoi(argv[1]);
	n[1] = atoi(argv[2]);
	n[2] = atoi(argv[3]);

	char* input_prefix = argv[4];
	char* output_file = argv[5];

	float** data;
	int my_bounds[3][2];

	read_RWG_files(n, input_prefix, &data, my_bounds);
	write_AWP_files(n, output_file, data, my_bounds);
	MPI_Finalize();
}
