#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

//Alters awp:  swaps XX and YY, XZ and YZ, multiplies XZ and YZ by -1
//Insert in order (RWG): XX YY ZZ XY XZ YZ

int my_id;
int num_procs;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	if (argc<5) {
		printf("Usage: %s <input AWP SGT> <timesteps> <num_sgt_points> <output AWP SGT> <-z>", argv[0]);
		exit(1);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	char* input_filename = argv[1];
        int timesteps = atoi(argv[2]);
	int num_sgt_points = atoi(argv[3]);
	char* output_filename = argv[4];

        int z_comp = 0;
        if (argc==6 && strcmp(argv[5], "-z")==0) {
                z_comp = 1;
                printf("Will double all SGTs since they're the Z source\n");
        }

        long sgt_point_size = timesteps*6*sizeof(float);
        int my_start_pt, my_end_pt, my_num_pts;
        my_start_pt = num_sgt_points/num_procs * my_id;
        my_end_pt = num_sgt_points/num_procs * (my_id+1);
        if (my_id==num_procs-1) {
                my_end_pt = num_sgt_points;
        }
        my_num_pts = my_end_pt - my_start_pt;

	//Everyone opens file
	MPI_File fp_in, fp_out;
	int err = MPI_File_open(MPI_COMM_WORLD, input_filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_in);

	MPI_Datatype sgt_pt_datatype;
        MPI_Type_contiguous(sgt_point_size, MPI_BYTE, &sgt_pt_datatype);
        MPI_Type_commit(&sgt_pt_datatype);

	MPI_Offset disp = my_start_pt * sgt_point_size;
	err = MPI_File_set_view(fp_in, disp, MPI_BYTE, sgt_pt_datatype, "native", MPI_INFO_NULL);

	//Set up output file
	err = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_out);
	err = MPI_File_set_view(fp_out, disp, MPI_BYTE, sgt_pt_datatype, "native", MPI_INFO_NULL);

	//Restrict reads to 1.5 GB
	long MAX_READ = (long)1024*1024*1024/2*3;
	int num_reads = ceil(((long)(my_num_pts)*sgt_point_size)/MAX_READ);
	int pts_per_read = ceil(my_num_pts/num_reads);
	float* buffer = malloc(sgt_point_size * pts_per_read);
	int pts_in_read = pts_per_read;
	int i, k;
	long j;

        int XX = 0;
        int YY = 1;
        int ZZ = 2;
        int XY = 3;
        int XZ = 4;
        int YZ = 5;

	MPI_Status status;
	for (i=0; i<num_reads; i++) {
		printf("%d) Read %d of %d.\n", my_id, i, num_reads);
		if (i==num_reads-1) {
			//Set up last read
			pts_in_read = my_num_pts - (num_reads-1)*pts_per_read;
		}
		err = MPI_File_read(fp_in, buffer, pts_in_read, sgt_pt_datatype, &status);
		float* single_pt_buffer = malloc(sgt_point_size);

		if (z_comp==1) {
			//double everything and negate for flipped source
			for (j=0; j<6*timesteps*pts_in_read; j++) {
				buffer[j] = -2.0*buffer[j];
			}
		}
		for (j=0; j<pts_in_read; j++) {
			//YY -> XX
			memcpy(single_pt_buffer+YY*timesteps, buffer+j*sgt_point_size, timesteps*sizeof(float));
			//XX -> YY
			memcpy(single_pt_buffer+XX*timesteps, buffer+j*sgt_point_size+timesteps, timesteps*sizeof(float));
			//ZZ -> ZZ
			memcpy(single_pt_buffer+ZZ*timesteps, buffer+j*sgt_point_size+2*timesteps, timesteps*sizeof(float));
			//XY -> XY
			memcpy(single_pt_buffer+XY*timesteps, buffer+j*sgt_point_size+3*timesteps, timesteps*sizeof(float));
			//apply the -1 to XZ and YZ
                	for (k=0; k<timesteps; k++) {
                	        buffer[j*sgt_point_size + XZ*timesteps+k] = -1.0*buffer[j*sgt_point_size + XZ*timesteps+k];
	                        buffer[j*sgt_point_size + YZ*timesteps+k] = -1.0*buffer[j*sgt_point_size + YZ*timesteps+k];
                	}
			//-YZ -> XZ
			memcpy(single_pt_buffer+YZ*timesteps, buffer+j*sgt_point_size+4*timesteps, timesteps*sizeof(float));
			//-XZ -> YZ
			memcpy(single_pt_buffer+XZ*timesteps, buffer+j*sgt_point_size+5*timesteps, timesteps*sizeof(float));
			//copy back into buffer
			memcpy(buffer+j*sgt_point_size, single_pt_buffer, sgt_point_size);
		}
		//Write back to output file
		err = MPI_File_write(fp_out, buffer, pts_in_read, sgt_pt_datatype, &status);
		//Set up next read
	}

	MPI_File_close(&fp_in);
	MPI_File_close(&fp_out);
	return 0;
}
