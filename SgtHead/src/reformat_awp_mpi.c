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
	fprintf(stderr, "%d) Initial file offset = %ld\n", my_id, disp);
	err = MPI_File_set_view(fp_in, disp, MPI_BYTE, sgt_pt_datatype, "native", MPI_INFO_NULL);
	//err = MPI_File_seek(fp_in, disp, MPI_SEEK_SET);

	//Set up output file
	err = MPI_File_open(MPI_COMM_WORLD, output_filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fp_out);
	err = MPI_File_set_view(fp_out, disp, MPI_BYTE, sgt_pt_datatype, "native", MPI_INFO_NULL);
	//err = MPI_File_seek(fp_out, disp, MPI_SEEK_SET);

	//Restrict reads to 7.5 GB
	float MAX_READ = 2.0*1024*1024*1024;
	int num_reads = ceil((((float)(my_num_pts))*sgt_point_size)/MAX_READ);
	fprintf(stderr, "%d) %f / %f\n", my_id, (((float)(my_num_pts))*sgt_point_size), MAX_READ);
	int pts_per_read = ceil(((float)my_num_pts)/((float)num_reads));
	fprintf(stderr, "%d) %d reads, %d pts per read.\n", my_id, num_reads, pts_per_read);
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
	float* single_pt_buffer = malloc(sgt_point_size);
	for (i=0; i<num_reads; i++) {
        //for(i=0; i<1; i++) {
		printf("%d) Read %d of %d.\n", my_id, i, num_reads);
		if (i==num_reads-1) {
			//Set up last read
			pts_in_read = my_num_pts - (num_reads-1)*pts_per_read;
		}
		fprintf(stderr, "%d) Reading from point %d to point %d (%d points, %ld bytes).\n", my_id, i*pts_per_read, i*pts_per_read + pts_in_read, pts_in_read, ((long)pts_in_read)*sgt_point_size);
		err = MPI_File_read(fp_in, buffer, pts_in_read, sgt_pt_datatype, &status);

		if (z_comp==1) {
			//double everything and negate for flipped source
			for (j=0; j<6*timesteps*pts_in_read; j++) {
				buffer[j] = -2.0*buffer[j];
			}
		}
		for (j=0; j<pts_in_read; j++) {
		//for (j=0; j<2; j++) {
			/*if (my_id==0) {
				for (k=0; k<timesteps; k++) {
					fprintf(stderr, "%d) buffer_XX[%d] = %e\n", my_id, k, buffer[j*6*timesteps+k]);
				}
			}*/
			//YY -> XX
			memcpy(single_pt_buffer+XX*timesteps, buffer+j*6*timesteps+YY*timesteps, timesteps*sizeof(float));
			//XX -> YY
			memcpy(single_pt_buffer+YY*timesteps, buffer+j*6*timesteps+XX*timesteps, timesteps*sizeof(float));
			//ZZ -> ZZ
			memcpy(single_pt_buffer+ZZ*timesteps, buffer+j*6*timesteps+ZZ*timesteps, timesteps*sizeof(float));
			//XY -> XY
			memcpy(single_pt_buffer+XY*timesteps, buffer+j*6*timesteps+XY*timesteps, timesteps*sizeof(float));
			//apply the -1 to XZ and YZ
                	for (k=0; k<timesteps; k++) {
                	        buffer[j*6*timesteps + XZ*timesteps+k] = -1.0*buffer[j*6*timesteps + XZ*timesteps+k];
	                        buffer[j*6*timesteps + YZ*timesteps+k] = -1.0*buffer[j*6*timesteps + YZ*timesteps+k];
                	}
			//-YZ -> XZ
			memcpy(single_pt_buffer+XZ*timesteps, buffer+j*6*timesteps+YZ*timesteps, timesteps*sizeof(float));
			//-XZ -> YZ
			memcpy(single_pt_buffer+YZ*timesteps, buffer+j*6*timesteps+XZ*timesteps, timesteps*sizeof(float));
			//copy back into buffer
			memcpy(buffer+j*6*timesteps, single_pt_buffer, sgt_point_size);
			/*if (my_id==0) {
				for (k=0; k<timesteps; k++) {
					fprintf(stderr, "%d) single_pt_buffer_XX[%d] = %e\n", my_id, k, single_pt_buffer[k]);
				}
			}*/
		}
		//Write back to output file
		if (((long)pts_in_read)*sgt_point_size > (long)2*1024*2024*1024) {
			//Break read into pieces
			int num_writes = ceil(MAX_READ/(2.0*1024.0*1024.0*1024.0));
			int num_pts_per_write = ceil(((float)pts_in_read)/((float)num_writes));
			fprintf(stderr, "%d) Preparing for %d writes.\n", my_id, num_writes);
			for (j=0; j<num_writes-1; j++) {
				MPI_Offset offset;
                        	MPI_File_get_position(fp_out, &offset);
                        	fprintf(stderr, "%d) Writing %d points at location %ld.\n", my_id, num_pts_per_write, offset);
				err = MPI_File_write(fp_out, buffer + j*num_pts_per_write*6*timesteps, num_pts_per_write, sgt_pt_datatype, &status);
				if (err!=MPI_SUCCESS) {
                 		       char error_string[1024];
                        		int len, e_class;
                        		MPI_Error_class(err, &e_class);
                        		MPI_Error_string(e_class, error_string, &len);
                        		fprintf(stderr, "%d) %s\n", my_id, error_string);
                        		MPI_Error_string(err, error_string, &len);
                        		fprintf(stderr, "%d) %s\n", my_id, error_string);
				}
               		}
			MPI_Offset offset;
                        MPI_File_get_position(fp_out, &offset);
                        fprintf(stderr, "%d) Writing %d points at location %ld.\n", my_id, pts_in_read - (num_writes-1)*num_pts_per_write, offset);
			err = MPI_File_write(fp_out, buffer + (num_writes-1)*num_pts_per_write, pts_in_read - (num_writes-1)*num_pts_per_write, sgt_pt_datatype, &status);
			if (err!=MPI_SUCCESS) {
                        	char error_string[1024];
                                int len, e_class;
                                MPI_Error_class(err, &e_class);
                                MPI_Error_string(e_class, error_string, &len);
                                fprintf(stderr, "%d) %s\n", my_id, error_string);
                                MPI_Error_string(err, error_string, &len);
                                fprintf(stderr, "%d) %s\n", my_id, error_string);
                	}
		} else {
			MPI_Offset offset;
			MPI_File_get_position(fp_out, &offset);
			fprintf(stderr, "%d) Writing %d points at location %ld.\n", my_id, pts_in_read, offset);
			err = MPI_File_write(fp_out, buffer, pts_in_read, sgt_pt_datatype, &status);
		}
		if (err!=MPI_SUCCESS) {
			char error_string[1024];
			int len, e_class;
			MPI_Error_class(err, &e_class);
			MPI_Error_string(e_class, error_string, &len);
			fprintf(stderr, "%d) %s\n", my_id, error_string);
			MPI_Error_string(err, error_string, &len);
			fprintf(stderr, "%d) %s\n", my_id, error_string);
		}
		//Set up next read
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_File_close(&fp_in);
	MPI_File_close(&fp_out);
	return 0;
}
