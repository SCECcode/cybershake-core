/* Parallel code to create a file describing the velocity model used at each point, for use with smoothing */

#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "mpi.h"

#include "ucvm.h"

int my_id;
int num_procs;

void create_surface(ucvm_point_t* points, char* model_coords, int starting_strip, int ending_strip, int nx, int num_points) {
	FILE* fp_in = fopen(model_coords, "r");
	int i, j, x, y;
	j = 0;
	float lon, lat;
	//Fast X
	printf("Reading %d points\n", num_points);
	for (i=0; i<num_points; i++) {
		fscanf(fp_in, " %f %f %d %d\n", &lon, &lat, &x, &y);
		if (y>=starting_strip && y<ending_strip) {
			points[j].coord[0] = lon;
			points[j].coord[1] = lat;
			points[j].coord[2] = 10000.0;
			j++;
		}
		if (y>ending_strip) {
			//We're done parsing the file
			break;
		}
	}
	fclose(fp_in);
	//Check that we've read in the expected # of points
	if (j!=(ending_strip-starting_strip)*nx) {
		printf("Expected to read in %d points (starting strip=%d, ending strip=%d, nx=%d), but read in %d points.  Aborting.\n", (ending_strip-starting_strip)*nx, starting_strip, ending_strip, nx, j);
		MPI_Finalize();
		exit(1);
	}
}

void query_model(ucvm_point_t* points, ucvm_data_t** data, char* velocity_models, int num_local_points) {
	if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
		fprintf(stderr, "Failed to setup ucvm.\n");
                exit(-1);
        }
	printf("vel model string = %s\n", velocity_models);
	if (strstr(velocity_models, ",")!=NULL) {
		char* save;
		char* tok = strtok_r(velocity_models, ",", &save);
		while (tok!=NULL) {
			printf("tok=%s\n", tok);
			if (strcmp(tok, "cvms5")==0) {
				printf("Adding cvms5.\n");
				if (ucvm_add_model("cvms5")!=UCVM_CODE_SUCCESS) {
                        		fprintf(stderr, "Error retrieving CVM-S5.\n");
                        		fflush(stderr);
                        		exit(-2);
                        	}
			} else if (strcmp(tok, "usgs")==0) {
				printf("Adding cencal.\n");
                                if (ucvm_add_model(UCVM_MODEL_CENCAL)!=UCVM_CODE_SUCCESS) {
                                        fprintf(stderr, "Error retrieving USGS Bay Area.\n");
                                        fflush(stderr);
                                        exit(-2);
                                }
			} else if (strcmp(tok, "1d")==0) {
				printf("Adding 1D.\n");
                                if (ucvm_add_model(UCVM_MODEL_1D)!=UCVM_CODE_SUCCESS) {
                                        fprintf(stderr, "Error retrieving 1D model.\n");
                                        fflush(stderr);
                                        exit(-2);
                                }
			} else if (strcmp(tok, "cca")==0) {
                                printf("Adding CCA.\n");
                                if (ucvm_add_model("cca")!=UCVM_CODE_SUCCESS) {
                                        fprintf(stderr, "Error retrieving CCA model.\n");
                                        fflush(stderr);
                                        exit(-2);
                                }
                        } else if (strcmp(tok, "cvmsi")==0) {
                                printf("Adding CVM-SI.\n");
                                if (ucvm_add_model(UCVM_MODEL_CVMSI)!=UCVM_CODE_SUCCESS) {
                                        fprintf(stderr, "Error retrieving CVMSI model.\n");
                                        fflush(stderr);
                                        exit(-2);
                                }
			}
			tok = strtok_r(NULL, ",", &save);
		}
	} else {
		//We just have 1 model
		if (strcmp(velocity_models, "cvms5")==0) {
                       printf("Adding cvms5.\n");
                       if (ucvm_add_model("cvms5")!=UCVM_CODE_SUCCESS) {
                       		fprintf(stderr, "Error retrieving CVM-S5.\n");
                        	fflush(stderr);
                        	exit(-2);
                	}
                } else if (strcmp(velocity_models, "usgs")==0) {
                        printf("Adding cencal.\n");
                        if (ucvm_add_model(UCVM_MODEL_CENCAL)!=UCVM_CODE_SUCCESS) {
                	        fprintf(stderr, "Error retrieving USGS Bay Area.\n");
                                fflush(stderr);
                                exit(-2);
                        }
                } else if (strcmp(velocity_models, "1d")==0) {
                	printf("Adding 1D.\n");
                        if (ucvm_add_model(UCVM_MODEL_1D)!=UCVM_CODE_SUCCESS) {
                        	fprintf(stderr, "Error retrieving 1D model.\n");
                                fflush(stderr);
                                exit(-2);
                        }
                } else if (strcmp(velocity_models, "cca")==0) {
                        printf("Adding CCA.\n");
                        if (ucvm_add_model("cca")!=UCVM_CODE_SUCCESS) {
                	        fprintf(stderr, "Error retrieving CCA model.\n");
                                fflush(stderr);
                                exit(-2);
                        }
                } else if (strcmp(velocity_models, "cvmsi")==0) {
                        printf("Adding CVM-SI.\n");
                        if (ucvm_add_model(UCVM_MODEL_CVMSI)!=UCVM_CODE_SUCCESS) {
                 	       fprintf(stderr, "Error retrieving CVMSI model.\n");
                               fflush(stderr);
                               exit(-2);
                        }
                } else if (strcmp(velocity_models, "cca1d")==0) {
			printf("Adding CCA 1D.\n");
                        if (ucvm_add_model(UCVM_MODEL_BBP1D) != UCVM_CODE_SUCCESS) {
                           fprintf(stderr, "Error retrieving CCA 1D model.\n");
                           fflush(stderr);
			   exit(-2);
                         }
		}

	}
        /*char label_test[64];
        ucvm_model_label(1, label_test, 64);
        printf("ID 1 goes with label %s.\n", label_test);
        ucvm_model_label(-1, label_test, 64);
        printf("ID -1 goes with label %s.\n", label_test);
        ucvm_model_label(-2, label_test, 64);
        printf("ID -2 goes with label %s.\n", label_test);*/

	printf("%d) Point 0 = (%f, %f, %f)\n", my_id, points[0].coord[0], points[0].coord[1], points[0].coord[2]);

	printf("Querying %d points.\n", num_local_points);
	if (ucvm_query(num_local_points, points, *data)!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Query UCVM failed.\n");
                exit(-3);
	}	
}

void write_models(ucvm_data_t* data, char* output_filename, int nx, int ny, int num_local_pts) {
	//Writing a text file, so everyone assembles their IDs and sends that to the master
	//Use IDs instead of labels, then we only have to send ints
	printf("%d) Preparing to write results.\n", my_id);
	int i,j;
	int* ids = malloc(sizeof(int)*num_local_pts);
	int* all_ids = NULL;
	int* recvcounts = NULL;
	int* displs = NULL;
	for (i=0; i<num_local_pts; i++) {
		ids[i] = data[i].crust.source;
	}
	if (my_id==0) {
		all_ids = malloc(sizeof(int)*nx*ny);
		//Determine recvcounts, displs
		recvcounts = malloc(sizeof(int)*num_procs);
		displs = malloc(sizeof(int)*num_procs);
		float avg_strips_per_proc = ((float)ny)/((float)num_procs);
		for (i=0; i<num_procs; i++) {
		        int starting_strip = (int)(i*avg_strips_per_proc);
		        int ending_strip = (int)((i+1)*avg_strips_per_proc);
		        if (ending_strip>ny) {
		                ending_strip = ny;
		        }
			recvcounts[i] = nx*(ending_strip-starting_strip);
			if (i==0) {
				displs[i] = 0;
			} else {
				displs[i] = displs[i-1] + recvcounts[i-1];
			}
		}
	}
	printf("%d) Entering gatherv\n", my_id);
	MPI_Gatherv(ids, num_local_pts, MPI_INT, all_ids, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	printf("%d) Finished gatherv\n", my_id);
	if (my_id==0) {
		FILE* fp_out = fopen(output_filename, "w");
		char id_to_model[100][64];
		int i, j;
		for (i=0; i<100; i++) {
			id_to_model[i][0] = '\0';
		}
		for (i=0; i<ny; i++) {
			for (j=0; j<nx; j++) {
				int id = all_ids[i*nx+j];
				//int id = data[i*nx+j].crust.source;
				if (id>=0) {
					if (id_to_model[id][0]=='\0') {
						ucvm_model_label(id, id_to_model[id], 64);
						printf("ID=%d, label=%s\n", id, id_to_model[id]);
					}
					fprintf(fp_out, "%s ", id_to_model[id]);
				} else {
					fprintf(fp_out, "%d ", id);
				}
			}
			fprintf(fp_out, "\n");
		}
		fflush(fp_out);
		fclose(fp_out);
	}
	free(all_ids);
	free(recvcounts);
	free(displs);
	free(ids);
}

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	if (argc<5) {
		printf("Usage: %s <gridout file> <model coords file> <comma-separated list of velocity models> <output file>", argv[0]);
		exit(1);
	}

	char* gridout_file = argv[1];
	char* model_coords = argv[2];
	char* velocity_models = argv[3];
	char* output_filename = argv[4];

	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int i;
	int dims[2];
	if (my_id==0) {
		//Determine nx, ny
		FILE* fp_in = fopen(gridout_file, "r");
		char line[256];
		//xlen=440.00000
		//nx=1100
		fgets(line, 256, fp_in);
		fscanf(fp_in, "nx=%d\n", &dims[0]);
		for (i=0; i<dims[0]; i++) {
			fgets(line, 256, fp_in);
		}
		//ylen
		fgets(line, 256, fp_in);
		printf("%s", line);
		fscanf(fp_in, "ny=%d", &dims[1]);
		fclose(fp_in);
		printf("nx=%d, ny=%d\n", dims[0], dims[1]);
	}
	//Broadcast size to everyone
	MPI_Bcast(dims, 2, MPI_INT, 0, MPI_COMM_WORLD);
	int nx = dims[0];
	int ny = dims[1];

	//Divide surface using 1-D decomp along Y, so everyone gets X-strips
	if (num_procs>ny) {
		printf("Trying to use too many processors.  There are only %d work units, so you can use no more than %d processors.\n", ny, ny);
		MPI_Finalize();
		exit(1);
	}
	float avg_strips_per_proc = ((float)ny)/((float)num_procs);
	int starting_strip = (int)(my_id*avg_strips_per_proc);
	int ending_strip = (int)((my_id+1)*avg_strips_per_proc);
	if (ending_strip>ny) {
		ending_strip = ny;
	}
	printf("%d) Responsible for strips %d - %d\n", my_id, starting_strip, ending_strip);
	int num_strips = ending_strip - starting_strip;
	int num_local_pts = num_strips*nx;

	ucvm_point_t* points = malloc(sizeof(ucvm_point_t)*num_local_pts);
	ucvm_data_t* data = malloc(sizeof(ucvm_data_t)*num_local_pts);

	create_surface(points, model_coords, starting_strip, ending_strip, nx, nx*ny);
	query_model(points, &data, velocity_models, num_local_pts);
        /*int j;
        for (i=1500; i<1510; i++) {
                for (j=500; j<510; j++) {
                        int index = i*nx + j;
                        printf("Queried point (%lf, %lf, 10km), ", points[index].coord[0], points[index].coord[1]);
                        printf("Returned vp=%lf, vs=%lf, rho=%lf from model ids %d, %d, %d\n", data[index].cmb.vp, data[index].cmb.vs, data[index].cmb.rho, data[index].crust.source, data[index].gtl.source, data[index].cmb.source);
                }
        }*/
	/*char label_test[64];
	ucvm_model_label(1, label_test, 64);
	printf("ID 1 goes with label %s.\n", label_test);
	ucvm_model_label(-1, label_test, 64);
	printf("ID -1 goes with label %s.\n", label_test);
        ucvm_model_label(-2, label_test, 64);
        printf("ID -2 goes with label %s.\n", label_test);*/

	write_models(data, output_filename, nx, ny, num_local_pts);

	free(points);
	free(data);
	MPI_Finalize();
}
