/* Code to create a file describing the velocity model used at each point, for use with smoothing */

#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "ucvm.h"


void create_surface(ucvm_point_t* points, char* model_coords, int num_points) {
	FILE* fp_in = fopen(model_coords, "r");
	int i, x, y;
	//Fast X
	printf("Reading %d points\n", num_points);
	for (i=0; i<num_points; i++) {
		fscanf(fp_in, " %lf %lf %d %d\n", &(points[i].coord[0]), &(points[i].coord[1]), &x, &y);
		points[i].coord[2] = 10000.0;
	}
	fclose(fp_in);
}

void query_model(ucvm_point_t* points, ucvm_data_t** data, char* velocity_models, int num_points) {
	if (ucvm_init("/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/conf/ucvm.conf") != UCVM_CODE_SUCCESS) {
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
		if (strcmp(velocity_model, "cvms5")==0) {
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
        /*char label_test[64];
        ucvm_model_label(1, label_test, 64);
        printf("ID 1 goes with label %s.\n", label_test);
        ucvm_model_label(-1, label_test, 64);
        printf("ID -1 goes with label %s.\n", label_test);
        ucvm_model_label(-2, label_test, 64);
        printf("ID -2 goes with label %s.\n", label_test);*/

	printf("Querying %d points.\n", num_points);
	if (ucvm_query(num_points, points, *data)!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Query UCVM failed.\n");
                exit(-3);
	}	
}

void write_models(ucvm_data_t* data, char* output_filename, int nx, int ny) {
	FILE* fp_out = fopen(output_filename, "w");
	char id_to_model[100][64];
	int i, j;
	for (i=0; i<100; i++) {
		id_to_model[i][0] = '\0';
	}
	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) {
			int id = data[i*nx+j].crust.source;
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

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <gridout file> <model coords file> <comma-separated list of velocity models> <output file>", argv[0]);
		exit(1);
	}

	char* gridout_file = argv[1];
	char* model_coords = argv[2];
	char* velocity_models = argv[3];
	char* output_filename = argv[4];

	int i;
	//Determine nx, ny
	FILE* fp_in = fopen(gridout_file, "r");
	char line[256];
	int nx, ny;
	//xlen=440.00000
	//nx=1100
	fgets(line, 256, fp_in);
	fscanf(fp_in, "nx=%d\n", &nx);
	for (i=0; i<nx; i++) {
		fgets(line, 256, fp_in);
	}
	//ylen
	fgets(line, 256, fp_in);
	printf("%s", line);
	fscanf(fp_in, "ny=%d", &ny);
	fclose(fp_in);

	printf("nx=%d, ny=%d\n", nx, ny);

	ucvm_point_t* points = malloc(sizeof(ucvm_point_t)*nx*ny);
	ucvm_data_t* data = malloc(sizeof(ucvm_data_t)*nx*ny);

	create_surface(points, model_coords, nx*ny);
	query_model(points, &data, velocity_models, nx*ny);
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

	write_models(data, output_filename, nx, ny);

	free(points);
	free(data);
}
