#include <stdlib.h>
#include <stdio.h>

#include "ucvm.h"
#include "string.h"

int ucvm_initialized = 0;

void initialize_ucvm(char* model) {
	if (ucvm_init("/lustre/atlas/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/conf/ucvm.conf")!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Failed to init UCVM, aborting.");
                exit(3);
        }

	char* model_string;
	//If more than one model is included, load them all
	if (strstr(model, ",")==NULL) {
		if (strcmp(model, "cvmh")==0) {
			model_string = UCVM_MODEL_CVMH;
		} else if (strcmp(model, "cvms")==0) {
			model_string = UCVM_MODEL_CVMS;
		} else if (strcmp(model, "1d")==0) {
			model_string = UCVM_MODEL_1D;
		} else if (strcmp(model, "cvmsi")==0) {
			model_string = UCVM_MODEL_CVMSI;
		} else if (strcmp(model, "usgs")==0) {
			model_string = UCVM_MODEL_CENCAL;
		} else if (strcmp(model, "cvms5")==0) {
			model_string = "cvms5";
		} else if (strcmp(model, "cca1d")==0) {
			model_string = UCVM_MODEL_BBP1D;
		} else if (strcmp(model, "cca")==0) {
			model_string = "cca";
		} else {
			printf("Model string %s didn't match any known models, aborting.", model_string);
			exit(2);
		}
		printf("Adding model %s.\n", model_string);
                if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
                        fprintf(stderr, "Error retrieving UCVM model %s.", model_string);
                        exit(3);
                }
	} else {
		char* save;
	        char* tok = strtok_r(model, ",", &save);
		while (tok!=NULL) {
			if (strcmp(tok, "cvmh")==0) {
                	        model_string = UCVM_MODEL_CVMH;
                	} else if (strcmp(tok, "cvms")==0) {
                	        model_string = UCVM_MODEL_CVMS;
                	} else if (strcmp(tok, "1d")==0) {
                	        model_string = UCVM_MODEL_1D;
                	} else if (strcmp(tok, "cvmsi")==0) {
                	        model_string = UCVM_MODEL_CVMSI;
                	} else if (strcmp(tok, "usgs")==0) {
                	        model_string = UCVM_MODEL_CENCAL;
                	} else if (strcmp(tok, "cvms5")==0) {
                	        model_string = "cvms5";
                	} else if (strcmp(tok, "cca1d")==0) {
                	        model_string = UCVM_MODEL_BBP1D;
                	} else if (strcmp(tok, "cca")==0) {
                	        model_string = "cca";
                	} else {
                	        printf("Model string %s didn't match any known models, aborting.", model_string);
                	        exit(2);
                	}
                	printf("Adding model %s.\n", model_string);
                	if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
                	        fprintf(stderr, "Error retrieving UCVM model %s.", model_string);
                	        exit(3);
                	}
			tok = strtok_r(NULL, ",", &save);
		}
	}
	ucvm_initialized = 1;
}


//Retrieves vs30 using Vs30 = 30 / sum( 1 / (Vs sampled from [0, 29] at 1 meter increments, for 30 values) )
float ucvm_vs30(float lon, float lat) {
        ucvm_point_t query_pt;
        query_pt.coord[0] = lon;
        query_pt.coord[1] = lat;
        ucvm_data_t query_data;
        int i;
        double vs_sum = 0.0;
        for (i=0; i<30; i++) {
                query_pt.coord[2] = i+.5;
                if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
                        fprintf(stderr, "UCVM query failed.\n");
                        exit(-3);
                }
                //printf("%d: %f\n", i, query_data.cmb.vs);
                vs_sum += 1.0/query_data.cmb.vs;
        }
        float vs = 30.0/vs_sum;
        return vs;
}

//Determines second crossing (if it exists) of Z-values, with 10 m resolution, down to 50 km
float ucvm_zvalue(float lon, float lat, float vs_value) {
	float resolution = 100.0;
	float max_depth = 50000.0;
	int num_pts = (int)(max_depth/resolution)+1;
	ucvm_point_t* query_pts = malloc(sizeof(ucvm_point_t)*num_pts);
	ucvm_data_t* data_pts = malloc(sizeof(ucvm_data_t)*num_pts);
	int i;
	for (i=0; i<num_pts; i++) {
		query_pts[i].coord[0] = lon;
		query_pts[i].coord[1] = lat;
		query_pts[i].coord[2] = i*resolution;
	}
	if (ucvm_query(num_pts, query_pts, data_pts)!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "UCVM query failed.\n");
                exit(-3);
        }
	float depth = -1.0;
	int flag = 0;
	int crossing_num = 0;
	for (i=0; i<num_pts; i++) {
		if (crossing_num==2) {
			//We have found the 2nd crossing, stop looking
			break;
		}
		if (flag==0 && data_pts[i].cmb.vs>=vs_value) {
			depth = query_pts[i].coord[2];
			flag = 1;
			crossing_num++;
		} else if (data_pts[i].cmb.vs<vs_value) {
			flag = 0;
		}
	}

	free(query_pts);
	free(data_pts);

	return depth;
}
	

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <lat> <lon> <comma-separated model string> <output file>\n", argv[0]);
		exit(1);
	}
	float lat = atof(argv[1]);
	float lon = atof(argv[2]);
	char* model_string = argv[3];
	char* output_file = argv[4];

	initialize_ucvm(model_string);
	float vs30 = ucvm_vs30(lon, lat);
	float z10 = ucvm_zvalue(lon, lat, 1000.0);
	float z25 = ucvm_zvalue(lon, lat, 2500.0);

	FILE* fp_out = fopen(output_file, "w");
	fprintf(fp_out, "%.1f\n", vs30);
	fprintf(fp_out, "%.1f\n", z10);
	fprintf(fp_out, "%.1f\n", z25);
	fflush(fp_out);
	fclose(fp_out);
	
	return 0;
}
