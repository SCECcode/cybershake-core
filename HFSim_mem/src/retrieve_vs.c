#include "include.h"

#include "ucvm.h"

int ucvm_initialized = 0;

const float MIN_VS = 500.0;

void initialize_ucvm(char* model) {
	printf("Initializing UCVM.\n");
	if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0_01302019/conf/ucvm.conf")!=UCVM_CODE_SUCCESS) {
		fprintf(stderr, "Failed to init UCVM, aborting.");
        exit(3);
    }
	char* model_string;
	if (strstr(model, ",")==NULL) {
		//No comma, only one model
		if (strcmp(model, "cvmh")==0) {
                model_string = UCVM_MODEL_CVMH;
            } else if (strcmp(model, "cvms")==0) {
                model_string = UCVM_MODEL_CVMS;
            } else if (strcmp(model, "1d")==0 || strcmp(model, "scec1d")==0) {
                model_string = UCVM_MODEL_1D;
            } else if (strcmp(model, "cvmsi")==0) {
                model_string = UCVM_MODEL_CVMSI;
            } else if (strcmp(model, "bbp1d")==0) {
                model_string = UCVM_MODEL_BBP1D;
            } else if (strcmp(model, "usgs")==0) {
                model_string = UCVM_MODEL_CENCAL;
            } else if (strcmp(model, "cca")==0) {
                model_string = "cca";
            } else {
                fprintf(stderr, "Don't recognize model %s.  Aborting.\n", model);
                exit(1);
            }
            printf("Adding model %s.\n", model_string);
            if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Error retrieving UCVM model %s.", model_string);
                exit(3);
            }

            if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
                fprintf(stderr, "Set query mode by depth failed.\n");
                exit(-2);
            }
	} else {
		//Model is a comma-separated list of models; load them all
		char* tok;
		printf("model string is :%s:\n", model);
		tok = strtok(model, ",");
		while (tok!=NULL) {
			printf("Searching for model %s.\n", tok);
	        if (strcmp(tok, "cvmh")==0) { 
	    	    model_string = UCVM_MODEL_CVMH;
	        } else if (strcmp(tok, "cvms")==0) {
	            model_string = UCVM_MODEL_CVMS;
	        } else if (strcmp(tok, "1d")==0 || strcmp(tok, "scec1d")==0) {
	            model_string = UCVM_MODEL_1D;
	        } else if (strcmp(tok, "cvmsi")==0) {
	            model_string = UCVM_MODEL_CVMSI;
	        } else if (strcmp(tok, "bbp1d")==0) {
	            model_string = UCVM_MODEL_BBP1D;
	        } else if (strcmp(tok, "usgs")==0) {
				model_string = UCVM_MODEL_CENCAL;
			} else if (strcmp(tok, "cca")==0) {
				model_string = "cca";
			} else {
				fprintf(stderr, "Don't recognize model %s.  Aborting.\n", tok);
				exit(1);
			}
	        printf("Adding model %s.\n", model_string);
	        if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
	        	fprintf(stderr, "Error retrieving UCVM model %s.", model_string);
	            exit(3);
	        }
	
	        if (ucvm_setparam(UCVM_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
	        	fprintf(stderr, "Set query mode by depth failed.\n");
	            exit(-2);
	        }
			tok = strtok(NULL, ",");
		}
	}
	printf("Initialization complete.\n");
}

//Retrieve discrete slowness-averaged Vs depth value from UCVM
//Applies Vs=500 m/s floor
float ucvm_vs_discrete(float lon, float lat, char* model, int depth, int step_size) {
	//The algorithm is:
	//VsDiscrete@depth = (depth/step_size) / (0.5/Vs(Z=0) + 1/Vs(Z=1*step_size)
	//	+ 1/Vs(Z=2*step_size) + ... + 0.5/Vs(Z=depth)
	if (!ucvm_initialized) {
		initialize_ucvm(model);
		ucvm_initialized = 1;
	}

        ucvm_point_t query_pt;
        query_pt.coord[0] = lon;
        query_pt.coord[1] = lat;
        ucvm_data_t query_data;
        int i;
        double vs_sum = 0.0;
	int steps = depth/step_size;
	for (i=0; i<=steps; i++) {
		query_pt.coord[2] = i*step_size;
                if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
                        fprintf(stderr, "UCVM query failed.\n");
                        exit(-3);
                }
		if (query_data.cmb.vs<MIN_VS) {
			query_data.cmb.vs = MIN_VS;
		}
		if (i==0 || i==steps) {
			vs_sum += 0.5/query_data.cmb.vs;
		} else {
			vs_sum += 1.0/query_data.cmb.vs;
		}
	}
	float vs = ((float)steps)/vs_sum;
	return vs;
}


//Retrieve slowness-averaged Vs depth value from UCVM, for the given model.
float ucvm_vsD(float lon, float lat, char* model, int depth) {
	//The algorithm is:
	//Vs30 = 30 / sum( 1 / (Vs sampled from [0.5, 29.5] at 1 meter increments, for 30 values) )
	if (!ucvm_initialized) {
		initialize_ucvm(model);
		ucvm_initialized = 1;
	}

	ucvm_point_t query_pt;
	query_pt.coord[0] = lon;
	query_pt.coord[1] = lat;
	ucvm_data_t query_data;
	int i;
	double vs_sum = 0.0;
	for (i=0; i<depth; i++) {
		query_pt.coord[2] = i+.5;
		if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "UCVM query failed.\n");
			exit(-3);
		}
		//printf("%d: %f\n", i, query_data.cmb.vs);
		vs_sum += 1.0/query_data.cmb.vs;
	}
	float vs = ((float)depth)/vs_sum;
	return vs;
}

float vs_at_site(float lon, float lat, char* model) {
	if (!ucvm_initialized) {
                initialize_ucvm(model);
                ucvm_initialized = 1;
        }
	        ucvm_point_t query_pt;
        query_pt.coord[0] = lon;
        query_pt.coord[1] = lat;
	query_pt.coord[2] = 0.0;
        ucvm_data_t query_data;
	if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
        	fprintf(stderr, "UCVM query failed.\n");
                exit(-3);
        }
	float vs = query_data.cmb.vs;
	if (vs<500.0) {
		vs = 500.0;
	}
	return vs;
}


//Need Vs30, VsD5H, Vs5H
//where H is the grid spacing
//Vs5H = slowness average over 5H, like Vs30
//VsD5H = discrete average over increments of H.
//So VsD500 = 5/(0.5/Vs(Z=0) + 1/Vs(Z=100) + ... + 1/Vs(Z=400) + 0.5/Vs(Z=500))
int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <lon> <lat> <model> <gridspacing (m)> <out filename>\n", argv[0]);
		return 1;
	}
	float lon = atof(argv[1]);
	float lat = atof(argv[2]);
	char* model = argv[3];
	int gridspacing = atoi(argv[4]);
	char* filename = argv[5];

	float vs30 = ucvm_vsD(lon, lat, model, 30);
	float vs5H = ucvm_vsD(lon, lat, model, 5*gridspacing);
	float vsD5H = ucvm_vs_discrete(lon, lat, model, 5*gridspacing, gridspacing);
	FILE* fp_out = fopen(filename, "w");
	fprintf(fp_out, "Vs30 = %.1f\nVs%d = %.1f\nVsD%d = %.1f\n", vs30, 5*gridspacing, vs5H, 5*gridspacing, vsD5H);
	fflush(fp_out);
	fclose(fp_out);
	return 0;
}

