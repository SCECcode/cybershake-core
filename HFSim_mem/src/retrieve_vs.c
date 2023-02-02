#include "include.h"

#include "ucvm.h"
#include "ucvm_interp.h"
#include "ucvm_model_elygtl.h"

int ucvm_initialized = 0;
int ifless_taper = 0;
double taper_depth = 700.0;
int ely_init = 0;

const float MIN_VS = 500.0;

void initialize_ucvm(char* model) {
	printf("Initializing UCVM.\n");
	if (ucvm_init("/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-22.7.0/conf/ucvm.conf")!=UCVM_CODE_SUCCESS) {
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
				//For now, turn on the ifless taper too
				ifless_taper = 1;
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
				//for now, turn on the ifless taper
				ifless_taper = 1;
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
float ucvm_vs_discrete(float lon, float lat, char* model, float surface_depth, int depth, int step_size) {
	//The algorithm is:
	//VsDiscrete@depth = (depth/step_size) / (0.5/Vs(Z=0) + 1/Vs(Z=1*step_size)
	//	+ 1/Vs(Z=2*step_size) + ... + 0.5/Vs(Z=depth)
	if (!ucvm_initialized) {
		initialize_ucvm(model);
		ucvm_initialized = 1;
	}

    int i;
    double vs_sum = 0.0;
	int steps = depth/step_size;
	ucvm_point_t* query_pts = malloc(sizeof(ucvm_point_t)*(steps+1));
	ucvm_data_t* query_data = malloc(sizeof(ucvm_data_t)*(steps+1));
	double* vs_results = malloc(sizeof(double)*(steps+1));
	for (i=0; i<=steps; i++) {
		query_pts[i].coord[0] = lon;
	    query_pts[i].coord[1] = lat;
		if (i==0) {
			query_pts[i].coord[2] = surface_depth;
		} else {
			query_pts[i].coord[2] = i*step_size;
		}
	}
	if (ucvm_query(steps+1, query_pts, query_data)!=UCVM_CODE_SUCCESS) {
	    fprintf(stderr, "UCVM query failed.\n");
        exit(-3);
    }

	for (i=0; i<=steps; i++) {
		vs_results[i] = query_data[i].cmb.vs;
	}

	if (ifless_taper!=0) {
        printf("Using ifless taper.\n");
        //If we're using the taper, we need to initialize Ely
        //Similar code to that in ucvm-single_mpi.c
        if (ely_init==0) {
            ucvm_modelconf_t conf;
            ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
            ely_init = 1;
        }
        ucvm_data_t* ely_props = malloc(sizeof(ucvm_data_t)*(steps+1));
        ucvm_point_t* ely_pts = malloc(sizeof(ucvm_point_t)*(steps+1));
        memcpy(ely_props, query_data, sizeof(ucvm_data_t)*(steps+1));
        memcpy(ely_pts, query_pts, sizeof(ucvm_point_t)*(steps+1));
        for (i=0; i<=steps; i++) {
            //Query this point, to get the crustal value at taper depth
            ely_pts[i].coord[2] = taper_depth;
        }
        if(ucvm_query(steps+1, ely_pts, ely_props)!=UCVM_CODE_SUCCESS) {
            fprintf(stderr, "UCVM query failed.");
            exit(-3);
        }
        for (i=0; i<steps+1; i++) {
            ely_props[i].domain = UCVM_DOMAIN_INTERP;
            ely_props[i].depth = query_pts[i].coord[2];
        }
        ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, UCVM_COORD_GEO_DEPTH, steps+1, ely_pts, ely_props);
        for (i=0; i<=steps; i++) {
            ucvm_interp_ely(0.0, taper_depth, UCVM_COORD_GEO_DEPTH, ely_pts+i, ely_props+i);
        }
        //Combine results
        for (i=0; i<steps+1; i++) {
            //printf("Ely prop value at depth %f = %f\n", ely_props[i].depth, ely_props[i].cmb.vs);
            if (ely_props[i].cmb.vs<vs_results[i]) {
                vs_results[i] = ely_props[i].cmb.vs;
            }
        }
        free(ely_props);
        free(ely_pts);
    }

	for (i=0; i<=steps; i++) {
		if (vs_results[i]<MIN_VS) {
			vs_results[i] = MIN_VS;
		}
		printf("%f: %f\n", query_pts[i].coord[2], vs_results[i]);
		if (i==0 || i==steps) {
			vs_sum += 0.5/vs_results[i];
		} else {
			vs_sum += 1.0/vs_results[i];
		}
	}
	float vs = ((float)steps)/vs_sum;
	free(query_pts);
	free(query_data);
	free(vs_results);
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

	//Query without taper	
	ucvm_point_t* query_pts = malloc(sizeof(ucvm_point_t)*depth);
	ucvm_data_t* query_data = malloc(sizeof(ucvm_data_t)*depth);
	double* vs_results = malloc(sizeof(double)*depth);
	int i;
	for (i=0; i<depth; i++) {
		query_pts[i].coord[0] = lon;
	    query_pts[i].coord[1] = lat;
		query_pts[i].coord[2] = i+.5;
	}
	if (ucvm_query(depth, query_pts, query_data)!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "UCVM query failed.\n");
			exit(-3);
	}
	for (i=0; i<depth; i++) {
		vs_results[i] = query_data[i].cmb.vs;
		printf("%f: %f\n", query_pts[i].coord[2], query_data[i].cmb.vs);
	}
	//If we need to include the taper
	if (ifless_taper!=0) {
		printf("Using ifless taper.\n");
	    //If we're using the taper, we need to initialize Ely
    	//Similar code to that in ucvm-single_mpi.c
    	if (ely_init==0) {
    	    ucvm_modelconf_t conf;
    	    ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
    	    ely_init = 1;
    	}
		ucvm_data_t* ely_props = malloc(sizeof(ucvm_data_t)*depth);
		ucvm_point_t* ely_pts = malloc(sizeof(ucvm_point_t)*depth);
		memcpy(ely_props, query_data, sizeof(ucvm_data_t)*depth);
		memcpy(ely_pts, query_pts, sizeof(ucvm_point_t)*depth);
		for (i=0; i<depth; i++) {
			//Query this point, to get the crustal value at taper depth
			ely_pts[i].coord[2] = taper_depth;
		}
		if(ucvm_query(depth, ely_pts, ely_props)!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "UCVM query failed.");
			exit(-3);
		}
		for (i=0; i<depth; i++) {
			ely_props[i].domain = UCVM_DOMAIN_INTERP;
			ely_props[i].depth = query_pts[i].coord[2];
		}
		ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, UCVM_COORD_GEO_DEPTH, depth, ely_pts, ely_props);
		for (i=0; i<depth; i++) {
			ucvm_interp_ely(0.0, taper_depth, UCVM_COORD_GEO_DEPTH, ely_pts+i, ely_props+i);
		}
		//Combine results
		for (i=0; i<depth; i++) {
			//printf("Ely prop value at depth %f = %f\n", ely_props[i].depth, ely_props[i].cmb.vs);
			if (ely_props[i].cmb.vs<vs_results[i]) {
				vs_results[i] = ely_props[i].cmb.vs;
			}
		}
		free(ely_props);
		free(ely_pts);
	}
	
	//Slowness average
	double vs_sum = 0.0;
	for (i=0; i<depth; i++) {
		vs_sum += 1.0/vs_results[i];
	}
	float vs = ((float)depth)/vs_sum;
	free(vs_results);
	free(query_data);
	free(query_pts);
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
	if (argc<7) {
		printf("Usage: %s <lon> <lat> <model> <gridspacing (m)> <depth to query for surface mesh point (m)> <out filename>\n", argv[0]);
		return 1;
	}
	float lon = atof(argv[1]);
	float lat = atof(argv[2]);
	char* model = argv[3];
	int gridspacing = atoi(argv[4]);
	float surface_depth = atof(argv[5]);
	char* filename = argv[6];

	float vs30 = ucvm_vsD(lon, lat, model, 30);
	float vs5H = ucvm_vsD(lon, lat, model, 5*gridspacing);
	float vsD5H = ucvm_vs_discrete(lon, lat, model, surface_depth, 5*gridspacing, gridspacing);
	FILE* fp_out = fopen(filename, "w");
	fprintf(fp_out, "Vs30 = %.1f\nVs%d = %.1f\nVsD%d = %.1f\n", vs30, 5*gridspacing, vs5H, 5*gridspacing, vsD5H);
	fflush(fp_out);
	fclose(fp_out);
	return 0;
}

