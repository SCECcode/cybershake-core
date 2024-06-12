#include <stdlib.h>
#include <stdio.h>

#include "ucvm.h"
#include "ucvm_interp.h"
#include "ucvm_model_elygtl.h"

#include "string.h"

int ucvm_initialized = 0;
int ifless_taper = 0;
double taper_depth = 700.0;
int ely_init = 0;

const char* UCVM_CONF = "/lustre/orion/geo156/proj-shared/CyberShake/software/UCVM/ucvm_22.7.0_withSFCVM/conf/ucvm.conf";

void initialize_ucvm(char* model) {
	if (ucvm_initialized) {
		return;
	}
	if (ucvm_init(UCVM_CONF)!=UCVM_CODE_SUCCESS) {
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
		} else if (strcmp(model, "bbp1d")==0) {
			model_string = "bbp1d";
		} else if (strcmp(model, "sfcvm")==0) {
			model_string = "sfcvm";
		} else {
			printf("Model string %s didn't match any known models, aborting.", model_string);
			exit(2);
		}
		printf("Adding model %s.\n", model_string);
                if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
                        fprintf(stderr, "Error retrieving UCVM model %s.\n", model_string);
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
					} else if (strcmp(tok, "sfcvm")==0) {
							model_string = "sfcvm";
					} else if (strcmp(tok, "bbp1d")==0) {
							model_string = UCVM_MODEL_BBP1D;
                	} else {
                	        printf("Model string %s didn't match any known models, aborting.\n", tok);
                	        exit(2);
                	}
                	printf("Adding model %s.\n", model_string);
                	if (ucvm_add_model(model_string)!=UCVM_CODE_SUCCESS) {
                	        fprintf(stderr, "Error retrieving UCVM model %s.\n", model_string);
                	        exit(3);
                	}
			tok = strtok_r(NULL, ",", &save);
		}
	}
	ucvm_initialized = 1;
}

int test_for_taper(float lon, float lat, char* taper_models) {
	ucvm_point_t query_pt;
	ucvm_data_t query_data;
	query_pt.coord[0] = lon;
	query_pt.coord[1] = lat;
	query_pt.coord[2] = 0.0;
	query_pt.coord[2] = 0.0;
	if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
		fprintf(stderr, "UCVM query failed.\n");
        exit(-3);
    }
	char pt_label[64];
	ucvm_model_label(query_data.crust.source, pt_label, 64);
	char* tok = strtok(taper_models, ",");
	char tmp[64];
	while (tok!=NULL) {
		strcpy(tmp, tok);
		if (strcmp(tmp, pt_label)==0) {
			//Yes, we will be needing the taper
			return 1;
		}
		tok = strtok(NULL, ",");
	}
	return 0;
}


//Retrieves vs30 using Vs30 = 30 / sum( 1 / (Vs sampled from [0.5, 29.5] at 1 meter increments, for 30 values) )
float ucvm_vs30(float lon, float lat, int taper_flag) {
		ucvm_point_t* query_pts = malloc(sizeof(ucvm_point_t)*30);
		ucvm_data_t* query_data = malloc(sizeof(ucvm_data_t)*30);
		double* vs_results = malloc(sizeof(double)*30);
		int i;
		for (i=0; i<30; i++) {
	        query_pts[i].coord[0] = lon;
	        query_pts[i].coord[1] = lat;
			query_pts[i].coord[2] = i+0.5;
		}
        if (ucvm_query(30, query_pts, query_data)!=UCVM_CODE_SUCCESS) {
        	fprintf(stderr, "UCVM query failed.\n");
            exit(-3);
        }
		for (i=0; i<30; i++) {
			vs_results[i] = query_data[i].cmb.vs;
            //printf("%d: %f\n", i, query_data[i].cmb.vs);
        }
		
		//Include the taper, if needed
	    if (taper_flag!=0) {
	        printf("Using ifless taper.\n");
    	    //If we're using the taper, we need to initialize Ely
    	    //Similar code to that in ucvm-single_mpi.c
    	    if (ely_init==0) {
    	        ucvm_modelconf_t conf;
    	        ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
    	        ely_init = 1;
    	    }
    	    ucvm_data_t* ely_props = malloc(sizeof(ucvm_data_t)*30);
    	    ucvm_point_t* ely_pts = malloc(sizeof(ucvm_point_t)*30);
    	    memcpy(ely_props, query_data, sizeof(ucvm_data_t)*30);
    	    memcpy(ely_pts, query_pts, sizeof(ucvm_point_t)*30);
    	    for (i=0; i<30; i++) {
    	        //Query this point, to get the crustal value at taper depth
    	        ely_pts[i].coord[2] = taper_depth;
	        }
	        if(ucvm_query(30, ely_pts, ely_props)!=UCVM_CODE_SUCCESS) {
	            fprintf(stderr, "UCVM query failed.");
	            exit(-3);
	        }
    	    for (i=0; i<30; i++) {
	            ely_props[i].domain = UCVM_DOMAIN_INTERP;
	            ely_props[i].depth = query_pts[i].coord[2];
	        }
    	    ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, UCVM_COORD_GEO_DEPTH, 30, ely_pts, ely_props);
        	for (i=0; i<30; i++) {
            	ucvm_interp_ely(0.0, taper_depth, UCVM_COORD_GEO_DEPTH, ely_pts+i, ely_props+i);
        	}
        	//Combine results
        	for (i=0; i<30; i++) {
	            //printf("Ely prop value at depth %f = %f\n", ely_props[i].depth, ely_props[i].cmb.vs);
	            if (ely_props[i].cmb.vs<vs_results[i]) {
	                vs_results[i] = ely_props[i].cmb.vs;
	            }
	        }
	        free(ely_props);
	        free(ely_pts);
		}

		double vs_sum = 0.0;
		for (i=0; i<30; i++) {
			vs_sum += 1.0/vs_results[i];
		}
	    float vs = 30.0/vs_sum;
		free(query_pts);
		free(query_data);
		free(vs_results);
        return vs;
}

//Determines second crossing (if it exists) of Z-values, with 10 m resolution, down to 50 km
float ucvm_zvalue(float lon, float lat, float vs_value, int taper_flag) {
	float resolution = 10.0;
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

	for (i=0; i<num_pts; i++) {
		vs_results[i] = data_pts[i].cmb.vs;
	}

    //If we need to include the taper
    if (taper_flag!=0) {
        printf("Using ifless taper.\n");
        //If we're using the taper, we need to initialize Ely
        //Similar code to that in ucvm-single_mpi.c
        if (ely_init==0) {
            ucvm_modelconf_t conf;
            ucvm_elygtl_model_init(UCVM_MAX_MODELS-1, &conf);
            ely_init = 1;
        }

		//We only need to use the taper down to the taper cutoff
		int num_taper_pts = (int)(taper_depth/resolution)+1;

        ucvm_data_t* ely_props = malloc(sizeof(ucvm_data_t)*num_taper_pts);
        ucvm_point_t* ely_pts = malloc(sizeof(ucvm_point_t)*num_taper_pts);
        memcpy(ely_props, data_pts, sizeof(ucvm_data_t)*num_taper_pts);
        memcpy(ely_pts, query_pts, sizeof(ucvm_point_t)*num_taper_pts);
        for (i=0; i<num_taper_pts; i++) {
            //Query this point, to get the crustal value at taper depth
            ely_pts[i].coord[2] = taper_depth;
        }
        if(ucvm_query(num_taper_pts, ely_pts, ely_props)!=UCVM_CODE_SUCCESS) {
            fprintf(stderr, "UCVM query failed.");
            exit(-3);
        }
        for (i=0; i<num_taper_pts; i++) {
            ely_props[i].domain = UCVM_DOMAIN_INTERP;
            ely_props[i].depth = query_pts[i].coord[2];
        }
        ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, UCVM_COORD_GEO_DEPTH, num_taper_pts, ely_pts, ely_props);
        for (i=0; i<num_taper_pts; i++) {
            ucvm_interp_ely(0.0, taper_depth, UCVM_COORD_GEO_DEPTH, ely_pts+i, ely_props+i);
        }
	        //Combine results
        for (i=0; i<num_taper_pts; i++) {
            //printf("Ely prop value at depth %f = %f\n", ely_props[i].depth, ely_props[i].cmb.vs);
            if (ely_props[i].cmb.vs<vs_results[i]) {
                vs_results[i] = ely_props[i].cmb.vs;
            }
        }
        free(ely_props);
        free(ely_pts);
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
		printf("Usage: %s <lat> <lon> <comma-separated model string> <output file> [comma-separated taper string]\n", argv[0]);
		exit(1);
	}
	float lat = atof(argv[1]);
	float lon = atof(argv[2]);
	char* model_string = argv[3];
	char* output_file = argv[4];

    initialize_ucvm(model_string);

	char* taper_models = NULL;
	int taper_flag = 0;
	if (argc>5) {
		taper_models = argv[5];
		//Since all the queries are for a single point, determine if this point is using a model that requires the taper
		taper_flag = test_for_taper(lon, lat, taper_models);
		printf("Taper flag = %d.\n", taper_flag);
	}

	initialize_ucvm(model_string);
	float vs30 = ucvm_vs30(lon, lat, taper_flag);
	float z10 = ucvm_zvalue(lon, lat, 1000.0, taper_flag);
	float z25 = ucvm_zvalue(lon, lat, 2500.0, taper_flag);

	FILE* fp_out = fopen(output_file, "w");
	fprintf(fp_out, "%.1f\n", vs30);
	fprintf(fp_out, "%.1f\n", z10);
	fprintf(fp_out, "%.1f\n", z25);
	fflush(fp_out);
	fclose(fp_out);
	
	return 0;
}
