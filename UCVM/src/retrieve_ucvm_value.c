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
	if (ucvm_init("/lustre/orion/geo156/proj-shared/CyberShake/software/UCVM/ucvm_22.7.0_withSFCVM/conf/ucvm.conf")!=UCVM_CODE_SUCCESS) {
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
				//ifless_taper = 1;
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

            if (ucvm_setparam(UCVM_MODEL_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
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
	
	        if (ucvm_setparam(UCVM_MODEL_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
	        	fprintf(stderr, "Set query mode by depth failed.\n");
	            exit(-2);
	        }
			tok = strtok(NULL, ",");
		}
	}
	printf("Initialization complete.\n");
}

//Retrieve point from UCVM using the API (and a possible taper).
void query(float lon, float lat, float depth, char* model) {
	if (!ucvm_initialized) {
		initialize_ucvm(model);
		ucvm_initialized = 1;
	}

	//Query without taper	
	ucvm_point_t query_pt;
	ucvm_data_t query_data;
	query_pt.coord[0] = lon;
	query_pt.coord[1] = lat;
	query_pt.coord[2] = depth;
	if (ucvm_query(1, &query_pt, &query_data)!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "UCVM query failed.\n");
			exit(-3);
	}
	double vp_result = query_data.cmb.vp;
	double vs_result = query_data.cmb.vs;
	double rho_result = query_data.cmb.rho;
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
		ucvm_data_t ely_prop;
		ucvm_point_t ely_pt;
		memcpy(&ely_prop, &query_data, sizeof(ucvm_data_t));
		memcpy(&ely_pt, &query_pt, sizeof(ucvm_point_t));
		ely_pt.coord[2] = taper_depth;
		if(ucvm_query(1, &ely_pt, &ely_prop)!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "UCVM query failed.");
			exit(-3);
		}
		ely_prop.domain = UCVM_DOMAIN_INTERP;
		ely_prop.depth = query_pt.coord[2];
		ucvm_elygtl_model_query(UCVM_MAX_MODELS-1, UCVM_COORD_GEO_DEPTH, 1, &ely_pt, &ely_prop);
		ucvm_interp_ely(0.0, taper_depth, UCVM_COORD_GEO_DEPTH, &ely_pt, &ely_prop);
		
		//printf("Ely prop value at depth %f = %f\n", ely_props[i].depth, ely_props[i].cmb.vs);
		if (ely_prop.cmb.vs<vs_result) {
			vs_result = ely_prop.cmb.vs;
			vp_result = ely_prop.cmb.vp;
			rho_result = ely_prop.cmb.rho;	
		}
	}
	printf("Vp = %f, Vs = %f, rho = %f\n", vp_result, vs_result, rho_result);	
}


//Need Vs30, VsD5H, Vs5H
//where H is the grid spacing
//Vs5H = slowness average over 5H, like Vs30
//VsD5H = discrete average over increments of H.
//So VsD500 = 5/(0.5/Vs(Z=0) + 1/Vs(Z=100) + ... + 1/Vs(Z=400) + 0.5/Vs(Z=500))
int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <lon> <lat> <depth> <model>\n", argv[0]);
		return 1;
	}
	float lon = atof(argv[1]);
	float lat = atof(argv[2]);
	float depth = atof(argv[3]);
	char* model = argv[4];
	query(lon, lat, depth, model);
	return 0;
}

