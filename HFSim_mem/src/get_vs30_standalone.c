#include "include.h"

#include "ucvm.h"

int ucvm_initialized = 0;

//Retrieve vs30 from UCVM, for the given model.
float ucvm_vs30(float lon, float lat, char* model) {
	//The algorithm is:
	//Vs30 = 30 / sum( 1 / (Vs sampled from [0.5, 29.5] at 1 meter increments, for 30 values) )
	if (!ucvm_initialized) {
		if (ucvm_init("/projects/sciteam/jmz/CyberShake/software/UCVM/ucvm_14.3.0/conf/ucvm.conf")!=UCVM_CODE_SUCCESS) {
			fprintf(stderr, "Failed to init UCVM, aborting.");
			exit(3);
		}

		char* model_string;
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
		ucvm_initialized = 1;
	}

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


int main(int argc, char** argv) {
	float lon = atof(argv[1]);
	float lat = atof(argv[2]);
	char* model = argv[3];

	float vs30 = ucvm_vs30(lon, lat, model);
	printf("Vs30 = %f\n", vs30);
	return 0;
}

