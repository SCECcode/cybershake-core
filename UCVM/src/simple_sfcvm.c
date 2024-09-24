#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>
#include <math.h>

#include "ucvm.h"

int main(int argc, char** argv) {
	char* config_filename = "/lustre/orion/proj-shared/geo156/CyberShake/software/UCVM/ucvm_22.7.0_withSFCVM/conf/ucvm.conf";

	if (ucvm_init(config_filename)!=UCVM_CODE_SUCCESS) {
		printf("Failed to set up ucvm.\n");
		exit(-1);
	}

	if (ucvm_add_model("sfcvm") != UCVM_CODE_SUCCESS) {
		printf("Error retrieving SFCVM model.\n");
		exit(-2);
	}
	/*if (ucvm_add_model("cca") != UCVM_CODE_SUCCESS) {
    	printf("Error retrieving CCA model.\n");
		exit(-2);
    }
	if (ucvm_add_model(UCVM_MODEL_CVMSI) != UCVM_CODE_SUCCESS) {
	    printf("Error retrieving CVM-S4.26 model.\n");
        exit(-3);
	}*/

	if (ucvm_setparam(UCVM_MODEL_PARAM_QUERY_MODE, UCVM_COORD_GEO_DEPTH)!=UCVM_CODE_SUCCESS) {
		printf("Set query mode by depth failed.\n");
		exit(-3);
	}

	struct timeval start_tv;
	struct timeval end_tv;

	int num_x = 1;
	int num_y = 1;
	int num_z = 2;
	int num_pts = num_x*num_y*num_z;
	ucvm_point_t* pts = malloc(sizeof(ucvm_point_t)*(size_t)num_pts);
	ucvm_data_t* props = malloc(sizeof(ucvm_data_t)*(size_t)num_pts);

	float starting_lon = -121.886002;
	float starting_lat = 36.7304;
	//float starting_lon = -125.0;
	//float starting_lat = 39.0;
	//float starting_lon = -119.0;
	//float starting_lat = 38.5;
	//float starting_lon = -117.5;
	//float starting_lat = 36.5;
	float starting_depth = 20.0;
	int i, j, k;
	for(i=0; i<num_x; i++) {
		for(j=0; j<num_y; j++) {
			for (k=0; k<num_z; k++) {
				pts[i*num_y*num_z + j*num_z + k].coord[0] = starting_lon + 0.001*i;
				pts[i*num_y*num_z + j*num_z + k].coord[1] = starting_lat - 0.001*j;
				pts[i*num_y*num_z + j*num_z + k].coord[2] = starting_depth + 60.0*k;
			}
		}
	}

	gettimeofday(&start_tv, NULL);
	if (ucvm_query(num_pts, pts, props)!=UCVM_CODE_SUCCESS) {
		printf("Error querying UCVM.\n");
		exit(-4);
	}
	gettimeofday(&end_tv, NULL);
	printf("Runtime: %f\n", (end_tv.tv_sec - start_tv.tv_sec) + (end_tv.tv_usec - start_tv.tv_usec)/1000000.0);

	float min_vs = 400.0;
	if (props[0].cmb.vs < min_vs) {
		float vpvs_ratio = props[0].cmb.vp/props[0].cmb.vs;
		if (fabs(pts[0].coord[2])<0.01 || fabs(pts[0].coord[2]-20.0)<0.01) {
			float next_vpvs_ratio = props[1].cmb.vp/props[1].cmb.vs;
			if (next_vpvs_ratio>4.0) {
				next_vpvs_ratio = 4.0;
			}
			vpvs_ratio = next_vpvs_ratio;
		}
		props[0].cmb.vs = min_vs;
		props[0].cmb.vp = min_vs*vpvs_ratio;
	}

	for (i=0; i<num_pts; i++) {
		printf("(%f, %f, %f) -> Vp=%f, Vs=%f, rho=%f.\n", pts[i].coord[0], pts[i].coord[1], pts[i].coord[2], props[i].cmb.vp, props[i].cmb.vs, props[i].cmb.rho);
	}

	free(pts);
	free(props);
}
