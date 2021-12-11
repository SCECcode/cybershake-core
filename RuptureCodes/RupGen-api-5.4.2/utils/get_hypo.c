#include <stdio.h>
#include <stdlib.h>

#include "rupgen_api.h"

int main(int argc, char** argv) {
	if (argc<2) {
		printf("Usage: %s <srf file>\n", argv[0]);
		exit(1);
	}
    char* srf_filename = argv[1];
	struct standrupformat srf;
	read_srf(&srf, srf_filename, 0);
	float hlon, hlat, hdep;
	float tmin;
	struct srf_apointvalues *apval_ptr1;
	tmin = 1.0e+15;
	apval_ptr1 = srf.srf_apnts.apntvals;
	for (int ip=0; ip<srf.srf_apnts.np; ip++) {
		if(apval_ptr1[ip].tinit < tmin) {
			hlon = apval_ptr1[ip].lon;
			hlat = apval_ptr1[ip].lat;
			hdep = apval_ptr1[ip].dep;
			tmin = apval_ptr1[ip].tinit;
		}
	}
	printf("%.5f\t%.5f\t%.5f\n",hlon,hlat,hdep);
    free_srf_ptrs(&srf);
}

