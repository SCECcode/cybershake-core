#include <stdio.h>
#include <stdlib.h>
//#include "getpar.h"
//#include "string.h"

//#include "structure.h"
//#include "functions.h"

int main(int argc, char** argv) {
	char seis_in_filename[256];
	char seis_out_filename[256];
	char site_response_module[128] = "cb2014";
	char vs30_model[128] = "cvmsi";
	int i;
	float lat;
	float lon;
	float vs30 = -1.0;
	/*
	setpar(argc, argv);
	mstpar("seis_in", "s", seis_in_filename);
	mstpar("seis_out", "s", seis_out_filename);
	mstpar("slat", "f", &lat);
	mstpar("slon", "f", &lon);
	getpar("module", "s", site_response_module);
	getpar("vs30", "f", &vs30);
	getpar("vs30_model", "s", vs30_model);
	endpar();*/

	/*if (vs30==-1.0) {
		ucvm_vs30(lat, lon, vs30_model);
	}*/
	/*
	//Read in single seismogram
	FILE* fp_in = fopen(seis_in_filename, "rb");
	struct seisheader in_header;
	fread(&in_header, sizeof(struct seisheader), 1, fp_in);
	printf("dt=%f, nt=%d\n", in_header.dt, in_header.nt);
	float** seis = check_malloc(sizeof(float*)*2);
	for(i=0; i<2; i++) {
		seis[i] = check_malloc(sizeof(float)*in_header.nt);
		fread(seis[i], sizeof(float), in_header.nt, fp_in);
	}
	fclose(fp_in);*/
	/*
	for (i=0; i<1; i++) {
		//to accel
		//integ_diff(0, seis[i], in_header.nt, in_header.dt);
		float pga = wcc_getpeak(seis[i], in_header.nt, in_header.dt)/981.0;
		printf("PGA=%f\n", pga);
		wcc_siteamp14(seis[i], in_header.nt, in_header.dt, pga, vs30, site_response_module);
		//back to velocity
		//integ_diff(1, seis[i], in_header.nt, in_header.dt);
	}*/
	/*
	FILE* fp_out = fopen(seis_out_filename, "wb");
	fwrite(&in_header, sizeof(struct seisheader), 1, fp_out);
	for (i=0; i<2; i++) {
		fwrite(seis[i], sizeof(float), in_header.nt, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
	*/
	/*for(i=0; i<2; i++) {
		free(seis[i]);
	}
	free(seis);*/
	return 0;
}


