#include "include.h"
#include "structure.h"
#include "functions.h"

int main(int argc, char** argv) {
	char seis_in_filename[256];
	char seis_out_filename[256];
	char site_response_module[128] = "cb2014";
	char vs30_model[128] = "cvmsi";
	int i;
	float lat;
	float lon;
	float vs30 = -1.0;
	float vref = 865.0;
	
	setpar(argc, argv);
	mstpar("seis_in", "s", seis_in_filename);
	mstpar("seis_out", "s", seis_out_filename);
	mstpar("slat", "f", &lat);
	mstpar("slon", "f", &lon);
	getpar("module", "s", site_response_module);
	getpar("vs30", "f", &vs30);
	getpar("vs30_model", "s", vs30_model);
	getpar("vref", "f", &vref);
	endpar();

	if (vs30==-1.0) {
		ucvm_vs30(lon, lat, vs30_model);
	}

	//Process all seismograms in the input file
	FILE* fp_in = fopen(seis_in_filename, "rb");
	struct seisheader in_header;
	FILE* fp_out = fopen(seis_out_filename, "wb");
	while (fread(&in_header, sizeof(struct seisheader), 1, fp_in)>0) {
		float** seis = check_malloc(sizeof(float*)*2);
		for(i=0; i<2; i++) {
			seis[i] = check_malloc(sizeof(float)*in_header.nt);
			fread(seis[i], sizeof(float), in_header.nt, fp_in);
		}
	
		for (i=0; i<2; i++) {
			//to accel
			integ_diff(0, seis[i], in_header.nt, in_header.dt);
			float pga = wcc_getpeak(seis[i], in_header.nt, in_header.dt)/981.0;
			printf("PGA=%f\n", pga);
			//must pad seismogram
			int nt_p2 = getnt_p2(in_header.nt);
			seis[i] = check_realloc(seis[i], nt_p2*sizeof(float));
			wcc_siteamp14(seis[i], in_header.nt, in_header.dt, pga, vs30, vref, site_response_module);
			//back to velocity
			integ_diff(1, seis[i], in_header.nt, in_header.dt);
		}
	
		fwrite(&in_header, sizeof(struct seisheader), 1, fp_out);
		for (i=0; i<2; i++) {
			fwrite(seis[i], sizeof(float), in_header.nt, fp_out);
		}
		fflush(fp_out);
		for(i=0; i<2; i++) {
                	free(seis[i]);
        	}
		free(seis);
	}
	fclose(fp_in);
	fclose(fp_out);	

	return 0;
}


