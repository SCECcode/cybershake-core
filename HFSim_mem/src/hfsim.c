#include <stdlib.h>
#include <stdio.h>

#include "string.h"
#include "structure.h"

/* This code acts as a wrapper for Rob's high frequency code, integrator, and merged into
* a single seismogram file in GRM format.
*/

//extern void hb_high__(char* stat, int stat_len, float* slon, float* slat, char* slipfile, int slipfile_len, char* local_vmod, int local_vmod_len, char* outfile, int outfile_len, float* tlen, float* dt, float* seis);

extern void hb_high_(char* stat, float* slon, float* slat, char* local_vmod, char* outfile, float* tlen, float* dt, float* seis, int* nseg, float* elon, float* elat, int* nx, int* ny, float* dx, float* dy, float* strike, float* dip, float* ravg, float* dtop, float* shypo, float* dhypo, float* sp, float* tr, float* ti, float* qfexp, int stat_len, int local_vmod_len, int outfile_len);

void integ_diff(int integ_notdiff, float* seis, int nt, float dt);
float wcc_getpeak(float* seis, int nt, float dt);
int wcc_siteamp09(float* seis, int nt, float dt, float pga, float vs30);
int wcc_rotate(float seis[][mmv], int nt, float dt, float rot);

void write_grm(char* outname, float seis[][mmv], int nt) {
	int i;
	FILE* fp_out = fopen(outname, "wb");
	if (fp_out==NULL) {
		fprintf(stderr, "Unable to open file %s for writing.\n", outname);
		exit(1);
	}
	for (i=0; i<3; i++) {
		fwrite(seis[i], sizeof(float), nt, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
}

void hfsim(char* stat, float slon, float slat, char* local_vmod, char* output, float vs30, float tlen, float dt, float modelrot, struct slipfile* sfile, int do_site_response) {
	//parse cmd-line args
	int stat_len;
	int local_vmod_len;
	int output_len;
	float seis[mmv][3]; //needs to agree with mmv in hb_high
	int i;

	int nt = (int)(tlen/dt+.5); //round in case of funny business
	
	//replace null terminator w/space for fortran
	stat_len = strlen(stat);
        stat[stat_len] = ' ';
	local_vmod_len = strlen(local_vmod);
	local_vmod[local_vmod_len] = ' ';
	output_len = strlen(output);
	output[output_len] = ' ';

	//Generate high-frequency seismogram
	printf("Calculating HF seismogram.\n");
	hb_high_(stat, &slon, &slat, local_vmod, output, &tlen, &dt, &(seis[0][0]), &(sfile->nseg), sfile->elon, sfile->elat, sfile->nx, sfile->ny, sfile->dx, sfile->dy, sfile->strike, sfile->dip, sfile->ravg, sfile->dtop, sfile->shypo, sfile->dhypo, sfile->sp, sfile->tr, sfile->ti, &(sfile->qfexp), stat_len, local_vmod_len, output_len);
	//put a null terminator back in output
	output[output_len] = '\0';
	//rearrange seismogram into C-friendly form
	//change order from 090, 000, ver to 000, 090, ver
	//While changing order, check that seismogram isn't all zeroes, w/o if statement
	float seisC[3][mmv];
	int j;
	float all_zeroes[3] = {0.0, 0.0, 0.0};
	//swap 090 and 000
	for (i=0; i<2; i++) {
		int index = (i+1)%2;
		for (j=0; j<mmv; j++) {
			all_zeroes[i] += seis[j][i];
			seisC[index][j] = seis[j][i];
		}
		if (all_zeroes[i]==0.0) {
			fprintf(stderr, "Error in high frequency seismogram generation for component %d; seismogram is all zeros\n", i);
			exit(2);
		}
	}
	//ver
	all_zeroes[2] = 0.0;
	for (i=0; i<mmv; i++) {
		all_zeroes[2] += seis[i][2];
		seisC[2][i] = seis[i][2];
	}
        if (all_zeroes[2]==0.0) {
                fprintf(stderr, "Error in high frequency seismogram generation for component %d.\n", i);
                exit(2);
        }

	write_grm("raw_hf", seisC, nt);
	/*for each component:
	1) calculate peak PGA
	2) calculate site response
	3) integrate */
	for (i=0; i<3; i++) {
		if (do_site_response) {
			//printf("Calculating peak for component %d.\n", i);
			float pga = wcc_getpeak(seisC[i], nt, dt)/981.0;
			//printf("Calculating site amp with pga=%f, vs30=%f.\n", pga, vs30);
			wcc_siteamp09(seisC[i], nt, dt, pga, vs30);
		}
		//printf("Integrating.\n");
		integ_diff(1, seisC[i], nt, dt);
	}
	//remember seisC is 000, 090, ver
	//so must rotate so that X has heading 90 + modelrot to get to CyberShake
	//This is no longer correct; per emails with Rob Graves, jbsim3d rotates the output seismograms to 000/090/ver, in the mech_sgt() function.
	//wcc_rotate(seisC, nt, dt, modelrot+90);

	write_grm(output, seisC, nt);
}
