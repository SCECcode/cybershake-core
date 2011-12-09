#include <stdio.h>
#include <stdlib.h>

#include "structure.h"
#include "function.h"
/* 1) Filter HF and LF seismograms
*  2) Resample HF to LF dt
*  3) Add seismograms together
*/

int wcc_tfilter (float** seis, int order, float fhi, float flo, int phase, int ncomps, int nt_in, float dt_in);
int wcc_resamp_arbdt(float*** seis, int ncomps, float newdt, float olddt, int oldnt);
int wcc_add(float** lf_seis, float** hf_seis, float** merged_seis, int ncomp, float dt, int nt);

void write_grm(char* outname, float** seis, int nt, int ncomp) {
        int i;
        FILE* fp_out = fopen(outname, "wb");
        if (fp_out==NULL) {
                fprintf(stderr, "Unable to open file %s for writing.\n", outname);
                exit(1);
        }
        for (i=0; i<ncomp; i++) {
                fwrite(seis[i], sizeof(float), nt, fp_out);
        }
        fflush(fp_out);
        fclose(fp_out);
}

int main(int argc, char** argv) {
	char lf_seis_name[256];
	char hf_seis_name[256];
	char outfile[256];
	int order = 4;
	float match_freq = 1.0;
	float flo = 1.0e+15;
	int phase = 0;
	float lf_dt = 0.1;
	int lf_nt = 3000;
	float hf_dt = 0.025;
	int hf_nt = 12000;
	int num_comps = 2;

	int i;

	setpar(argc, argv);
	getpar("freq", "f", &match_freq);
	mstpar("lf_seis", "s", lf_seis_name);
	mstpar("hf_seis", "s", hf_seis_name);
	mstpar("outfile", "s", outfile);
	getpar("lf_dt", "f", &lf_dt);
	getpar("lf_nt", "d", &lf_nt);
	getpar("hf_dt", "f", &hf_dt);
	getpar("hf_nt", "d", &hf_nt);
	getpar("comps", "d", &num_comps);
	endpar();

	//get data from HF and LF seismograms
	float** lf_seis = (float**)check_malloc(sizeof(float*)*num_comps);
	for (i=0; i<num_comps; i++) {
		lf_seis[i] = (float*)check_malloc(sizeof(float)*lf_nt);
	}
	float** hf_seis = (float**)check_malloc(sizeof(float*)*num_comps);
	for (i=0; i<num_comps; i++) {
		hf_seis[i] = (float*)check_malloc(sizeof(float)*hf_nt);
	}
	float** merged_seis = (float**)check_malloc(sizeof(float*)*num_comps);
	for (i=0; i<num_comps; i++) {
		merged_seis[i] = (float*)check_malloc(sizeof(float)*hf_nt);
	}

	FILE* fp_in = fopen(lf_seis_name, "rb");
	for (i=0; i<num_comps; i++) {
		fread(lf_seis[i], sizeof(float), lf_nt, fp_in);
	}
	fclose(fp_in);
	fp_in = fopen(hf_seis_name, "rb");
	for (i=0; i<num_comps; i++) {
		fread(hf_seis[i], sizeof(float), hf_nt, fp_in);
	}
	fclose(fp_in);


	wcc_tfilter(hf_seis, order, match_freq, flo, phase, num_comps, hf_nt, hf_dt);
	//resample LF to HF
	//need to pass lf_seis by reference b/c reallocating pointer
	wcc_resamp_arbdt(&lf_seis, num_comps, hf_dt, lf_dt, lf_nt);
	wcc_add(lf_seis, hf_seis, merged_seis, num_comps, hf_dt, hf_nt);

	write_grm(outfile, merged_seis, hf_nt, num_comps);

	for (i=0; i<num_comps; i++) {
		free(lf_seis[i]);
		free(hf_seis[i]);
		free(merged_seis[i]);
	}
	free(lf_seis);
	free(hf_seis);
	free(merged_seis);
	return 0;
}
