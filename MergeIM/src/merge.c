#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "structure.h"
#include "function.h"
/* 1) Filter HF and LF seismograms
*  2) Resample HF to LF dt
*  3) Add seismograms together
*/

int wcc_tfilter (float** seis, int order, float fhi, float flo, int phase, int ncomps, int nt_in, float dt_in);
int wcc_resamp_arbdt(float*** seis, int ncomps, float newdt, float olddt, int oldnt);
int wcc_add(float** lf_seis, float** hf_seis, float** merged_seis, int ncomp, float dt, int nt);

void write_grm_head(char* outname, float** seis, int ncomp, struct seisheader* header) {
        int i;
        FILE* fp_out = fopen(outname, "wb");
        if (fp_out==NULL) {
                fprintf(stderr, "Unable to open file %s for writing.\n", outname);
                exit(1);
        }
        fwrite(header, sizeof(struct seisheader), 1, fp_out);
        for (i=0; i<ncomp; i++) {
                fwrite(seis[i], sizeof(float), header->nt, fp_out);
        }
        fflush(fp_out);
        fclose(fp_out);
}

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

void merge(char* lf_seis_name, char* hf_seis_name, char* outfile, float match_freq, int num_comps, float*** merged_seis, struct seisheader* merged_header) {
	int order = 4;
	int phase = 0;
	float flo = 1.0e+15;	

	int i;

	//get data from HF and LF seismograms
	float** lf_seis = (float**)check_malloc(sizeof(float*)*num_comps);
	float** hf_seis = (float**)check_malloc(sizeof(float*)*num_comps);

	FILE* fp_in = fopen(lf_seis_name, "rb");
	struct seisheader lf_header;
	fread(&lf_header, 1, sizeof(struct seisheader), fp_in);
	for (i=0; i<num_comps; i++) {
                lf_seis[i] = (float*)check_malloc(sizeof(float)*lf_header.nt);
		fread(lf_seis[i], sizeof(float), lf_header.nt, fp_in);
	}
	fclose(fp_in);

	fp_in = fopen(hf_seis_name, "rb");
        struct seisheader hf_header;
        fread(&hf_header, 1, sizeof(struct seisheader), fp_in);
	for (i=0; i<num_comps; i++) {
		hf_seis[i] = (float*)check_malloc(sizeof(float)*hf_header.nt);
		fread(hf_seis[i], sizeof(float), hf_header.nt, fp_in);
	}
	fclose(fp_in);

	//Check header versions
	if (strcmp(lf_header.version, "12.10")!=0) {
		fprintf(stderr, "Header of seismogram %s is version %s.  We can only operate on version 12.10, aborting.\n", lf_seis_name, lf_header.version);
		exit(1);
	}
	if (strcmp(hf_header.version, "12.10")!=0) {
                fprintf(stderr, "Header of seismogram %s is version %s.  We can only operate on version 12.10, aborting.\n", hf_seis_name, hf_header.version);
                exit(1);
        }

	//Make sure these seismograms are for the same site
	if (strcmp(lf_header.site_name, hf_header.site_name)!=0) {
		fprintf(stderr, "Sites of seismograms don't match: %s is for site %s but %s is for site %s, aborting.\n", lf_seis_name, lf_header.site_name, hf_seis_name, hf_header.site_name);
		exit(2);
	}
	//Check source, rupture, rup_var
	if (lf_header.source_id!=hf_header.source_id || lf_header.rupture_id!=hf_header.rupture_id || lf_header.rup_var_id!=hf_header.rup_var_id) {
		fprintf(stderr, "Events don't match: %s is for source %d, rupture %d, rup var %d; %s is for source %d, rupture %d, rup var %d, aborting.\n", lf_seis_name, lf_header.source_id, lf_header.rupture_id, lf_header.rup_var_id, hf_header.source_id, hf_header.rupture_id, hf_header.rup_var_id);
		exit(3);
	}

	if (lf_header.det_max_freq<match_freq) {
		fprintf(stderr, "Seismograms are supposed to be combined at %f Hz, but the low frequency seismogram is only valid to %f Hz.  Aborting.\n", match_freq, lf_header.det_max_freq);
		exit(4);
	}

	*merged_seis = check_malloc(sizeof(float*)*num_comps);
	for (i=0; i<num_comps; i++) {
		(*merged_seis)[i] = check_malloc(sizeof(float)*hf_header.nt);
	}

	wcc_tfilter(hf_seis, order, match_freq, flo, phase, num_comps, hf_header.nt, hf_header.dt);
	//resample LF to HF
	//need to pass lf_seis by reference b/c reallocating pointer
	wcc_resamp_arbdt(&lf_seis, num_comps, hf_header.dt, lf_header.dt, lf_header.nt);
	wcc_add(lf_seis, hf_seis, *merged_seis, num_comps, hf_header.dt, hf_header.nt);

	//Set up merged header
	strcpy(merged_header->version, lf_header.version);
	strcpy(merged_header->site_name, lf_header.site_name);
	memset(merged_header->padding, 0, 8);
	merged_header->source_id = lf_header.source_id;
	merged_header->rupture_id = lf_header.rupture_id;
	merged_header->rup_var_id = lf_header.rup_var_id;
	merged_header->dt = hf_header.dt;
	merged_header->nt = hf_header.nt;
	merged_header->comps = lf_header.comps;
	merged_header->det_max_freq = match_freq;
	merged_header->stoch_max_freq = hf_header.stoch_max_freq;

	write_grm_head(outfile, *merged_seis, num_comps, merged_header);

	for (i=0; i<num_comps; i++) {
		free(lf_seis[i]);
		free(hf_seis[i]);
	}
	free(lf_seis);
	free(hf_seis);
}

void merge_data(float** lf_seis, struct seisheader lf_header, float** hf_seis, struct seisheader hf_header, float match_freq, int num_comps, float** merged_seis, struct seisheader* merged_header) {
        int order = 4;
        int phase = 0;
        float flo = 1.0e+15;

	wcc_tfilter(hf_seis, order, match_freq, flo, phase, num_comps, hf_header.nt, hf_header.dt);
        //resample LF to HF
        //need to pass lf_seis by reference b/c reallocating pointer
        wcc_resamp_arbdt(&lf_seis, num_comps, hf_header.dt, lf_header.dt, lf_header.nt);
        wcc_add(lf_seis, hf_seis, merged_seis, num_comps, hf_header.dt, hf_header.nt);

        //Set up merged header
        strcpy(merged_header->version, lf_header.version);
        strcpy(merged_header->site_name, lf_header.site_name);
        memset(merged_header->padding, 0, 8);
        merged_header->source_id = lf_header.source_id;
        merged_header->rupture_id = lf_header.rupture_id;
        merged_header->rup_var_id = lf_header.rup_var_id;
        merged_header->dt = hf_header.dt;
        merged_header->nt = hf_header.nt;
        merged_header->comps = lf_header.comps;
        merged_header->det_max_freq = match_freq;
        merged_header->stoch_max_freq = hf_header.stoch_max_freq;
}

