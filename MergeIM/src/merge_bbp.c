#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "structure.h"
#include "function_bbp.h"
#include "defs.h"
#include "include.h"

#include "bbp_wcc_include.h"
/* 1) Filter LF seismogram
*  2) Resample LF to HF dt
*  3) Add seismograms together
*/

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

void merge_data(float*** lf_seis, struct seisheader lf_header, float** hf_seis, struct seisheader hf_header, float match_freq, int num_comps, float** merged_seis, struct seisheader* merged_header, int debug) {
	int num_params=0;
	int i;
	char** param_string = check_malloc(sizeof(char*) * MAX_PARAMS);
    for (i=0; i<MAX_PARAMS; i++) {
        param_string[i] = check_malloc(sizeof(char) * MAX_PARAM_LENGTH);
    }
	char pstring[MAX_PARAM_LENGTH];
	//We no longer need to high-pass filter the stochastic, since this is now done during HFSim
	//wcc_tfilter(hf_seis, order, match_freq, flo, phase, num_comps, hf_header.nt, hf_header.dt);
	
	FILE* filt_fp_out, *resamp_lf_fp_out;
	filt_fp_out = resamp_lf_fp_out = NULL;
	if (debug) {
		char filtered_filename[256];
		sprintf(filtered_filename, "%s_LF_filtered.grm", lf_header.site_name);
		filt_fp_out = fopen(filtered_filename, "wb");
		fwrite(&lf_header, sizeof(struct seisheader), 1, filt_fp_out);
		char resamp_filename[256];
		sprintf(resamp_filename, "%s_LF_resamp.grm", lf_header.site_name);
		resamp_lf_fp_out = fopen(resamp_filename, "wb");
		fwrite(&hf_header, sizeof(struct seisheader), 1, resamp_lf_fp_out);
	}

	for (i=0; i<num_comps; i++) {
		//Low-pass filter low frequency seismogram
		struct statdata lf_head;
		lf_head.nt = lf_header.nt;
		lf_head.dt = lf_header.dt;
		struct statdata hf_head;
		hf_head.nt = hf_header.nt;
		hf_head.dt = hf_header.dt;
		lf_head.hr = hf_head.hr = 0;
		lf_head.min = hf_head.min = 0;
		lf_head.sec = hf_head.sec = 0;
		lf_head.edist = hf_head.edist = 0.0;
		lf_head.az = hf_head.az = 0.0;
		lf_head.baz = hf_head.baz = 0.0;
		int order = 4;
		sprintf(pstring, "order=%d", order);
		add_param_reset(param_string, pstring, 1);
		float val = -1.0 / (2.0 * order);
		float flo = match_freq*exp(val*log(sqrt(2.0)-1.0));
		sprintf(pstring, "flo=%f", flo);
		printf("flo=%f\n", flo);
		num_params = add_param(param_string, pstring);
		int j;
		for (j=1; j<num_params; j++) {
			printf("param %d: %s\n", j, param_string[j]);
		}
		wcc_tfilter(num_params, param_string, (*lf_seis)[i], &lf_head);
	
		if (debug) {
			//Write out filtered LF seismograms
			fwrite((*lf_seis)[i], sizeof(float), lf_header.nt, filt_fp_out);
		}
	
	    //resample LF to HF
	    //need to pass lf_seis by reference b/c reallocating pointer

		sprintf(pstring, "newdt=%f", hf_header.dt);
		num_params = add_param_reset(param_string, pstring, 1);
		wcc_resamp_arbdt(num_params, param_string, &((*lf_seis)[i]), &lf_head);

		printf("Resampled LF has dt=%f, nt=%d\n", lf_head.dt, lf_head.nt);
	
		if (fabs(lf_head.dt-hf_head.dt)>0.0001) {
			printf("Error: resampled lf has dt=%f, but hf has dt=%f, aborting.\n", lf_head.dt, hf_head.dt);
			exit(1);
		}

		if (debug) {
			//Write out resampled LF seismogram
			fwrite((*lf_seis)[i], sizeof(float), lf_head.nt, resamp_lf_fp_out);
		}

		float* hs = hf_seis[i];
		float* ms = merged_seis[i];
		struct statdata merged_head;

        printf("ms addr in merge_bbp: %x\n", ms);
		printf("lf_head.nt=%d, hf_head.nt=%d\n", lf_head.nt, hf_head.nt);

		num_params = 0;
		wcc_add(num_params, param_string, (*lf_seis)[i], &lf_head, hs, &hf_head, ms, &merged_head);
	}

	if (debug) {
		fflush(filt_fp_out);
		fclose(filt_fp_out);
		fflush(resamp_lf_fp_out);
		fclose(resamp_lf_fp_out);
	}

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

	for (i=0; i<MAX_PARAMS; i++) {
        free(param_string[i]);
    }
	free(param_string);
}

