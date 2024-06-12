#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdlib.h>
#include <stdio.h>

#include "string.h"
#include "include.h"
#include "structure.h"
#include "function_bbp.h"
#include "defs.h"
#include "bbp_wcc_include.h"

/* This code acts as a wrapper for Rob's high frequency code, integrator, and merged into
* a single seismogram file in GRM format.
*/

extern void cs_params_();

extern void hb_high_(char* stat, float* slon, float* slat, char* local_vmod, char* outfile, float* tlen, float* dt, float* seis, int* nseg, float* elon, float* elat, float* targ_mag, int* nlskip, int* nx, int* ny, float* dx, float* dy, float* strike, float* dip, float* ravg, float* dtop, float* shypo, float* dhypo, float* sp, float* tr, float* ti, float* qfexp, char* velfile, float* stress_average, float* rvfac, int* seed, int* msite, int stat_len, int local_vmod_len, int outfile_len, int velfile_len);

float wcc_getpeak(int param_string_len, char** param_string, float* seis, struct statdata* head1);

int wcc_siteamp09(float* seis, int nt, float dt, float pga, float vs30);
int wcc_rotate(float seis[][mmv], int nt, float dt, float rot);

float get_rvfac(double mean_rvfac, double range_rvfac, int seed) {
	//char* filename = "/gpfs/alpine/proj-shared/geo112/CyberShake/software/bbp/bbp-19.8.0-python3/bbp/comps/hfsims_cfg.py";
	char* filename = "/work2/00349/scottcal/frontera/CyberShake/software/bbp/bbp-22.4.0/bbp/comps/hfsims_cfg.py";
	FILE* hfsims_cfg_fp = fopen(filename, "r");
	Py_Initialize();
	PyRun_SimpleFile(hfsims_cfg_fp, filename);
	PyObject* main_module = PyImport_AddModule("__main__");
	PyObject* global_dict = PyModule_GetDict(main_module);
	PyObject* expression = PyDict_GetItemString(global_dict, "calculate_rvfac");
	PyObject* args = PyTuple_New(3);
	PyTuple_SetItem(args, 0, PyFloat_FromDouble(mean_rvfac));
	PyTuple_SetItem(args, 1, PyFloat_FromDouble(range_rvfac));
	PyTuple_SetItem(args, 2, PyLong_FromLong((long)seed));
	PyObject* py_rvfac = PyObject_CallObject(expression, args);

	double rvfac = PyFloat_AsDouble(py_rvfac);
	Py_Finalize();
	printf("internal rvfac = %lf\n", rvfac);

	fclose(hfsims_cfg_fp);
	return (float)rvfac;
}

void write_grm(char* outname, float** seis, int nt) {
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

void write_grm_comp(char* outname, float* seis, int nt) {
	FILE* fp_out = fopen(outname, "wb");
	if (fp_out==NULL) {
        fprintf(stderr, "Unable to open file %s for writing.\n", outname);
        exit(1);
    }
	fwrite(seis, sizeof(float), nt, fp_out);
	fflush(fp_out);
	fclose(fp_out);
}

void write_grm_head(FILE* fp_out, float** seis, struct seisheader* header, int num_comps) {
    int i;
	fwrite(header, sizeof(struct seisheader), 1, fp_out);
    for (i=0; i<num_comps; i++) {
        fwrite(seis[i], sizeof(float), header->nt, fp_out);
    }
}

void hfsim(float** seisC, char* stat, float slon, float slat, char* local_vmod, FILE* output_fp, FILE* pga_output_fp, float vs30, struct seisheader* header, float modelrot, struct slipfile* sfile, int num_comps, int do_site_response, float vref, float vpga, int seed, int debug) {
	//parse cmd-line args
	int stat_len;
	int local_vmod_len;
	int output_len;
	float seis[mmv][3]; //needs to agree with mmv in hb_high
	int i;
	char output[256] = "output.tmp";
	char** param_string;
	//Must be provided
	float targ_mag = -1;
	float stress_average = 50.0;
	int msite = 1;

	int nt = header->nt;
	float dt = header->dt;
	float tlen = header->nt*header->dt;
	int nlskip = -99;

	param_string = check_malloc(sizeof(char*) * MAX_PARAMS);
	for (i=0; i<MAX_PARAMS; i++) {
		param_string[i] = check_malloc(sizeof(char) * MAX_PARAM_LENGTH);
	}

	//replace null terminator w/space for fortran
	stat_len = strlen(stat);
    stat[stat_len] = ' ';
	local_vmod_len = strlen(local_vmod);
	local_vmod[local_vmod_len] = ' ';
	output_len = strlen(output);
	output[output_len] = ' ';

	char velfile[256];
	memset(velfile, ' ', 256);
	velfile[0] = '-';
	velfile[1] = '1';
	int velfile_len = 256;

	//Reorder sp, tr, and ti arrays to be column-major for Fortran
	float *tmp_sp, *tmp_tr, *tmp_ti;
	tmp_sp = check_malloc(sizeof(float) * NP*NQ*LV);
	tmp_tr = check_malloc(sizeof(float) * NP*NQ*LV);
	tmp_ti = check_malloc(sizeof(float) * NP*NQ*LV);

	//Copy over, clear, then swap
	memcpy(tmp_sp, sfile->sp, sizeof(float)*NP*NQ*LV);
    memcpy(tmp_tr, sfile->tr, sizeof(float)*NP*NQ*LV);
    memcpy(tmp_ti, sfile->ti, sizeof(float)*NP*NQ*LV);
	memset(sfile->sp, 0, sizeof(float)*NP*NQ*LV);
	memset(sfile->tr, 0, sizeof(float)*NP*NQ*LV);
	memset(sfile->ti, 0, sizeof(float)*NP*NQ*LV);

	int ix, iy;

	//In fortran, fastest index is y, then x, then nseg
	//sddp(nseg,x,y) with dims sddp(lv,nq,np)
	//Also, column-major for fortran
	for(i=0; i<sfile->nseg; i++) {
		//Do separate loops for hopefully better caching
		printf("ny=%d\n", sfile->ny[i]);
		for(iy=0; iy<sfile->ny[i]; iy++) {
			for(ix=0; ix<sfile->nx[i]; ix++) {
				sfile->sp[iy*NQ*LV + ix*LV + i] = tmp_sp[i*NQ*NP + iy*sfile->nx[i] + ix];
			}
		}
		for(iy=0; iy<sfile->ny[i]; iy++) {
            for(ix=0; ix<sfile->nx[i]; ix++) {
				sfile->tr[iy*NQ*LV + ix*LV + i] = tmp_tr[i*NQ*NP + iy*sfile->nx[i] + ix];
            }
        }
		for(iy=0; iy<sfile->ny[i]; iy++) {
            for(ix=0; ix<sfile->nx[i]; ix++) {
				sfile->ti[iy*NQ*LV + ix*LV + i] = tmp_ti[i*NQ*NP + iy*sfile->nx[i] + ix];
            }
        }
	}
	free(tmp_sp);
	free(tmp_tr);
	free(tmp_ti);

	//Determine rvfac
	double mean_rvfac = 0.8;
	double range_rvfac = 0.05;
	sfile->rvfac = get_rvfac(mean_rvfac, range_rvfac, seed);
	if (debug) printf("rvfac = %f\n", sfile->rvfac);

	//Generate high-frequency seismogram
	if (debug) printf("Calculating HF seismogram.\n");
	/*stat_len = 64;
	local_vmod_len = 256;
	output_len = 256;
	velfile_len = 256;*/
	if (debug) printf("stat len=%d, local_vmod_len=%d, output_len=%d, velfile_len=%d\n", stat_len, local_vmod_len, output_len, velfile_len);
	//extern void hb_high_(char* stat, float* slon, float* slat, char* local_vmod, char* outfile, float* tlen, float* dt, float* seis, int* nseg, float* elon, float* elat, float* targ_mag, int* nx, int* ny, float* dx, float* dy, float* strike, float* dip, float* ravg, float* dtop, float* shypo, float* dhypo, float* sp, float* tr, float* ti, float* qfexp, float* stress_average, int* msite, int stat_len, int local_vmod_len, int outfile_len);
	//Initializes a bunch of hb_high variables to their CyberShake values
	cs_params_();
	hb_high_(stat, &slon, &slat, local_vmod, output, &tlen, &dt, &(seis[0][0]), &(sfile->nseg), sfile->elon, sfile->elat, &targ_mag, &nlskip, sfile->nx, sfile->ny, sfile->dx, sfile->dy, sfile->strike, sfile->dip, sfile->ravg, sfile->dtop, sfile->shypo, sfile->dhypo, sfile->sp, sfile->tr, sfile->ti, &(sfile->qfexp), velfile, &stress_average, &(sfile->rvfac), &seed, &msite, stat_len, local_vmod_len, output_len, velfile_len);
	if (debug) printf("hb_high completed.\n");

	//put a null terminator back in strings for C
	stat[stat_len] = '\0';
	local_vmod[local_vmod_len] = '\0';
	output[output_len] = '\0';
	//rearrange seismogram into C-friendly form
	//change order from 090, 000, ver to 000, 090, ver
	//While changing order, check that seismogram isn't all zeroes, w/o if statement
	//float seisC[3][mmv];
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

	if (debug) {
		write_grm("raw_hf", seisC, header->nt);
		FILE* acc_fp = fopen("accseis", "wb");
		write_grm_head(acc_fp, seisC, header, num_comps);
		fclose(acc_fp);
	}
	
	
	/*for each component:
	1) calculate peak PGA
	2) calculate site response
	3) filter
	4) integrate */
	struct statdata head1;
	head1.nt = nt;
	head1.dt = dt;
	char pstring[MAX_PARAM_LENGTH];
	int num_params = 0;
	char pga_string[256] = {'\0'};
	for (i=0; i<3; i++) {
		if (do_site_response) {
			//printf("Calculating peak for component %d.\n", i);
			/*float pga = wcc_getpeak(seisC[i], header->nt, header->dt)/981.0;
			if (debug) printf("Calculating site amp with pga=%f, vs30=%f.\n", pga, vs30);
			wcc_siteamp14(seisC[i], header->nt, header->dt, pga, vs30);*/
			float pga = wcc_getpeak(num_params, param_string, seisC[i], &head1)/981.0;
			sprintf(pga_string, "%s%.6f ", pga_string, pga);
			if (debug) printf("Calculating site amp with pga=%f, vs30=%f.\n", pga, vs30);
			sprintf(pstring, "vref=%f", 500.0);
			add_param_reset(param_string, pstring, 1);
			sprintf(pstring, "vsite=%f", vs30);
			add_param(param_string, pstring);
			sprintf(pstring, "vpga=%f", 500.0);
			add_param(param_string, pstring);
			add_param(param_string, "model='bssa2014'");
			add_param(param_string, "flowcap=0.0");
			add_param(param_string, "fmax=50.0");
			add_param(param_string, "fhightop=20.0");
			add_param(param_string, "fmidbot=0.1");
			add_param(param_string, "fmin=0.05");
			sprintf(pstring, "pga=%f", pga);
			num_params = add_param(param_string, pstring);
			int k;
			for (k=0; k<num_params; k++) {
				printf("Param %d: %s\n", k, param_string[k]);
			}
			wcc_siteamp14(num_params, param_string, &(seisC[i]), &head1);
			if (debug) {
				char filename[256];
				sprintf(filename, "raw_hf_with_amp_%d", i);
				write_grm_comp(filename, seisC[i], nt);
			}
			int order = 4;
			sprintf(pstring, "order=%d", order);
			add_param_reset(param_string, pstring, 1);
			//Same as BBP logic to calculate high-pass filter frequency
			float val = 1.0 / (2.0 * order);
			float fhi = exp(val*log(sqrt(2.0) - 1.0));
			sprintf(pstring, "fhi=%f", fhi);
			num_params = add_param(param_string, pstring);
			printf("Filtering with fhi=%f\n", fhi);
			wcc_tfilter(num_params, param_string, seisC[i], &head1);
		}
		printf("Integrating.\n");
		add_param_reset(param_string, "integ=1", 1);
		num_params = add_param(param_string, "diff=0");
		integ_diff(num_params, param_string, seisC[i], &head1);
		//integ_diff(1, seisC[i], nt, dt);
	}
	if (strlen(pga_string)>0) {
		//Write PGA string to file
		fprintf(pga_output_fp, "%d %s", header->rup_var_id, pga_string);
	}
	//remember seisC is 000, 090, ver
	//so must rotate so that X has heading 90 + modelrot to get to CyberShake
	//This is no longer correct; per emails with Rob Graves, jbsim3d rotates the output seismograms to 000/090/ver, in the mech_sgt() function.
	//wcc_rotate(seisC, nt, dt, modelrot+90);

	write_grm_head(output_fp, seisC, header, num_comps);

	for (i=0; i<MAX_PARAMS; i++) {
		free(param_string[i]);
	}
	free(param_string);
}
