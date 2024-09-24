#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "structure.h"
#include "function_bbp.h"
#include "defs.h"
#include "include.h"

const int X_COMP_FLAG = 1;
const int Y_COMP_FLAG = 2;
const int Z_COMP_FLAG = 4;

extern void spectrad_(struct seisheader* header, int *nx, int *ny, int *npts, float *dt, char* seis_units, char* output_units, char* output_type, char* period, float *filter_high_hz, char* byteswap, char* input_file, char* output_file, float* seis, int* output_option, float* psa_data, int seis_units_len, int output_units_len, int output_type_len, int period_len, int byteswap_len, int input_file_len, int output_file_len);

int add_param_reset(char** param_string, char* param, int reset) {
    //getpar expects the arguments to start with the 1st parameter
    static int num_params = 1;
    if (reset==1) {
        num_params = 1;
    }
    if (strlen(param)>MAX_PARAM_LENGTH) {
        fprintf(stderr, "Parameter '%s' is %d long, which is longer than the maximum permitted parameter length of %d.  Change MAX_PARAM_LENGTH in defs.h.\n", param, strlen(param), MAX_PARAM_LENGTH);
        exit(2);
    }
    if (num_params>=MAX_PARAMS) {
        fprintf(stderr, "More parameters are being passed than are permitted.  Increase MAX_PARAMS in defs.h.\n");
        exit(3);
    }
    strcpy(param_string[num_params], param);
    num_params++;
    return num_params;
}

int add_param(char** param_string, char* param) {
    return add_param_reset(param_string, param, 0);
}

int main(int argc, char**argv) {
    char lf_seis_name[256];
    char hf_seis_name[256];
    char outfile[256];
    float match_freq = 1.0;
    float flo = 1.0e+15;
    int num_comps = 2;
	int i, j, rv, rc;
	int num_rup_vars = 1;
	FILE* seis_out;
	int debug = 0;

	setpar(argc, argv);
	//merging
    mstpar("lf_seis", "s", lf_seis_name);
    mstpar("hf_seis", "s", hf_seis_name);
    mstpar("seis_out", "s", outfile);
    getpar("freq", "f", &match_freq);
    getpar("comps", "d", &num_comps);
	getpar("num_rup_vars", "d", &num_rup_vars);
	getpar("debug", "d", &debug);

	//Modifications to support merging of rupture files
	//Iterate through LF file, looking for HF file match
	FILE* lf_in = fopfile(lf_seis_name, "rb");
	FILE* hf_in = fopfile(hf_seis_name, "rb");
	
	float*** lf_seis = NULL;
	float*** hf_seis = NULL;

	struct seisheader* lf_headers = check_malloc(sizeof(struct seisheader)*num_rup_vars);
	struct seisheader* hf_headers = check_malloc(sizeof(struct seisheader)*num_rup_vars);

	for (rv=0; rv<num_rup_vars; rv++) {
		//read and set up LF seismogram
		rc = fread(lf_headers+rv, sizeof(struct seisheader), 1, lf_in);
		if (rc!=1) {
			fprintf(stderr, "Error reading from file %s: supposed to read %d bytes, read %d.\n", lf_seis_name, sizeof(struct seisheader), rc);
			perror("");
			exit(10);
		}
		if (lf_seis==NULL) {
			lf_seis = check_malloc(sizeof(float**)*num_rup_vars);
			for (i=0; i<num_rup_vars; i++) {
				lf_seis[i] = check_malloc(sizeof(float*)*num_comps);
				for (j=0; j<num_comps; j++) {
					lf_seis[i][j] = check_malloc(sizeof(float)*lf_headers[rv].nt);
				}
			}
		}
		for (i=0; i<num_comps; i++) {
			rc = fread(lf_seis[rv][i], sizeof(float), lf_headers[rv].nt, lf_in);
			if (rc!=lf_headers[rv].nt) {
				fprintf(stderr, "Error reading from file %s: supposed to read %d bytes, read %d.\n", lf_seis_name, sizeof(float)*lf_headers[rv].nt, rc);
				perror("fread");
	                        exit(11);
			}
		}
		//Do HF seismogram
		rc = fread(hf_headers+rv, sizeof(struct seisheader), 1, hf_in);
        if (rc!=1) {
            fprintf(stderr, "Error reading from header from file %s.\n", hf_seis_name);
			perror("");
            exit(12);
        }
		printf("HF NT: %d\n", hf_headers[rv].nt);
        printf("HF num_comps: %d\n", num_comps);
        if (hf_seis==NULL) {
        	hf_seis = check_malloc(sizeof(float**)*num_rup_vars);
            for (i=0; i<num_rup_vars; i++) {
            	hf_seis[i] = check_malloc(sizeof(float*)*num_comps);
                for (j=0; j<num_comps; j++) {
                	hf_seis[i][j] = check_malloc(sizeof(float)*hf_headers[rv].nt);
                }
            }
         }
         for (i=0; i<num_comps; i++) {
         	rc = fread(hf_seis[rv][i], sizeof(float), hf_headers[rv].nt, hf_in);
            if (rc!=hf_headers[rv].nt) {
            	fprintf(stderr, "Error reading data from file %s.\n", hf_seis_name);
				perror("");
                exit(11);
            }
         }
		 //HF 3-component seismograms are in order 000,090,ver
		 //LF 3-component seismograms are in order ver,000,090
		 //If HF seismograms are 3-component, swap components to match LF
		 if (num_comps==3) {
			float* tmp = hf_seis[rv][2];
			hf_seis[rv][2] = hf_seis[rv][1];
			hf_seis[rv][1] = hf_seis[rv][0];
			hf_seis[rv][0] = tmp;
		}
	}
	fclose(lf_in);
	fclose(hf_in);

	seis_out = fopfile(outfile, "wb");

	//Get shared PSA arguments
    //PSA
    int nx, ny, npts;
    float dt;
    char seis_units[120];
    char output_units[120];
    char output_type[120];
    char period[120];
    char input_file[120];
    float filter_high_hz;
    char byteswap[120];
    char psa_output_file[256];
    int seis_units_len, output_units_len, output_type_len, period_len, byteswap_len, input_file_len, output_file_len;

    memset(seis_units, ' ', 120);
    memset(input_file, ' ', 120);
    memset(output_units, ' ', 120);
    memset(output_type, ' ', 120);
    memset(period, ' ', 120);
    memset(byteswap, ' ', 120);
    memset(psa_output_file, ' ', 256);
    input_file[0] = '\0';

    mstpar("simulation_out_pointsX","d",&nx);
    mstpar("simulation_out_pointsY","d",&ny);
    mstpar("simulation_out_timesamples","d",&npts);
    mstpar("simulation_out_timeskip","f",&dt);
    mstpar("surfseis_rspectra_seismogram_units","s",seis_units);
    mstpar("surfseis_rspectra_output_units","s",output_units);
    mstpar("surfseis_rspectra_output_type","s",output_type);
    mstpar("surfseis_rspectra_period","s",period);
    mstpar("surfseis_rspectra_apply_filter_highHZ","f",&filter_high_hz);
    mstpar("surfseis_rspectra_apply_byteswap","s",byteswap);
    mstpar("out","s",&psa_output_file);

	FILE* psa_out = fopen(psa_output_file, "wb");

        //replace null terminator with space
        //fortran expects all the strings to be space-terminated
        seis_units_len = strlen(seis_units);
        seis_units[seis_units_len] = ' ';
        input_file_len = strlen(input_file);
        input_file[input_file_len] = ' ';
        output_units_len = strlen(output_units);
        output_units[output_units_len] = ' ';
        output_type_len = strlen(output_type);
        output_type[output_type_len] = ' ';
	period_len = strlen(period);
	period[period_len] = ' ';

        output_file_len = strlen(psa_output_file);
        psa_output_file[output_file_len] = ' ';

        int output_option = PSA_NO_WRITE;
        float* psa_data = check_malloc(sizeof(float)*MAXPERIODS);

        //Put terminator back
        psa_output_file[output_file_len] = '\0';

        //RotD
        int run_rotd = 1;
        char rotd_filename[256];
        getpar("run_rotd", "d", &run_rotd);
        getpar("rotd_out", "s", rotd_filename);
	FILE* rotd_fp = NULL;
	int run_vert_rsp = 0;
	char vert_rsp_filename[256];
	FILE* vert_rsp_fp = NULL;
	if (run_rotd) {
		rotd_fp = fopfile(rotd_filename, "wb");
		if (num_comps==3) {
			getpar("run_vert_rsp", "d", &run_vert_rsp);
			getpar("vert_rsp_out", "s", vert_rsp_filename);
			vert_rsp_fp = fopfile(vert_rsp_filename, "wb");
		}
	}

	//Duration
	int run_duration = 1;
	char duration_filename[256];
	getpar("run_duration", "d", &run_duration);
	getpar("duration_out", "s", duration_filename);
	FILE* duration_fp = NULL;
	if (run_duration) {
		duration_fp = fopfile(duration_filename, "wb");
	}

	//Period-dependent duration
	int run_period_duration = 0;
	char period_duration_filename[256];
    getpar("run_period_durations", "d", &run_period_duration);
	getpar("period_duration_out", "s", period_duration_filename);
	endpar();
	FILE* period_duration_fp = NULL;
	if (run_period_duration) {
		period_duration_fp = fopfile(period_duration_filename, "wb");
	}
	
	
	//Now, operate on each LF rupture variation
	float** merged_seis = check_malloc(sizeof(float*)*num_comps);
	for (i=0; i<num_comps; i++) {
		merged_seis[i] = check_malloc(sizeof(float)*hf_headers[0].nt);
		//printf("allocated ms addr: %x\n", merged_seis[i]);
	}
	float* seis = NULL;

	for (rv=0; rv<num_rup_vars; rv++) {
		for (i=0; i<num_comps; i++) {
			memset(merged_seis[i], 0, sizeof(float)*hf_headers[0].nt);
		}
		int hf_index = rv;
		while (lf_headers[rv].source_id!=hf_headers[hf_index].source_id || lf_headers[rv].rupture_id!=hf_headers[hf_index].rupture_id || lf_headers[rv].rup_var_id!=hf_headers[hf_index].rup_var_id) {
			hf_index = (hf_index+1)%num_rup_vars;
			if (hf_index==rv) {
				fprintf(stderr, "Searched entire file, but can't find source %d, rupture %d, rup var %d in file %s, aborting.\n", lf_headers[rv].source_id, lf_headers[rv].rupture_id, lf_headers[rv].rup_var_id, hf_seis_name);
				exit(15);
			}
		}
		struct seisheader merged_header;
		merge_data(&(lf_seis[rv]), lf_headers[rv], hf_seis[hf_index], hf_headers[hf_index], match_freq, num_comps, merged_seis, &merged_header, debug);

        long bytes = fwrite(&merged_header, sizeof(struct seisheader), 1, seis_out);
        for (i=0; i<num_comps; i++) {
        	bytes = fwrite(merged_seis[i], sizeof(float), merged_header.nt, seis_out);
        }

		//Copy to 1-D array for use with PSA
		if (seis==NULL) {
			seis = check_malloc(sizeof(float)*num_comps*merged_header.nt);
		} else {
			memset(seis, 0, sizeof(float)*num_comps*merged_header.nt);
		}
        	for (i=0; i<num_comps; i++) {
        	    memcpy(seis+(i*merged_header.nt), merged_seis[i], sizeof(float)*merged_header.nt);
		}
		spectrad_(&merged_header, &nx, &ny, &npts, &dt, seis_units, output_units, output_type, period, &filter_high_hz, byteswap, input_file, psa_output_file, seis, &output_option, psa_data, seis_units_len, output_units_len, output_type_len, period_len, byteswap_len, input_file_len, output_file_len);
		bytes = fwrite(&merged_header, sizeof(struct seisheader), 1, psa_out);
	    bytes = fwrite(psa_data, sizeof(float), nx*ny*NUM_SCEC_PERIODS, psa_out);
		
		if (run_rotd) {
	    	int rc = rotd(&merged_header, seis, rotd_fp);
	        if (rc!=0) {
	        	fprintf(stderr, "Error in rotd code, aborting.\n");
	            exit(rc);
	        }
			if (run_vert_rsp) {
				rc = vert_rsp(&merged_header, seis, vert_rsp_fp);
				if (rc!=0) {
                	fprintf(stderr, "Error in vert_rsp code, aborting.\n");
                	exit(rc);
				}
			}
	    }
	
		if (run_duration) {
			int rc = duration(merged_header, merged_seis, duration_fp);
			if (rc!=0) {
				fprintf(stderr, "Error in duration code, aborting.\n");
				exit(rc);
			}
		}

		if (run_period_duration) {
			int rc = period_durations(merged_header, merged_seis, period_duration_fp);
            if (rc!=0) {
                fprintf(stderr, "Error in period-dependent duration code, aborting.\n");
                exit(rc);
            }
		}
	}

	for (i=0; i<num_comps; i++) {
		free(merged_seis[i]);
	}
	free(merged_seis);

/*
	//pass ref to merge_seis, since allocated in merge
        merge(lf_seis_name, hf_seis_name, outfile, match_freq, num_comps, &merged_seis, &merged_header);
*/
	fflush(seis_out);
	fclose(seis_out);
	fflush(psa_out);
	fclose(psa_out);
	if (run_rotd) {
		fflush(rotd_fp);
		fclose(rotd_fp);
		if (run_vert_rsp) {
			fflush(vert_rsp_fp);
			fclose(vert_rsp_fp);
		}
	}
	if (run_duration) {
		fflush(duration_fp);
		fclose(duration_fp);
	}

    if (run_period_duration) {
		fflush(period_duration_fp);
		fclose(period_duration_fp);
        python_finalize();
    }

	
	free(psa_data);

	for(rv=0; rv<num_rup_vars; rv++) {
		for(i=0; i<num_comps; i++) {
			free(lf_seis[rv][i]);
			free(hf_seis[rv][i]);
		}
		free(lf_seis[rv]);
		free(hf_seis[rv]);
	}
	free(lf_seis);
	free(hf_seis);

	free(lf_headers);
	free(hf_headers);
	free(seis);
	return 0;
}
