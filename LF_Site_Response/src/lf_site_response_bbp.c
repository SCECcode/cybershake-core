#include "include.h"
#include "structure.h"
#include "functions_bbp.h"
#include "defs.h"

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


int main(int argc, char** argv) {
	char seis_in_filename[256];
	char seis_out_filename[256];
	char pga_filename[256];
	pga_filename[0] = '\0';
	char site_response_module[128] = "bssa2014";
	char vs30_model[128] = "cvmsi";
	int i;
	float lat;
	float lon;
	float vs30 = -1.0;
	float vref = 865.0;
	float vpga = vref;
	float pga_param[2] = {-1, -1};
	FILE* fp_in, *hf_fp_in;

    char pstring[MAX_PARAM_LENGTH];
	char** param_string = check_malloc(sizeof(char*) * MAX_PARAMS);
    for (i=0; i<MAX_PARAMS; i++) {
		param_string[i] = check_malloc(sizeof(char) * MAX_PARAM_LENGTH);
    }
    int num_params = 0;
	
	setpar(argc, argv);
	mstpar("seis_in", "s", seis_in_filename);
	mstpar("seis_out", "s", seis_out_filename);
	mstpar("slat", "f", &lat);
	mstpar("slon", "f", &lon);
	getpar("module", "s", site_response_module);
	getpar("vs30", "f", &vs30);
	getpar("vs30_model", "s", vs30_model);
	getpar("vref", "f", &vref);
	getpar("vpga", "f", &vpga);
	//If this is supplied, don't run wcc_getpeak()
	getpar("nspga", "f", &(pga_param[0]));
	getpar("ewpga", "f", &(pga_param[1]));
	//Pass in the HF seis to be used for determining PGA
	getpar("pga_infile", "s", pga_filename);
	endpar();

	if (vs30==-1.0) {
		ucvm_vsD(lon, lat, vs30_model, 30);
	}

	//If specified, parse PGA infile
	int pga_rv_count = 0;
	struct pga_entry* pga_entries;
	if (strlen(pga_filename)>0) {
		FILE* pga_fp_in = fopen(pga_filename, "r");
		int max_num_rvs = 100;
		pga_entries = check_malloc(sizeof(struct pga_entry)*max_num_rvs);
		while (!feof(pga_fp_in)) {
			fscanf(pga_fp_in, "%d %f %f %f\n", &(pga_entries[pga_rv_count].rup_var_id), &(pga_entries[pga_rv_count].pga[0]), &(pga_entries[pga_rv_count].pga[1]), &(pga_entries[pga_rv_count].pga[2]));
			pga_rv_count += 1;
			if (pga_rv_count>max_num_rvs) {
				max_num_rvs *= 2;
				pga_entries = check_realloc(pga_entries, sizeof(struct pga_entry)*max_num_rvs);
			}
		}
		fclose(pga_fp_in);
	}


	//Process all seismograms in the input file
	fp_in = fopen(seis_in_filename, "rb");
	struct seisheader in_header;
	struct statdata sdata;
	FILE* fp_out = fopen(seis_out_filename, "wb");
	while (fread(&in_header, sizeof(struct seisheader), 1, fp_in)>0) {
		float** seis = check_malloc(sizeof(float*)*2);
		for(i=0; i<2; i++) {
			seis[i] = check_malloc(sizeof(float)*in_header.nt);
			fread(seis[i], sizeof(float), in_header.nt, fp_in);
		}
		sdata.nt = in_header.nt;
		sdata.dt = in_header.dt;
	
		for (i=0; i<2; i++) {
			//to accel
	        add_param_reset(param_string, "diff=1", 1);
	        num_params = add_param(param_string, "integ=0");
	        integ_diff(num_params, param_string, seis[i], &sdata);
			num_params = 0;
			float pga = 0;
			if (strlen(pga_filename)>0) {
				//Get the PSA value for the corresponding HF
				int rv_id = in_header.rup_var_id;
				int pga_index = 0;
				while (pga_index<pga_rv_count) {
					if (pga_entries[pga_index].rup_var_id == rv_id) {
						break;
					}
					pga_index++;
				}
				if (pga_index==pga_rv_count) {
					printf("Couldn't find rupture variation ID %d in PGA file %s.  Aborting.\n", rv_id, pga_filename);
					exit(1);
				}
				pga = pga_entries[pga_index].pga[i];
			} else if (pga_param[i]<0) {
				pga = wcc_getpeak(num_params, param_string, seis[i], &sdata)/981.0;
			} else {
				pga = pga_param[i];
			}
			printf("PGA=%f\n", pga);
			sprintf(pstring, "vref=%f", vref);
            add_param_reset(param_string, pstring, 1);
            sprintf(pstring, "vsite=%f", vs30);
            add_param(param_string, pstring);
            sprintf(pstring, "vpga=%f", vpga);
            add_param(param_string, pstring);
            add_param(param_string, "model='bssa2014'");
            add_param(param_string, "flowcap=0.0");
            add_param(param_string, "fmax=50.0");
            add_param(param_string, "fhightop=20.0");
            add_param(param_string, "fmidbot=0.1");
            add_param(param_string, "fmin=0.05");
            sprintf(pstring, "pga=%f", pga);
            num_params = add_param(param_string, pstring);
            wcc_siteamp14(num_params, param_string, &(seis[i]), &sdata);
			//back to velocity
			add_param_reset(param_string, "diff=0", 1);
			num_params = add_param(param_string, "integ=1");
			integ_diff(num_params, param_string, seis[i], &sdata);
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

    for (i=0; i<MAX_PARAMS; i++) {
        free(param_string[i]);
    }
	free(param_string);

	if (strlen(pga_filename)>0) {
		free(pga_entries);
	}

	return 0;
}


