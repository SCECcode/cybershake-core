#include "structure.h"
#include "include.h"
#include "function_bbp.h"
#include "defs.h"
#include "bbp_srf_include.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

const int X_COMP_FLAG = 1;
const int Y_COMP_FLAG = 2;
const int Z_COMP_FLAG = 4;

void parse_rup_vars(char* rup_var_string, int num_rup_vars, struct rupture_variation* rup_vars);

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

void write_slipfile(struct slipfile sfile, char* outfile) {
	//Included for debugging
	FILE* fp_out = fopfile(outfile, "w");
	int i, j, k;
	fprintf(fp_out,"%d\n", sfile.nseg);
	for (i=0; i<sfile.nseg; i++) {
		fprintf(fp_out,"%10.4f %10.4f %5d %5d %8.2f %8.2f\n",sfile.elon[i],sfile.elat[i],sfile.nx[i],sfile.ny[i],sfile.dx[i],sfile.dy[i]);
		fprintf(fp_out,"%4.0f %4.0f %4.0f %8.2f %8.2f %8.2f\n",sfile.strike[i],sfile.dip[i],sfile.ravg[i],sfile.dtop[i],sfile.shypo[i],sfile.dhypo[i]);
		for (j=0; j<sfile.ny[i]; j++) {
			for (k=0; k<sfile.nx[i]; k++) {
				fprintf(fp_out, "%13.5e", sfile.sp[i*NQ*NP + j*NP + k]);
			}
			fprintf(fp_out, "\n");
		}
        for (j=0; j<sfile.ny[i]; j++) {
        	for (k=0; k<sfile.nx[i]; k++) {
            	fprintf(fp_out, "%13.5e", sfile.tr[i*NQ*NP + j*NP + k]);
            }
			fprintf(fp_out, "\n");
        }
        for (j=0; j<sfile.ny[i]; j++) {
        	for (k=0; k<sfile.nx[i]; k++) {
            	fprintf(fp_out, "%13.5e", sfile.ti[i*NQ*NP + j*NP + k]);
            }
			fprintf(fp_out, "\n");
        }
	}
	fclose(fp_out);
}

int main(int argc, char** argv) {
	char slipfile[256];
	char outfile[256];
	char infile[256];
	char** param_string;
	float avgstk = -1.0e+15;
	//for hfsim
	char stat[12];
	float slon, slat;
	char local_vmod[256];
	float vs30 = -1.0;
	float vref = 500.0;
	float vpga = 500.0;
	float tlen = 102.4;
	float dt = 0.025;
	float modelrot = -55;

	char rup_geom_file[256] = {0};
	int slip_id;
	int hypo_id;
	int i;

	int source_id = -1;
	int rupture_id = -1;
	int rup_var_id = -1;
	int det_max_freq = -1.0;
	int stoch_max_freq = 10.0;

	int do_site_response = 1;

	int num_rup_vars = -1;
	char rup_var_string[16384] = {0};
	int num_comps = 2;

	int debug = 0;

	param_string = check_malloc(sizeof(char*) * MAX_PARAMS);
	for (i=0; i<MAX_PARAMS; i++) {
		param_string[i] = check_malloc(sizeof(char) * MAX_PARAM_LENGTH);
	}

	memset(stat, ' ', 12);
	memset(local_vmod, ' ', 256);
	memset(outfile, ' ', 256);	

	setpar(argc, argv);
	getpar("debug", "d", &debug);
	//for srf2stoch
	//getpar("infile", "s", infile);
	//Can use rupture geometry file or SRF
	getpar("rup_geom_file", "s", rup_geom_file);
	getpar("infile", "s", infile);
	getpar("avgstk", "f", &avgstk);
	//for hfsim
    mstpar("stat", "s", stat);
    mstpar("slon", "f", &slon);
    mstpar("slat", "f", &slat);
    mstpar("vmod", "s", local_vmod);
    mstpar("outfile", "s", outfile);
    getpar("vs30", "f", &vs30);
	getpar("vref", "f", &vref);
	getpar("vpga", "f", &vpga);
    getpar("tlen", "f", &tlen);
    getpar("dt", "f", &dt);
    getpar("modelrot", "f", &modelrot);
	getpar("do_site_response", "d", &do_site_response);
	//Things to populate CyberShake header
    mstpar("source_id", "d", &source_id);
    mstpar("rupture_id", "d", &rupture_id);
    getpar("det_max_freq", "f", &det_max_freq);
    getpar("stoch_max_freq", "f", &stoch_max_freq);
	getpar("num_comps", "d", &num_comps);

	//If vs30 is not given, then use UCVM to determine it
	if (vs30==-1.0) {
		if (debug) printf("Retrieving Vs30 from UCVM.\n");
		char vs30_model[10];
		strcpy(vs30_model, "cvmsi");
		getpar("vs30_model", "s", vs30_model);
		vs30 = ucvm_vs30(slon, slat, vs30_model);
		if (debug) printf("vs30=%f\n", vs30);
	}

	//Support for multiple rupture variations
	getpar("num_rup_vars", "d", &num_rup_vars);
	printf("num_rvs = %d\n", num_rup_vars);
	struct rupture_variation* rup_vars;
	if (num_rup_vars==-1) {
		mstpar("rup_var_id", "d", &rup_var_id);
		mstpar("slip", "d", &slip_id);
		mstpar("hypo", "d", &hypo_id);
		num_rup_vars = 1;
		rup_vars = check_malloc(sizeof(struct rupture_variation)*num_rup_vars);
		rup_vars[0].rup_var_id = rup_var_id;
		rup_vars[0].slip_id = slip_id;
		rup_vars[0].hypo_id = hypo_id;
	} else {
		getpar("rup_vars", "s", rup_var_string);
		rup_vars = check_malloc(sizeof(struct rupture_variation)*num_rup_vars);
		//If rup_vars aren't defined, then we do all the rupture variations for this rupture
		//Figure out how many from rupgen_get_num_rv
		if (rup_var_string[0]==0) {
			if (rup_geom_file[0]==0) {
				fprintf(stderr, "Rupture geometry file was not specified, aborting.\n");
				exit(3);
			}
			rg_stats_t stats;
			set_memcached_server("localhost");
			rupgen_get_num_rv(rup_geom_file, &stats, RUPGEN_UNIFORM_HYPO);
			for (i=0; i<num_rup_vars; i++) {
				rup_vars[i].rup_var_id = i;
				rup_vars[i].slip_id = i;
				rup_vars[i].hypo_id = 0;
			}
		} else {
			parse_rup_vars(rup_var_string, num_rup_vars, rup_vars);
			set_memcached_server("localhost");
		}
	}
	endpar();

	FILE* fp_out = fopen(outfile, "wb");

	int rv;
	struct standrupformat srf;
	//Generate them sequentially
	for (rv=0; rv<num_rup_vars; rv++) {
		//float seis[3][mmv];
		float** seis = check_malloc(sizeof(float*) * 3);
		for (i=0; i<3; i++) {
			seis[i] = check_malloc(sizeof(float)*mmv);
		}
	    struct seisheader header;
	    strcpy(header.version, "12.10");
	    strcpy(header.site_name, stat);
	    header.source_id = source_id;
	    header.rupture_id = rupture_id;
	    header.rup_var_id = rup_vars[rv].rup_var_id;
	    header.det_max_freq = det_max_freq;
	    header.stoch_max_freq = stoch_max_freq;
		header.dt = dt;
		header.nt = (int)(tlen/dt+0.5);
		header.comps = X_COMP_FLAG | Y_COMP_FLAG;

		rg_stats_t stats;
		int rupture_seed = rupgen_get_rupture_seed(header.source_id, header.rupture_id);
		rupgen_genslip_seed(rup_geom_file, rv, 0, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05, rupture_seed);
		if (debug) {
			char srf_filename[128];
			sprintf(srf_filename, "%s.srf", rup_geom_file);
			write_srf(&srf, srf_filename, 0);
		}

	    struct slipfile sfile;
		//Set up expexted arguments
		int param_string_len;
		char pstring[MAX_PARAM_LENGTH];
		add_param(param_string, "target_dx=1.0");
		add_param(param_string, "target_dy=1.0");
		sprintf(pstring, "avgstk=%e", avgstk);
		param_string_len = add_param(param_string, pstring);

		srf2stoch(param_string_len, param_string, &srf, &sfile, debug);

		sprintf(slipfile, "%s.slip", infile);
		if (debug) write_slipfile(sfile, slipfile);
		//Need to do some rounding
		int i;
		printf("elon value from srf2stoch is %f\n", sfile.elon[0]);
		//Note that casting to int is not the same as floor for negative #s
		//So we are using floorf instead
		for(i=0; i<sfile.nseg; i++) {
			//sfile.elon[i] = ((int)(10000.0*sfile.elon[i]+0.5))/10000.0;
			sfile.elon[i] = (float)(floor(10000.0*sfile.elon[i]+0.5)/10000.0);
			//sfile.elat[i] = ((int)(10000.0*sfile.elat[i]+0.5))/10000.0;
			sfile.elat[i] = (float)(floor(10000.0*sfile.elat[i]+0.5)/10000.0);
		}
		printf("rounded elon value from srf2stoch is %f\n", sfile.elon[0]);
		for(i=0; i<sfile.nseg; i++) {
			/*sfile.ravg[i] = (int)(sfile.ravg[i]+0.5);
			sfile.dtop[i] = ((int)(100.0*sfile.dtop[i]+0.5))/100.0;
			sfile.shypo[i] = ((int)(100.0*sfile.shypo[i]+0.5))/100.0;
			sfile.dhypo[i] = ((int)(100.0*sfile.dhypo[i]+0.5))/100.0;*/
			sfile.ravg[i] = (float)(floor(sfile.ravg[i]+0.5));
			sfile.dtop[i] = (float)(floor(100.0*sfile.dtop[i]+0.5)/100.0);
			sfile.shypo[i] = (float)(floor(100.0*sfile.shypo[i]+0.5)/100.0);
			sfile.dhypo[i] = (float)(floor(100.0*sfile.dhypo[i]+0.5)/100.0);
		}
		for(i=0; i<NP*NQ*LV; i++) {
			sfile.sp[i] = (int)(100000.0*sfile.sp[i]+0.5)/100000.0;
			sfile.tr[i] = ((int)(100000.0*sfile.tr[i]+0.5))/100000.0;
			sfile.ti[i] = ((int)(100000.0*sfile.ti[i]+0.5))/100000.0;
		}
		for(i=0; i<LV; i++) {
			sfile.dx[i]=((int)(100.0*sfile.dx[i]+0.5))/100.0;
			sfile.dy[i]=((int)(100.0*sfile.dy[i]+0.5))/100.0;
		}
		sprintf(slipfile, "%s.slip.postround", infile);
		if (debug) write_slipfile(sfile, slipfile);
		//This is the default for LA Basin; other velocity models might be different
		sfile.qfexp = 0.6;

		//Calculate the variation seed out here, since we have all the parameters
		int variation_seed = rupgen_get_variation_seed(rup_geom_file, &stats, RUPGEN_UNIFORM_HYPO, header.source_id, header.rupture_id, rv, 0); 

		hfsim(seis, stat, slon, slat, local_vmod, fp_out, vs30, &header, modelrot, &sfile, num_comps, do_site_response, vref, vpga, variation_seed, debug);

		free(sfile.sp);
       	free(sfile.tr);
       	free(sfile.ti);
		for (i=0; i<3; i++) {
			free(seis[i]);
		}
		free(seis);
	}
	fflush(fp_out);
	fclose(fp_out);

	for (i=0; i<MAX_PARAMS; i++) {
		free(param_string[i]);
	}
	free(param_string);
}

void parse_rup_vars(char* rup_var_string, int num_rup_vars, struct rupture_variation* rup_vars) {
        //rup_var_string is in form (<rv_id>,<slip_id>,<hypo_id>);(....)
        int i;
        char* tok, *inner_tok;
        char* outer_save, *inner_save;
        tok = strtok_r(rup_var_string, "(;)", &outer_save);
        for (i=0; i<num_rup_vars; i++) {
                //tok points to rv<rv_id>,s<slip_id>,h<hypo_id>
                inner_tok = strtok_r(tok, ",", &inner_save);
                //parse rv
                rup_vars[i].rup_var_id = atoi(inner_tok);
                inner_tok = strtok_r(NULL, ",", &inner_save);
                //parse slip
                rup_vars[i].slip_id = atoi(inner_tok);
                //parse hypo
                inner_tok = strtok_r(NULL, ",", &inner_save);
                rup_vars[i].hypo_id = atoi(inner_tok);
                //Advance outer
                tok = strtok_r(NULL, "(;)", &outer_save);
        }
}
