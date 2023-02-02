#include "structure.h"
#include "include.h"
#include "function_bbp.h"
#include "defs.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

const int X_COMP_FLAG = 1;
const int Y_COMP_FLAG = 2;
const int Z_COMP_FLAG = 4;

const int RUP_GEOM_MODE = 1;
const int SRF_MODE = 2;

void parse_rup_vars(char* rup_var_string, int num_rup_vars, struct rupture_variation* rup_vars, int rvfrac_seed);

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
				fprintf(fp_out, "%13.5e", sfile.sp[i*NQ*NP + j*sfile.nx[i] + k]);
			}
			fprintf(fp_out, "\n");
		}
        for (j=0; j<sfile.ny[i]; j++) {
        	for (k=0; k<sfile.nx[i]; k++) {
            	fprintf(fp_out, "%13.5e", sfile.tr[i*NQ*NP + j*sfile.nx[i] + k]);
            }
			fprintf(fp_out, "\n");
        }
        for (j=0; j<sfile.ny[i]; j++) {
        	for (k=0; k<sfile.nx[i]; k++) {
            	fprintf(fp_out, "%13.5e", sfile.ti[i*NQ*NP + j*sfile.nx[i] + k]);
            }
			fprintf(fp_out, "\n");
        }
	}
	fclose(fp_out);
}

int main(int argc, char** argv) {
	char slipfile[256];
	char outfile[256];
	char infile[256] = {'\0'};
	//Stores PGAs needed for site response
	char pga_outfile[256] = {'\0'};
	char** param_string;
	float avgstk = -1.0e+15;
	//for hfsim
	char stat[64];
	float slon, slat;
	char local_vmod[256];
	float vs30 = -1.0;
	float vref = 500.0;
	float vpga = 500.0;
	float tlen = 102.4;
	float dt = 0.025;
	float modelrot = -55;

	char rup_geom_file[256] = {'\0'};
	int slip_id;
	int hypo_id;
	int i;

	int source_id = -1;
	int rupture_id = -1;
	int rup_var_id = -1;
	int det_max_freq = -1.0;
	int stoch_max_freq = 10.0;

	float target_dx = 1.0;
	float target_dy = 1.0;

	int do_site_response = 1;

	int num_rup_vars = -1;
	char rup_var_string[16384] = {0};
	int num_comps = 2;

	int debug = 0;

	int mode = RUP_GEOM_MODE;
	int srf_seed;
	int rvfrac_seed = 0;

	param_string = check_malloc(sizeof(char*) * MAX_PARAMS);
	for (i=0; i<MAX_PARAMS; i++) {
		param_string[i] = check_malloc(sizeof(char) * MAX_PARAM_LENGTH);
	}

	memset(stat, ' ', 64);
	memset(local_vmod, ' ', 256);
	memset(outfile, ' ', 256);	

	setpar(argc, argv);
	getpar("debug", "d", &debug);
	//for srf2stoch
	//getpar("infile", "s", infile);
	//Can use rupture geometry file or SRF
	//Specify rup_geom_file OR infile
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
	getpar("pga_outfile", "s", &pga_outfile);
	getpar("target_dx", "f", &target_dx);
	getpar("target_dy", "f", &target_dy);
	//If we're specifying rvfrac and seed
	getpar("rvfrac_seed_given", "d", &rvfrac_seed);

	//If vs30 is not given, then use UCVM to determine it
	if (vs30==-1.0) {
		if (debug) printf("Retrieving Vs30 from UCVM.\n");
		char vs30_model[10];
		strcpy(vs30_model, "cvmsi");
		getpar("vs30_model", "s", vs30_model);
		vs30 = ucvm_vs30(slon, slat, vs30_model);
		if (debug) printf("vs30=%f\n", vs30);
	}

	struct rupture_variation* rup_vars;
	if (rup_geom_file[0]!='\0') {
		//If rup_geom_file is passed, assume we're generating the SRFs
		//Support for multiple rupture variations
		mode = RUP_GEOM_MODE;
		getpar("num_rup_vars", "d", &num_rup_vars);
		printf("num_rvs = %d\n", num_rup_vars);
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
				parse_rup_vars(rup_var_string, num_rup_vars, rup_vars, rvfrac_seed);
				set_memcached_server("localhost");
			}
		}
	} else if (infile[0]!='\0') {
		//SRF file was specified
		printf("Running in SRF mode.\n");
		mode = SRF_MODE;
		num_rup_vars = 1;
		mstpar("rup_var_id", "d", &rup_var_id);
		mstpar("srf_seed", "d", &srf_seed);
		rup_vars = check_malloc(sizeof(struct rupture_variation)*num_rup_vars);
        rup_vars[0].rup_var_id = rup_var_id;
		//Force calculation of rvfrac in hfsims
		rup_vars[0].rvfrac = -1.0;
	} else {
		printf("Error: must specify either rup_geom_file or infile.  Aborting.\n");
		exit(2);
	}
	endpar();

	FILE* fp_out = fopen(outfile, "wb");
	FILE* pga_fp_out = NULL;
	if (strlen(pga_outfile)>0) { 
		pga_fp_out = fopen(pga_outfile, "w");
	} else {
		sprintf(pga_outfile, "%s.pga", outfile);
		pga_fp_out = fopen(pga_outfile, "w");
	}

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
	
			if (mode==RUP_GEOM_MODE) {
				rg_stats_t stats;
				if (rvfrac_seed==1) {
					//Use rvfrac and seed from arguments
					int seed = rup_vars[rv].seed;
					float rvfrac = rup_vars[rv].rvfrac;
					char params[256];
					sprintf(params, "seed=%d rvfrac=%f use_unmodified_seed=1", seed, rvfrac);
					rupgen_genslip_with_params(rup_geom_file, rup_vars[rv].slip_id, rup_vars[rv].hypo_id, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05, params);
				} else {
					//No seed/rvfrac passed
					int rupture_seed = rupgen_get_rupture_seed(header.source_id, header.rupture_id);
					rupgen_genslip_seed(rup_geom_file, rup_vars[rv].slip_id, rup_vars[rv].hypo_id, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05, rupture_seed);
				}
				if (debug) {
					char srf_filename[256];
					sprintf(srf_filename, "%s.%d.srf", rup_geom_file, rup_vars[rv].rup_var_id);
					_write_srf(&srf, srf_filename, 0);
				}
			} else if (mode==SRF_MODE) {
				_read_srf(&srf, infile, 0);
			}
	
		    struct slipfile sfile;
			//Set up expexted arguments
			int param_string_len;
			char pstring[MAX_PARAM_LENGTH];
			sprintf(pstring, "target_dx=%.1f\n", target_dx);
			add_param(param_string, pstring);
			sprintf(pstring, "target_dy=%.1f\n", target_dy);
			add_param(param_string, pstring);
			sprintf(pstring, "avgstk=%e", avgstk);
			param_string_len = add_param(param_string, pstring);
	
			srf2stoch(param_string_len, param_string, &srf, &sfile, debug);
	
			sprintf(slipfile, "%s.slip", infile);
		if (debug) write_slipfile(sfile, slipfile);
		//Need to do some rounding
		int i;
		//Note that casting to int is not the same as floor for negative #s
		//So we are using floorf instead
		for(i=0; i<sfile.nseg; i++) {
			//sfile.elon[i] = ((int)(10000.0*sfile.elon[i]+0.5))/10000.0;
			sfile.elon[i] = (float)(floor(10000.0*sfile.elon[i]+0.5)/10000.0);
			//sfile.elat[i] = ((int)(10000.0*sfile.elat[i]+0.5))/10000.0;
			sfile.elat[i] = (float)(floor(10000.0*sfile.elat[i]+0.5)/10000.0);
		}
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
		int variation_seed = 0;
		if (mode==RUP_GEOM_MODE) {
			//See if we have the seed already
			if (rup_vars[rv].seed==-1) {
				rg_stats_t stats;
				variation_seed = rupgen_get_variation_seed(rup_geom_file, &stats, RUPGEN_UNIFORM_HYPO, header.source_id, header.rupture_id, rv, 0);
			} else {
				variation_seed = rup_vars[rv].seed;
			}
		} else if (mode==SRF_MODE) {
			variation_seed = srf_seed;
		}

		hfsim(seis, stat, slon, slat, local_vmod, fp_out, pga_fp_out, vs30, &header, modelrot, &sfile, num_comps, do_site_response, vref, vpga, variation_seed, rup_vars[rv].rvfrac, debug);

		free(sfile.sp);
       	free(sfile.tr);
       	free(sfile.ti);
		for (i=0; i<3; i++) {
			free(seis[i]);
		}
		free(seis);
		free_srf_ptrs(&srf);
	}
	fflush(fp_out);
	fclose(fp_out);

	fflush(pga_fp_out);
	fclose(pga_fp_out);

	free(rup_vars);

	for (i=0; i<MAX_PARAMS; i++) {
		free(param_string[i]);
	}
	free(param_string);
}

void parse_rup_vars(char* rup_var_string, int num_rup_vars, struct rupture_variation* rup_vars, int rvfrac_seed) {
        int i;
        char* tok, *inner_tok;
        char* outer_save, *inner_save;
		if (rvfrac_seed==0) {
			//rup_var_string is in form (<rv_id>,<slip_id>,<hypo_id>);(....)
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
                //Values not specified
                rup_vars[i].rvfrac = -1.0;
				rup_vars[i].seed = -1;
				//Advance outer
                tok = strtok_r(NULL, "(;)", &outer_save);
	        }
		} else {
			//rup_var_string is in form (<rv_id>,<slip_id>,<hypo_id>,<rvfrac>,<seed>);(....)
			tok = strtok_r(rup_var_string, "(;)", &outer_save);
			for (i=0; i<num_rup_vars; i++) {
				//tok points to <rv_id>,<slip_id>,<hypo_id>,<rvfrac>,<seed>
				inner_tok = strtok_r(tok, ",", &inner_save);
				//parse rv
				rup_vars[i].rup_var_id = atoi(inner_tok);
				inner_tok = strtok_r(NULL, ",", &inner_save);
				//slip
				rup_vars[i].slip_id = atoi(inner_tok);
				inner_tok = strtok_r(NULL, ",", &inner_save);
				//hypo
				rup_vars[i].hypo_id = atoi(inner_tok);
				inner_tok = strtok_r(NULL, ",", &inner_save);
				//rvfrac
				rup_vars[i].rvfrac = atof(inner_tok);
				inner_tok = strtok_r(NULL, ",", &inner_save);
				//seed
				rup_vars[i].seed = atoi(inner_tok);
				//Advance outer
				tok = strtok_r(NULL, "(;)", &outer_save);
			}
		}
tok = strtok_r(NULL, "(;)", &outer_save);}
