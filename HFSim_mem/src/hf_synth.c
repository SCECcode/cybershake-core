#include "structure.h"
#include "include.h"
#include "function.h"
#include "defs.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

const int X_COMP_FLAG = 1;
const int Y_COMP_FLAG = 2;
const int Z_COMP_FLAG = 4;

void parse_rup_vars(char* rup_var_string, int num_rup_vars, struct rupture_variation* rup_vars);

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
				fprintf(fp_out, " %5.0f", sfile.sp[j*NQ*LV + k*LV + i]);
			}
			fprintf(fp_out, "\n");
		}
                for (j=0; j<sfile.ny[i]; j++) {
                        for (k=0; k<sfile.nx[i]; k++) {
                                fprintf(fp_out, " %5.2f", sfile.tr[j*NQ*LV + k*LV + i]);
                        }
			fprintf(fp_out, "\n");
                }
                for (j=0; j<sfile.ny[i]; j++) {
                        for (k=0; k<sfile.nx[i]; k++) {
                                fprintf(fp_out, " %5.2f", sfile.ti[j*NQ*LV + k*LV + i]);
                        }
			fprintf(fp_out, "\n");
                }
	}
	fclose(fp_out);
}

int main(int argc, char** argv) {
	//for srf2stoch
	char infile[256] = {'\0'};
	char outfile[256];
	char slipfile[256];
	float dx;
	float dy;
	float avgstk = -1.0e+15;
	int inbin = 0;
	//for hfsim
	char stat[12];
	float slon, slat;
	char local_vmod[256];
	float vs30 = -1.0;
	float tlen = 102.4;
	float dt = 0.025;
	float modelrot = -55;

	char rup_geom_file[256];
	int slip_id;
	int hypo_id;

	int source_id = -1;
	int rupture_id = -1;
	int rup_var_id = -1;
	int det_max_freq = -1.0;
	int stoch_max_freq = 10.0;

	int do_site_response = 1;

	int num_rup_vars = -1;
	char rup_var_string[16384];

	int debug = 0;

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
	mstpar("dx", "f", &dx);
	mstpar("dy", "f", &dy);
	getpar("inbin", "d", &inbin);
	getpar("avgstk", "f", &avgstk);
	//for hfsim
        mstpar("stat", "s", stat);
        mstpar("slon", "f", &slon);
        mstpar("slat", "f", &slat);
        mstpar("vmod", "s", local_vmod);
        mstpar("outfile", "s", outfile);
        getpar("vs30", "f", &vs30);
        getpar("tlen", "f", &tlen);
        getpar("dt", "f", &dt);
        getpar("modelrot", "f", &modelrot);
	getpar("do_site_response", "d", &do_site_response);
	//Things to populate CyberShake header
        mstpar("source_id", "d", &source_id);
        mstpar("rupture_id", "d", &rupture_id);
        getpar("det_max_freq", "f", &det_max_freq);
        getpar("stoch_max_freq", "f", &stoch_max_freq);

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
		mstpar("rup_vars", "s", rup_var_string);
		rup_vars = check_malloc(sizeof(struct rupture_variation)*num_rup_vars);
		parse_rup_vars(rup_var_string, num_rup_vars, rup_vars);
	}
	endpar();

	FILE* fp_out = fopen(outfile, "wb");

	int rv = 0;
	//Generate them sequentially
	for (rv=0; rv<num_rup_vars; rv++) {
		float seis[3][mmv];
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
		header.comps = X_COMP_FLAG | Y_COMP_FLAG | Z_COMP_FLAG;

	        struct slipfile sfile;
	        sfile.sp = check_malloc(NP*NQ*LV*sizeof(float));
	        sfile.tr = check_malloc(NP*NQ*LV*sizeof(float));
	        sfile.ti = check_malloc(NP*NQ*LV*sizeof(float));

		srf2stoch(rup_geom_file, rup_vars[rv].slip_id, rup_vars[rv].hypo_id, infile, dx, dy, inbin, avgstk, &sfile, dt, debug);
		sprintf(slipfile, "%s.slip", infile);
		if (debug) write_slipfile(sfile, slipfile);
		//Need to int the contents of sfile.sp
		//Round the others to 2 digits
		int i;
		for(i=0; i<sfile.nseg; i++) {
			sfile.ravg[i] = (int)(sfile.ravg[i]+0.5);
			sfile.dtop[i] = ((int)(100.0*sfile.dtop[i]+0.5))/100.0;
			sfile.shypo[i] = ((int)(100.0*sfile.shypo[i]+0.5))/100.0;
			sfile.dhypo[i] = ((int)(100.0*sfile.dhypo[i]+0.5))/100.0;
		}
		for(i=0; i<NP*NQ*LV; i++) {
			sfile.sp[i] = (int)(sfile.sp[i]+0.5);
			sfile.tr[i] = ((int)(100.0*sfile.tr[i]+0.5))/100.0;
			sfile.ti[i] = ((int)(100.0*sfile.ti[i]+0.5))/100.0;
		}
		for(i=0; i<LV; i++) {
			sfile.dx[i]=((int)(100.0*sfile.dx[i]+0.5))/100.0;
			sfile.dy[i]=((int)(100.0*sfile.dy[i]+0.5))/100.0;
		}
		//This is the default for LA Basin; other velocity models might be different
		sfile.qfexp = 0.6;
		hfsim(seis, stat, slon, slat, local_vmod, outfile, vs30, &header, modelrot, &sfile, do_site_response, debug);

		//Write to file
	        fwrite(&header, sizeof(struct seisheader), 1, fp_out);
	        for (i=0; i<3; i++) {
	                fwrite(seis[i], sizeof(float), header.nt, fp_out);
	        }
		fflush(fp_out);

		free(sfile.sp);
        	free(sfile.tr);
        	free(sfile.ti);
	}
	fclose(fp_out);
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
