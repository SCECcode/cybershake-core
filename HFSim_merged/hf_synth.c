#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951


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
	char infile[256];
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
	float vs30 = -1;
	float tlen = 102.4;
	float dt = 0.025;
	float modelrot = -55;

	memset(stat, ' ', 12);
	memset(local_vmod, ' ', 256);
	memset(outfile, ' ', 256);	

	setpar(argc, argv);
	//for srf2stoch
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
        mstpar("vs30", "f", &vs30);
        getpar("tlen", "f", &tlen);
        getpar("dt", "f", &dt);
        getpar("modelrot", "f", &modelrot);
	endpar();


	struct slipfile sfile;
	sfile.sp = check_malloc(NP*NQ*LV*sizeof(float));
        sfile.tr = check_malloc(NP*NQ*LV*sizeof(float));
        sfile.ti = check_malloc(NP*NQ*LV*sizeof(float));

	srf2stoch(infile, dx, dy, inbin, avgstk, &sfile);
	sprintf(slipfile, "%s.slip", infile);
	write_slipfile(sfile, slipfile);
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
	hfsim(stat, slon, slat, local_vmod, outfile, vs30, tlen, dt, modelrot, &sfile);

	free(sfile.sp);
        free(sfile.tr);
        free(sfile.ti);
}
