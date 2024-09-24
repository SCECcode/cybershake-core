#include <stdio.h>
#include <stdlib.h>

#include "srf_structure.h"

struct pointsource *read_gsfpars(char *file,struct pointsource *psrc,struct generic_slip *gslip,float *dx,float *dy,float *dtop,float *dip);

int main(int argc, char** argv) {
	if (argc<8) {
		printf("This code extracts the rupture geometry from a GSF file, before fault roughness is added.\n");
		printf("Usage: %s <text GSF file> <probability> <magnitude> <grid spacing (km)> <num along dip> <num along strike> <output rupture geometry file>\n", argv[0]);
		exit(1);
	}

	char* gsf_file = argv[1];
	float probability = atof(argv[2]);
	float magnitude = atof(argv[3]);
	float grid_spacing = atof(argv[4]);
	int ndip = atoi(argv[5]);
	int nstk = atoi(argv[6]);
	char* output_file = argv[7];
	//Text format
	int bflag = 0;

	FILE* fp_out = fopen(output_file, "w");

	//struct pointsource *read_gsfpars(char *file,struct pointsource *psrc,struct generic_slip *gslip,float *dx,float *dy,float *dtop,float *dip)
	struct pointsource* psrc_orig = NULL;
	struct generic_slip gslip;
	gslip.np = -1;
	gslip.spar = NULL;
	float dstk, ddip, dtop, avgstk, avgdip;
	psrc_orig = read_gsfpars(gsf_file, psrc_orig, &gslip, &dstk, &ddip, &dtop, &avgdip);

	fprintf(fp_out, "Probability = %e\n", probability);
	fprintf(fp_out, "Magnitude = %.2f\n", magnitude);
	fprintf(fp_out, "GridSpacing = %0.2f\n", grid_spacing);
	fprintf(fp_out, "NumRows = %d\n", ndip);
	fprintf(fp_out, "NumCols = %d\n", nstk);
	fprintf(fp_out, "#   Lat         Lon         Depth      Rake    Dip     Strike\n");
	int i;
	for(i=0; i<gslip.np; i++) {
		struct slippars sp = gslip.spar[i];
		fprintf(fp_out, "%.6f    %.6f    %.7f    %.1f    %.1f    %.5f\n", sp.lat, sp.lon, sp.dep, sp.rake, sp.dip, sp.stk);
	}
	fflush(fp_out);
	fclose(fp_out);
}

