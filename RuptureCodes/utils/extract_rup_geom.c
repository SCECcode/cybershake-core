#include <stdio.h>
#include <stdlib.h>

#include "srf_structure.h"

void read_srf(struct standrupformat *srf,char *file,int bflag);

int main(int argc, char** argv) {
	if (argc<6) {
		printf("This code extracts the rupture geometry from an SRF file.\n");
		printf("Usage: %s <text SRF file> <probability> <magnitude> <grid spacing (km)> <output rupture geometry file>\n", argv[0]);
		exit(1);
	}

	char* srf_file = argv[1];
	float probability = atof(argv[2]);
	float magnitude = atof(argv[3]);
	float grid_spacing = atof(argv[4]);
	char* output_file = argv[5];
	//Text format
	int bflag = 0;

	FILE* fp_out = fopen(output_file, "w");

	struct standrupformat srf;
	read_srf(&srf, srf_file, bflag);

	//Single-segment
	if (srf.srf_prect.nseg==1) {

		fprintf(fp_out, "Probability = %e\n", probability);
		fprintf(fp_out, "Magnitude = %.2f\n", magnitude);
		fprintf(fp_out, "GridSpacing = %0.2f\n", grid_spacing);
		fprintf(fp_out, "NumRows = %d\n", srf.srf_prect.prectseg[0].ndip);
		fprintf(fp_out, "NumCols = %d\n", srf.srf_prect.prectseg[0].nstk);
		fprintf(fp_out, "#   Lat         Lon         Depth      Rake    Dip     Strike\n");
		int i;
		for(i=0; i<srf.srf_apnts.np; i++) {
			struct srf_apointvalues apv = srf.srf_apnts.apntvals[i];
			fprintf(fp_out, "%.6f    %.6f    %.7f    %.1f    %.1f    %.5f\n", apv.lat, apv.lon, apv.dep, apv.rake, apv.dip, apv.stk);
		}
	} else {
		long long num_points = 0;
		int i, j;
		for (i=0; i<srf.srf_prect.nseg; i++) {
			num_points += srf.np_seg[i];
		}
		//Multi-segment
        fprintf(fp_out, "Probability = %e\n", probability);
        fprintf(fp_out, "Magnitude = %.2f\n", magnitude);
        fprintf(fp_out, "GridSpacing = %0.2f\n", grid_spacing);
		fprintf(fp_out, "NumPoints = %ld\n", num_points);
        fprintf(fp_out, "#   Lat         Lon         Depth      Rake    Dip     Strike\n");
		struct srf_apointvalues apv;
		int pt_offset = 0;
		for (i=0; i<srf.srf_prect.nseg; i++) {
			for (j=0; j<srf.np_seg[i]; j++) {
				apv = srf.srf_apnts.apntvals[pt_offset+j];
				fprintf(fp_out, "%.6f    %.6f    %.7f    %.1f    %.1f    %.5f\n", apv.lat, apv.lon, apv.dep, apv.rake, apv.dip, apv.stk);
			}
			pt_offset += srf.np_seg[i];
		}
	}
	fflush(fp_out);
	fclose(fp_out);
}

