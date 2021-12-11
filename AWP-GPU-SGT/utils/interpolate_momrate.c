#include <stdio.h>
#include <stdlib.h>

#include "math.h"

int main(int argc, char** argv) {
	if (argc<7) {
		printf("Usage: %s <momrate file> <npts> <nt> <dt> <interp factor> <momrate out>\n", argv[0]);
		exit(1);
	}

	char* momrate_file = argv[1];
	int npts = atoi(argv[2]);
	int nt = atoi(argv[3]);
	float dt = atof(argv[4]);
	int interp_factor = atoi(argv[5]);
	char* momrate_out = argv[6];
	float new_dt = dt/interp_factor;
	int new_nt = nt*interp_factor;
	printf("new_dt=%f, new_nt=%d\n", new_dt, new_nt);

	FILE* fp_in = fopen(momrate_file, "rb");
	FILE* fp_out = fopen(momrate_out, "wb");
	int i, j, k;
	int coords[3];
	//Buffers are fast component, slow time
	float* in_buffer = malloc(sizeof(float)*6*nt);
	float* out_buffer = malloc(sizeof(float)*6*new_nt);
	int time_index[2];
	int out_buffer_index, in_buffer_index;
	float m, y[2], interp_value;
	float ZERO_CUTOFF = 1e-30;
	for (i=0; i<npts; i++) {
	//for (i=0; i<1; i++) {
		if (i%1000==0) {
			printf("Processing point %d of %d.\n", i, npts);
		}
		fread(coords, sizeof(float), 3, fp_in);
		fwrite(coords, sizeof(float), 3, fp_out);
		fread(in_buffer, sizeof(float), 6*nt, fp_in);
		for (j=0; j<6; j++) {
			for (k=0; k<new_nt; k++) {
				out_buffer_index = k*6 + j;
				if (k%interp_factor==0) {
					//Falls on a old point, keep
					in_buffer_index = (k/interp_factor)*6 + j;
					out_buffer[out_buffer_index] = in_buffer[in_buffer_index];
					continue;
				}
				//Use linear interpolation
				//printf("k=%d, interp_factor=%d\n", k, interp_factor);
				time_index[0] = (int)(k/interp_factor);
				time_index[1] = time_index[0]+1;
				y[0] = in_buffer[time_index[0]*6 + j];
				y[1] = in_buffer[time_index[1]*6 + j];
				//Don't interpolate if this is before the first non-zero value, or after the last non-zero value.
				if (y[0]==0.0 || y[1]==0.0) {
					out_buffer[out_buffer_index] = 0.0;
					continue;
				}
				m = (y[1]-y[0])/(dt*(time_index[1]-time_index[0]));
				interp_value = m*(k*new_dt-time_index[0]*dt)+y[0];
				/*if (k==2121) {
					printf("time_index[0]=%d, time_index[1]=%d, y[0]=%e, y[1]=%e, m=%e, interp_value=%e, out_buffer_index=%d\n", time_index[0], time_index[1], y[0], y[1], m, interp_value, out_buffer_index);
				}*/
				out_buffer[out_buffer_index] = interp_value;
			}
		}
		fwrite(out_buffer, sizeof(float), 6*new_nt, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);
}
