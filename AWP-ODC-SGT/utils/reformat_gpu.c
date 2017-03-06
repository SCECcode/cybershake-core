#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <input AWP GPU SGT> <timesteps> <num pts> <output AWP SGT>\n", argv[0]);
		exit(1);
	}
	
	/*
	Reordering from fast component, station, slow timestep
	to fast timestep, component, slow station
	*/

	FILE* fp_in;
	FILE* fp_out;
	int i, j, k;
	long out_offset;
	char* filename_in = argv[1];
	int timesteps = atoi(argv[2]);
	int num_pts = atoi(argv[3]);
	char* filename_out = argv[4];

	float* sgt_in_data = malloc(sizeof(float)*6*(long)timesteps*num_pts);
	float* sgt_out_data = malloc(sizeof(float)*6*(long)timesteps*num_pts);

	fp_in = fopen(filename_in, "r");
	fp_out = fopen(filename_out, "w");

	printf("Reading data.\n");
	fread(sgt_in_data, sizeof(float), 6*(long)timesteps*(long)num_pts, fp_in);
	fclose(fp_in);

	for (i=0; i<timesteps; i++) {
		if (i%100==0) {
			printf("%d of %d timesteps.\n", i, timesteps);
		}
		for (j=0; j<num_pts; j++) {
			for (k=0; k<6; k++) {
				out_offset = j*timesteps*6 + k*timesteps + i;
				sgt_out_data[out_offset] = sgt_in_data[i*num_pts*6 + j*6 + k];
			}
		}
	}

	fwrite(sgt_out_data, sizeof(float), 6*(long)timesteps*(long)num_pts, fp_out);
	fflush(fp_out);
	fclose(fp_out);
	free(sgt_in_data);
	free(sgt_out_data);
}
