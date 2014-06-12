#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* This program checks a full SGT file for NaNs.*/

int main(int argc, char** argv) {
	if (argc<2) {
		printf("%s <full SGT file>\n", argv[0]);
		return 1;
	}
	char* sgt_filename = argv[1];
	
	FILE* fp_in = fopen(sgt_filename, "rb");
	//Read in chunks
	if (fp_in==NULL) {
		fprintf(stderr, "Error opening file %s.\n", sgt_filename);
		exit(1);
	}
	long buffer_size = 1*1024*1024*1024;
	long i;
	float* buffer = malloc(buffer_size);
	long floats_read = fread(buffer, sizeof(float), buffer_size/sizeof(float), fp_in);
	long tot_read = floats_read;
	//Also, check for points where all the values are <1e-20
	int consecZeros = 0;
	while (floats_read>0) {
		for (i=0; i<floats_read; i++) {
			if (isnan(buffer[i])) {
				printf("NaN in %s at byte offset %ld.\n", sgt_filename, (tot_read-floats_read+i)*sizeof(float));
				return 2;
			}
			if (fabs(buffer[i])<1e-20) {
				consecZeros++;
			} else {
				consecZeros = 0;
			}
			if (consecZeros>=12000) {
				printf("Found %d consecutive values with an abs less than 1e-20, aborting, at byte offset %ld, index %ld.\n", consecZeros, (tot_read-floats_read+i)*sizeof(float), tot_read-floats_read+i);
				return 3;
			}
		}
		floats_read = fread(buffer, sizeof(float), buffer_size/sizeof(float), fp_in);
		tot_read += floats_read;
	}
	printf("Read %ld floats.\n", tot_read);
	fclose(fp_in);
	return 0;
}
