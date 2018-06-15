#include <stdio.h>
#include <stdlib.h>

#include "math.h"

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <velocity model 1> <velocity model 2> <num pts>\n", argv[0]);
		exit(1);
	}

	char* filename1 = argv[1];
	char* filename2 = argv[2];
	long num_pts = atol(argv[3]);

	FILE* fp_in1 = fopen(filename1, "rb");
	FILE* fp_in2 = fopen(filename2, "rb");

	int BUF_SIZE = 1000000000;
	int PTS_PER_BUF = BUF_SIZE/sizeof(float);
	float* buf1 = malloc(BUF_SIZE);
	float* buf2 = malloc(BUF_SIZE);

	float max_diff = 0.0;
	float max_per_diff = 0.0;
	float diff;
	float per_diff;
	long max_diff_index = -1;
	long max_per_diff_index = -1;

	long index = 0;
	long i;
	while (index<num_pts) {
		printf("Index %ldK of %ldK.\n", index/1000, num_pts/1000);
		long num_pts_to_read = PTS_PER_BUF;
		if (index + num_pts_to_read > num_pts) {
			num_pts_to_read = num_pts - index;
		}
		fread(buf1, num_pts_to_read, sizeof(float), fp_in1);
                fread(buf2, num_pts_to_read, sizeof(float), fp_in2);
		for (i=0; i<num_pts_to_read; i++) {
			diff = buf2[i] - buf1[i];
			per_diff = (buf2[i]-buf1[i])/buf1[i];
			if (fabs(diff)>max_diff) {
				max_diff = diff;
				max_diff_index = index+i;
			}
			if (fabs(per_diff)>max_per_diff) {
				max_per_diff = per_diff;
				max_per_diff_index = index+i;
			}
		}
		index += num_pts_to_read;
	}
	printf("Largest difference: %.2f, at index %ld.\n", max_diff, max_diff_index);
	printf("Largest percent difference: %.2f%%, at index %ld\n", max_per_diff*100.0, max_per_diff_index);

	fclose(fp_in1);
	fclose(fp_in2);	

	free(buf1);
	free(buf2);	
	return 0;
}
