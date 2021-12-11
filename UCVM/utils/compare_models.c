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
	float abs_diff;
	float abs_per_diff;
	double tot_abs_diff = 0.0;
	double tot_abs_per_diff = 0.0;
	long max_diff_index = -1;
	long max_per_diff_index = -1;

	long index = 0;
	long i;
	int j;

	//Track histogram, too
        int num_bins = 11;
        float cutoffs[10] = {0.0001, 0.001, 0.01, 0.1, 0.3, 1, 3, 10, 100, 1000};
        long bins[11] = {0};

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
			abs_diff = fabs(diff);
			abs_per_diff = fabs(per_diff);
			if (abs_diff>max_diff) {
				max_diff = abs_diff;
				max_diff_index = index+i;
			}
			if (abs_per_diff>max_per_diff) {
				max_per_diff = abs_per_diff;
				max_per_diff_index = index+i;
			}
			tot_abs_diff += abs_diff;
			tot_abs_per_diff += abs_per_diff;
                        int flag = 0;
                        for (j=0; j<num_bins-1; j++) {
	                        if (abs_per_diff<cutoffs[j]) {
                                        bins[j]++;
                                        flag = 1;
                                        break;
                                }
                        }
                        if (flag==0) {
                                bins[num_bins-1]++;
                        }
		}
		index += num_pts_to_read;
	}
	printf("Average absolute difference: %.f\n", (tot_abs_diff/num_pts));
	printf("Average absolute percentage difference: %f%%\n", (100.0*tot_abs_per_diff/num_pts));
	printf("Largest difference: %.2f, at index %ld.\n", max_diff, max_diff_index);
	printf("Largest percent difference: %.2f%%, at index %ld\n", max_per_diff*100.0, max_per_diff_index);
	printf("Absolute percentage difference histogram:\n");
        printf("   ");
        for (i=0; i<num_bins-1; i++) {
                printf("%.4f   ", cutoffs[i]);
        }
        printf("\n");
        for (i=0; i<num_bins; i++) {
                printf("%ld   ", bins[i]);
        }
        printf("\n ");
        for (i=0; i<num_bins; i++) {
                printf("%.2f%%   ", ((float)bins[i])/num_pts*100.0);
        }
        printf("\n");
        float tot_per = 0.0;
        for (i=0; i<num_bins; i++) {
                tot_per += ((float)bins[i])/num_pts*100.0;
                printf("%.2f%%   ", tot_per);
        }
        printf("\n");

	fclose(fp_in1);
	fclose(fp_in2);	

	free(buf1);
	free(buf2);	
	return 0;
}
