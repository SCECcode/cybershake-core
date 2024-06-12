#include <stdio.h>
#include <stdlib.h>
#include "math.h"

//Code to compare two raw (AWP-format) SGTs

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <reference file> <test file> <num pts> <timesteps>\n", argv[0]);
		exit(1);
	}

	char* ref_file = argv[1];
	char* test_file = argv[2];
	long long num_pts = atol(argv[3]);
	int timesteps = atoi(argv[4]);
	printf("Reading %lld pts, %d timesteps\n", num_pts, timesteps);

	FILE* ref_fp = fopen(ref_file, "rb");
	FILE* test_fp = fopen(test_file, "rb");

	int i, j, k;
	//Read in ~2 GB at a time
	int pts_per_read = 2000000000/(6*timesteps*sizeof(float));
	long long tot_pts_read = 0;
	int pts_to_read;
	float* ref_buf = malloc(pts_per_read * 6 * timesteps * sizeof(float));
	float* test_buf = malloc(pts_per_read * 6 * timesteps * sizeof(float));
	float diff, percent_diff, abs_percent_diff;
	float max_diff = 0;
	float max_percent_diff = 0;
	long long max_diff_index = -1L;
	long long max_percent_diff_index = -1L;
	double tot_diff;
	double tot_percent_diff;
	double tot_abs_per_diff;
	long long percent_num_vals = 0;
	float PERCENT_CUTOFF = 1e-10;
	//histogram of abs percent differences
	int num_bins = 11;
	float cutoffs[10] = {0.0001, 0.001, 0.01, 0.1, 0.3, 1, 3, 10, 100, 1000};
	long bins[11] = {0};
	printf("Calculating percent differences for values larger than %e.\n", PERCENT_CUTOFF);
	while (tot_pts_read<num_pts) {
		pts_to_read = pts_per_read;
		if (tot_pts_read + pts_to_read > num_pts) {
			pts_to_read = num_pts - tot_pts_read;
		}
		printf("%lld of %lld pts read.\n", tot_pts_read, num_pts);
		fread(ref_buf, sizeof(float), pts_to_read*6*timesteps, ref_fp);
		fread(test_buf, sizeof(float), pts_to_read*6*timesteps, test_fp);
		for (i=0; i<pts_to_read*6*timesteps; i++) {
			diff = test_buf[i]-ref_buf[i];
			tot_diff += diff;
			if (diff>max_diff) {
				max_diff = diff;
				max_diff_index = (tot_pts_read*6*timesteps) + i;
				printf("New max diff: index %d, %e vs %e\n", i, ref_buf[i], test_buf[i]);
			}
			if (ref_buf[i]>PERCENT_CUTOFF) {
				percent_diff = diff/ref_buf[i]*100.0;
				abs_percent_diff = fabs(percent_diff);
				tot_percent_diff += percent_diff;
				tot_abs_per_diff += abs_percent_diff;
				percent_num_vals++;
				if (abs_percent_diff>max_percent_diff) {
					max_percent_diff = fabs(percent_diff);
					max_percent_diff_index = (tot_pts_read*6*timesteps) + i;
					printf("New max percent diff: index %d, %e vs %e\n", i, ref_buf[i], test_buf[i]);
				}
				int flag = 0;
				for (j=0; j<num_bins-1; j++) {
					if (abs_percent_diff<cutoffs[j]) {
						bins[j]++;
						flag = 1;
						break;
					}
				}
				if (flag==0) {
					bins[num_bins-1]++;
				}
			}
		}
		tot_pts_read += pts_to_read;
	}
	float avg_diff = tot_diff/((long long)num_pts*6*timesteps);
	float avg_percent_diff = tot_percent_diff/percent_num_vals;
	float avg_abs_per_diff = tot_abs_per_diff/percent_num_vals;
	printf("Average diff = %e.\n", avg_diff);
	printf("Using percent cutoff of %f, average percent diff = %f%%, average absolute percent diff = %f%%\n", PERCENT_CUTOFF, avg_percent_diff,avg_abs_per_diff);
	printf("Largest diff of %f at float index %lld.\n", max_diff, max_diff_index);
	printf("Largest percent diff of %f%% at float index %lld.\n", max_percent_diff, max_percent_diff_index);
	printf("Absolute percentage difference histogram:\n");
	printf("Bins:         ");
	for (i=0; i<num_bins-1; i++) {
		printf("%.4f      ", cutoffs[i]);
	}
	printf("\n");
	printf("Values:");
	for (i=0; i<num_bins; i++) {
		printf("%ld   ", bins[i]);
	}
	printf("\n");
	printf("%% in bin: ");
	for (i=0; i<num_bins; i++) {
		printf("%.2f%%     ", ((float)bins[i])/percent_num_vals*100.0);
	}
	printf("\n");
	printf("CDist:    ");
	float tot_per = 0.0;
	for (i=0; i<num_bins; i++) {
		tot_per += ((float)bins[i])/percent_num_vals*100.0;
		printf("%.2f%%     ", tot_per);
	}
	printf("\n");
	free(ref_buf);
	free(test_buf);
	return 0;
}
