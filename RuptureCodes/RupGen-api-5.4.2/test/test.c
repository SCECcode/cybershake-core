#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "rupgen_api.h"

int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <erf> <src> <rup> <rv start> <rv end>\n", argv[0]);
		exit(1);
	}
	int erf_id = atoi(argv[1]);
	int src = atoi(argv[2]);
	int rup = atoi(argv[3]);
	int rv_start = atoi(argv[4]);
	int rv_end = atoi(argv[5]);	

        set_memcached_server("localhost");
        char* rupture_file = malloc(sizeof(char)*256);
        sprintf(rupture_file, "/gpfs/alpine/proj-shared/geo112/CyberShake/ruptures/Ruptures_erf%d/%d/%d/%d_%d.txt", erf_id, src, rup, src, rup);
        char outfile[256];
        rg_stats_t stats;
        int slip, hypo;
        rupgen_get_num_rv(rupture_file, &stats, RUPGEN_UNIFORM_HYPO);
        printf("%d slips, %d hypos.\n", stats.numslip, stats.numhypo);
        int i, j;
        struct standrupformat srf;
	struct timespec ts_start, ts_end;
        //for (i=0; i<stats.numslip; i++) {
	for (i=rv_start; i<rv_end; i++) {
                for (j=0; j<stats.numhypo; j++) {
			clock_gettime(CLOCK_MONOTONIC, &ts_start);
                        rupgen_genslip(rupture_file, i, j, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05);
			clock_gettime(CLOCK_MONOTONIC, &ts_end);
			sprintf(outfile, "%d_%d_%d_uniform.srf", src, rup, i);
			write_srf(&srf, outfile, 0);
                        free_srf_ptrs(&srf);
			long elapsed_time_ns = 1000000000L*(ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec);
			printf("Elapsed time: %ld.%09ld\n", (elapsed_time_ns/1000000000L), (elapsed_time_ns%1000000000));
                }
        }
	free(rupture_file);
}

