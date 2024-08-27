#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "rupgen_api.h"

int main(int argc, char** argv) {
	if (argc<7) {
		printf("Usage: %s <erf> <src> <rup> <rv> <rvfrac> <seed>\n", argv[0]);
		exit(1);
	}
	int erf_id = atoi(argv[1]);
	int src = atoi(argv[2]);
	int rup = atoi(argv[3]);
	int rv = atoi(argv[4]);
	float rvfrac = atof(argv[5]);
	int seed = atoi(argv[6]);	

        set_memcached_server("localhost");
        char* rupture_file = malloc(sizeof(char)*256);
        sprintf(rupture_file, "/work2/00349/scottcal/frontera/CyberShake/ruptures/Ruptures_erf%d/%d/%d/%d_%d.txt", erf_id, src, rup, src, rup);
        char outfile[256];
        rg_stats_t stats;
        int slip, hypo;
        struct standrupformat srf;
	struct timespec ts_start, ts_end;
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	char params[256];
	sprintf(params, "seed=%d rvfrac=%f use_unmodified_seed=1", seed, rvfrac);
    rupgen_genslip_with_params(rupture_file, rv, 0, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05, params);
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	sprintf(outfile, "%d_%d_%d_uniform.srf", src, rup, rv);
	_write_srf(&srf, outfile, 0);
    free_srf_ptrs(&srf);
	long elapsed_time_ns = 1000000000L*(ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec);
	printf("Elapsed time: %ld.%09ld\n", (elapsed_time_ns/1000000000L), (elapsed_time_ns%1000000000));
	free(rupture_file);
}

