#include <stdio.h>
#include <stdlib.h>

#include "rupgen_api.h"

const char* RUPTURE_ROOT = "/work2/00349/scottcal/frontera/CyberShake/ruptures";

int main(int argc, char** argv) {
	set_memcached_server("localhost");
	if (argc<4) {
		printf("Usage: %s <erf id> <src id> <rup id>\n", argv[0]);
		exit(1);
	}
	int erf = atoi(argv[1]);
	int src = atoi(argv[2]);
	int rup = atoi(argv[3]);
	char rupture_file[512];
	sprintf(rupture_file, "%s/Ruptures_erf%d/%d/%d/%d_%d.txt", RUPTURE_ROOT, erf, src, rup, src, rup);
	rg_stats_t stats;
	int slip, hypo;
	rupgen_get_num_rv(rupture_file, &stats, RUPGEN_UNIFORM_HYPO);
	printf("%d\n", stats.numslip*stats.numhypo);
	return 0;
}
