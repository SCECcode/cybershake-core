#include <stdio.h>
#include <stdlib.h>

#include "rupgen_api.h"

const char* RUPTURE_ROOT = "/gpfs/alpine/proj-shared/geo112/CyberShake/ruptures";

int main(int argc, char** argv) {
	if (argc<5) {
		printf("Usage: %s <erf id> <src id> <rup id> <output file>\n", argv[0]);
		exit(1);
	}
	set_memcached_server("localhost");
        char rupture_file[512];
        int erf = atoi(argv[1]);
	int src = atoi(argv[2]);
        int rup = atoi(argv[3]);
	char* output_file = argv[4];
        sprintf(rupture_file, "%s/Ruptures_erf%d/%d/%d/%d_%d.txt", RUPTURE_ROOT, erf, src, rup, src, rup);
        rg_stats_t stats;
        int slip, hypo;
        float dt = 0.05;
        rupgen_get_num_rv(rupture_file, &stats, RUPGEN_UNIFORM_HYPO);
        struct standrupformat srf;
        int i, j, ip, id;
	float hlon, hlat, hdep;
	float tmin;
	struct srf_apointvalues *apval_ptr1;
	FILE* fp_out = fopen(output_file, "w");
        for (i=0; i<stats.numslip; i++) {
                for (j=0; j<stats.numhypo; j++) {
			rupgen_genslip(rupture_file, i, j, &stats, &srf, RUPGEN_UNIFORM_HYPO, dt);
			id = i*stats.numhypo + j;
			tmin = 1.0e+15;
			apval_ptr1 = srf.srf_apnts.apntvals;
			for (ip=0; ip<srf.srf_apnts.np; ip++) {
				if(apval_ptr1[ip].tinit < tmin) {
					hlon = apval_ptr1[ip].lon;
					hlat = apval_ptr1[ip].lat;
					hdep = apval_ptr1[ip].dep;
					tmin = apval_ptr1[ip].tinit;
				}
			}
			fprintf(fp_out, "%d\t%.5f\t%.5f\t%.5f\n",id,hlon,hlat,hdep);
                        free_srf_ptrs(&srf);
                }
        }
	fflush(fp_out);
	fclose(fp_out);
}

