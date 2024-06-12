#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "rupgen_api.h"

int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <erf> <src> <rup> <rv> <seed> <rvfrac>\n", argv[0]);
		exit(1);
	}
	int erf_id = atoi(argv[1]);
	int src = atoi(argv[2]);
	int rup = atoi(argv[3]);
	int rv = atoi(argv[4]);
	int seed = atoi(argv[5]);
	float rvfrac = atof(argv[6]);
	/*int erf_id=60;
	int src=0;
	int rup=0;
	int rv=0;
	int seed=2379646;*/

    set_memcached_server("localhost");
    char* rupture_file = malloc(sizeof(char)*256);
    sprintf(rupture_file, "/gpfs/alpine/proj-shared/geo112/CyberShake/ruptures/Ruptures_erf%d/%d/%d/%d_%d.txt", erf_id, src, rup, src, rup);
    //sprintf(rupture_file, "m6.73-0.10x0.10.gsf");
    char outfile[256];
    struct standrupformat srf;
	struct rg_stats_t stats;
	char params[2048];
	sprintf(params, "seed=%d", seed);
	sprintf(params, "%s rvfrac=%f", params, rvfrac);
	sprintf(params, "%s use_unmodified_seed=1", params);
	/*sprintf(params, "%s read_erf=0 read_gsf=1 calc_shypo=0 shypo=6.00 dhypo=19.40", params);
	sprintf(params, "%s mag=6.73 nstk=200 ndip=270 ns=1 nh=1 rvfrac=0.811046", params);
	sprintf(params, "%s kmodel=2 slip_sigma=0.75 circular_average=0", params);
	sprintf(params, "%s modified_corners=0 shal_vrup_dep=6.500000 shal_vrup_deprange=1.500000", params);
	sprintf(params, "%s shal_vrup_dep=6.500000 shal_vrup_deprange=1.500000 shal_vrup=0.600000 side_taper=0.02 bot_taper=0.0", params);
	sprintf(params, "%s top_taper=0.0 dt=0.100000 risetime_coef=1.600000 plane_header=1 risetimefac=2.000000", params);
	sprintf(params, "%s risetimedep=6.500000 risetimedep_range=1.500000 rt_scalefac=1.000000 slip_water_level=-1.000000", params);
	sprintf(params, "%s deep_risetimedep=17.500000 deep_risetimedep_range=2.500000 deep_risetimefac=2.000000", params);
	sprintf(params, "%s flen_max=-1.000000 rupture_delay=0.000000 moment_fraction=-1.000000 srf_version=2.0 rake_sigma=15.0", params);
	sprintf(params, "%s fdrup_time=1 deep_vrup=0.6 use_gaus=1 alpha_rough=0.01 lambda_min=0.08", params);
	sprintf(params, "%s tsfac_coef=1.1 tsfac1_sigma=1.0 tsfac1_scor=0.8 rtime1_sigma=0.750000 rtime1_scor=0.5 tsfac_bzero=-0.1 tsfac_slope=-0.5", params);
	sprintf(params, "%s velfile=/gpfs/alpine/proj-shared/geo112/CyberShake/software/bbp/bbp_data/indata/4277838/nr02-vs500.fk1d", params);*/
	rupgen_genslip_with_params(rupture_file, rv, 0, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05, params);
	sprintf(outfile, "%d_%d_%d_seed%d_uniform.srf", src, rup, rv, seed);
	write_srf(&srf, outfile, 0);
    free_srf_ptrs(&srf);
	free(rupture_file);
}

