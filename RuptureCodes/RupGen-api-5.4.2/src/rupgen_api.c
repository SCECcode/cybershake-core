#include <stdlib.h>
#include <stdio.h>
#include <libgen.h>
#include "string.h"
#include "math.h"

#include "rupgen_defs.h"
#include "rupgen_api.h"
#include "GenRandV5.0/function.h"

#ifdef _USE_MEMCACHED
char mc_server[128] = {'\0'};

int set_memcached_server(char* memcached_server) {
	strcpy(mc_server, memcached_server);
	return 0;
}
#endif

int rupgen_get_num_rv(char* rup_geom_file, rg_stats_t *stats, int hypo_dist) {
	char params[1] = "\0";
	return rupgen_get_num_rv_with_params(rup_geom_file, stats, hypo_dist, params);
}

int rupgen_get_num_rv_with_params(char* rup_geom_file, rg_stats_t *stats, int hypo_dist, char* params) {
	int writeout = 0;

	int j, rgargc;
	char **rgargv = NULL;

#ifdef _USE_MEMCACHED
	if (strcmp(mc_server, "")==0) {
		fprintf(stderr, "Error:  need to set memcached server with set_memcached_server()\n");
		exit(1);
	}
#endif

	/* Create pseudo-program args */
        int MAX_ARGS = 30;
        int index = 0;
        rgargc = MAX_ARGS;
        rgargv = malloc(rgargc * sizeof(char*));
        for (j = 0; j < MAX_ARGS; j++) {
                rgargv[j] = malloc(MAX_STR_ARGV);
                memset(rgargv[j], 0, MAX_STR_ARGV);
        }
	sprintf(rgargv[0], "%s", "rupgen");
	sprintf(rgargv[1], "infile=%s", rup_geom_file);
        if (hypo_dist==RUPGEN_RANDOM_HYPO) {
                sprintf(rgargv[2], "uniformgrid_hypo=0");
                sprintf(rgargv[3], "random_hypo=1");
        } else if (hypo_dist==RUPGEN_UNIFORM_HYPO) {
                sprintf(rgargv[2], "uniformgrid_hypo=1");
                sprintf(rgargv[3], "random_hypo=0");
        } else {
                fprintf(stderr, "Error, did not specify a valid hypocenter location distribution, aborting.\n");
                exit(2);
        }
        //Now, add args from passed-in parameters
        index = 4;
        char* tok = strtok(params, " ");
        while (tok!=NULL) {
        	sprintf(rgargv[index], tok);
                index++;
                tok = strtok(NULL, " ");
        } 

	/* Run rupture generator */
	struct standrupformat srf;
#ifdef _USE_MEMCACHED
	mc_genslip(rgargc, rgargv, stats, &srf, GET_STATS, mc_server);
#else
	genslip(rgargc, rgargv, stats, &srf, GET_STATS);
#endif
  
	/* Free pseudo-program args */
	for (j = 0; j < MAX_ARGS; j++) {
	    free(rgargv[j]);
	}

  	free(rgargv);
	return(0);

}



//Includes string "params" containing getpar-friendly key=value space-delimited pairs, passed directly through to genslip args
int rupgen_genslip_with_params(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf, int hypo_dist, float dt, char* params) {
        int write_srf = 1;

        int i, j, rgargc;
        char **rgargv = NULL;

#ifdef _USE_MEMCACHED
        if (strcmp(mc_server, "")==0) {
                fprintf(stderr, "Error:  need to set memcached server with set_memcached_server()\n");
                exit(1);
        }
#endif

        char* srf_out_file = malloc(strlen(rup_geom_file)+20);
        sprintf(srf_out_file, "%s_s%d_h%d.srf", basename(rup_geom_file), slip, hypo);

        /* Create pseudo-program args */
	int MAX_ARGS = 100;
	int index = 0;
        rgargc = MAX_ARGS;
        rgargv = malloc(rgargc * sizeof(char*));
        for (j = 0; j < MAX_ARGS; j++) {
                rgargv[j] = malloc(MAX_STR_ARGV);
                memset(rgargv[j], 0, MAX_STR_ARGV);
        }
        sprintf(rgargv[0], "%s", "rupgen");
        sprintf(rgargv[1], "infile=%s", rup_geom_file);
        sprintf(rgargv[2], "doslip=%d", slip);
        sprintf(rgargv[3], "dohypo=%d", hypo);
        sprintf(rgargv[4], "write_srf=%d", write_srf);
        sprintf(rgargv[5], "outfile=%s", srf_out_file);
        sprintf(rgargv[6], "dt=%f", dt);
        if (hypo_dist==RUPGEN_RANDOM_HYPO) {
                sprintf(rgargv[7], "uniformgrid_hypo=0");
                sprintf(rgargv[8], "random_hypo=1");
        } else if (hypo_dist==RUPGEN_UNIFORM_HYPO) {
                sprintf(rgargv[7], "uniformgrid_hypo=1");
                sprintf(rgargv[8], "random_hypo=0");
        } else if (hypo_dist==RUPGEN_PROVIDE_HYPO) { 
				sprintf(rgargv[7], "uniformgrid_hypo=0");
				sprintf(rgargv[8], "random_hypo=0");
		} else {
                fprintf(stderr, "Error, did not specify a valid hypocenter location distribution, aborting.\n");
                exit(2);
        }
	//Now, add args from passed-in parameters
	index = 9;
	char* tok = strtok(params, " ");
	while (tok!=NULL) {
		sprintf(rgargv[index], tok);
		index++;
		tok = strtok(NULL, " ");
	}

        /* Run rupture generator */
#ifdef _USE_MEMCACHED
        mc_genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP, mc_server);
#else
        genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP);
#endif

        //Round strike, dip, rake to integers
        //Round flen, fwid, dtop, shyp, dhyp, slip1, slip2, slip3 to 4 places
        for (i=0; i<srf->srf_prect.nseg; i++) {
                srf->srf_prect.prectseg[i].stk = floor(srf->srf_prect.prectseg[i].stk+0.5);
                srf->srf_prect.prectseg[i].dip = floor(srf->srf_prect.prectseg[i].dip+0.5);
                srf->srf_prect.prectseg[i].flen = floor(srf->srf_prect.prectseg[i].flen*10000+0.5)/10000.0;
                srf->srf_prect.prectseg[i].fwid = floor(srf->srf_prect.prectseg[i].fwid*10000+0.5)/10000.0;
                srf->srf_prect.prectseg[i].dtop = floor(srf->srf_prect.prectseg[i].dtop*10000+0.5)/10000.0;
                srf->srf_prect.prectseg[i].shyp = floor(srf->srf_prect.prectseg[i].shyp*10000+0.5)/10000.0;
                srf->srf_prect.prectseg[i].dhyp = floor(srf->srf_prect.prectseg[i].dhyp*10000+0.5)/10000.0;
        }
        for (i=0; i<srf->srf_apnts.np; i++) {
                srf->srf_apnts.apntvals[i].stk = floor(srf->srf_apnts.apntvals[i].stk+0.5);
                srf->srf_apnts.apntvals[i].dip = floor(srf->srf_apnts.apntvals[i].dip+0.5);
                srf->srf_apnts.apntvals[i].rake = floor(srf->srf_apnts.apntvals[i].rake+0.5);
                srf->srf_apnts.apntvals[i].slip1 = floor(srf->srf_apnts.apntvals[i].slip1*10000+0.5)/10000.0;
                srf->srf_apnts.apntvals[i].slip2 = floor(srf->srf_apnts.apntvals[i].slip2*10000+0.5)/10000.0;
                srf->srf_apnts.apntvals[i].slip3 = floor(srf->srf_apnts.apntvals[i].slip3*10000+0.5)/10000.0;
        }

        /* Free pseudo-program args */
        for (j = 0; j < MAX_ARGS; j++) {
            free(rgargv[j]);
        }

        free(rgargv);
        free(srf_out_file);
        return(0);
}

/* Version of rupgen_genslip that takes a seed.
 */

int rupgen_genslip_seed(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf, int hypo_dist, float dt, int seed) {
	char* seed_string = malloc(sizeof(char) * MAX_STR_ARGV);
	sprintf(seed_string, "seed=%d", seed);
	rupgen_genslip_with_params(rup_geom_file, slip, hypo, stats, srf, hypo_dist, dt, seed_string);
	free(seed_string);
}


int rupgen_genslip(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf, int hypo_dist, float dt) {
	int write_srf_param = 1;

        int i, j, rgargc;
        char **rgargv = NULL;

#ifdef _USE_MEMCACHED
        if (strcmp(mc_server, "")==0) {
                fprintf(stderr, "Error:  need to set memcached server with set_memcached_server()\n");
                exit(1);
        }
#endif

	char* srf_out_file = malloc(strlen(rup_geom_file)+20);
	sprintf(srf_out_file, "%s_s%d_h%d.srf", basename(rup_geom_file), slip, hypo);

    /* Create pseudo-program args */
    rgargc = NUM_GENSLIP_ARGS;
    rgargv = malloc(rgargc * sizeof(char*));
    for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
	    rgargv[j] = malloc(MAX_STR_ARGV);
		memset(rgargv[j], 0, MAX_STR_ARGV);
    }
    sprintf(rgargv[0], "%s", "rupgen");
    sprintf(rgargv[1], "infile=%s", rup_geom_file);
	sprintf(rgargv[2], "doslip=%d", slip);
	sprintf(rgargv[3], "dohypo=%d", hypo);
	sprintf(rgargv[4], "write_srf=%d", write_srf_param);
	sprintf(rgargv[5], "outfile=%s", srf_out_file);
	sprintf(rgargv[6], "dt=%f", dt);
	if (hypo_dist==RUPGEN_RANDOM_HYPO) {
		sprintf(rgargv[7], "uniformgrid_hypo=0");
		sprintf(rgargv[8], "random_hypo=1");
	} else if (hypo_dist==RUPGEN_UNIFORM_HYPO) {
                sprintf(rgargv[7], "uniformgrid_hypo=1");
                sprintf(rgargv[8], "random_hypo=0");
	} else {
		fprintf(stderr, "Error, did not specify a valid hypocenter location distribution, aborting.\n");
		exit(2);
	}

	for (i=0; i<9; i++) {
		printf("rupgen arg %d: %s\n", i, rgargv[i]);
	}

        /* Run rupture generator */
#ifdef _USE_MEMCACHED
        mc_genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP, mc_server);
#else
	genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP);
#endif

	//Round strike, dip, rake to integers
        //Round flen, fwid, dtop, shyp, dhyp, slip1, slip2, slip3 to 4 places
	for (i=0; i<srf->srf_prect.nseg; i++) {
		srf->srf_prect.prectseg[i].stk = floor(srf->srf_prect.prectseg[i].stk+0.5);
                srf->srf_prect.prectseg[i].dip = floor(srf->srf_prect.prectseg[i].dip+0.5);
		srf->srf_prect.prectseg[i].flen = floor(srf->srf_prect.prectseg[i].flen*10000+0.5)/10000.0;
		srf->srf_prect.prectseg[i].fwid = floor(srf->srf_prect.prectseg[i].fwid*10000+0.5)/10000.0;
		srf->srf_prect.prectseg[i].dtop = floor(srf->srf_prect.prectseg[i].dtop*10000+0.5)/10000.0;
		srf->srf_prect.prectseg[i].shyp = floor(srf->srf_prect.prectseg[i].shyp*10000+0.5)/10000.0;
                srf->srf_prect.prectseg[i].dhyp = floor(srf->srf_prect.prectseg[i].dhyp*10000+0.5)/10000.0;
	}
	for (i=0; i<srf->srf_apnts.np; i++) {
		srf->srf_apnts.apntvals[i].stk = floor(srf->srf_apnts.apntvals[i].stk+0.5);
                srf->srf_apnts.apntvals[i].dip = floor(srf->srf_apnts.apntvals[i].dip+0.5);
                srf->srf_apnts.apntvals[i].rake = floor(srf->srf_apnts.apntvals[i].rake+0.5);
		srf->srf_apnts.apntvals[i].slip1 = floor(srf->srf_apnts.apntvals[i].slip1*10000+0.5)/10000.0;
		srf->srf_apnts.apntvals[i].slip2 = floor(srf->srf_apnts.apntvals[i].slip2*10000+0.5)/10000.0;
		srf->srf_apnts.apntvals[i].slip3 = floor(srf->srf_apnts.apntvals[i].slip3*10000+0.5)/10000.0;
	}

        /* Free pseudo-program args */
        for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
            free(rgargv[j]);
        }

        free(rgargv);
	free(srf_out_file);
        return(0);
}

/* This function returns the initial seed used for a (source_id, rupture_id) tuple.
 * This seed can then be passed into genslip using the genslip_with_seed() function,
 * and genslip will automatically determine the correct seed to use with each rupture 
 * variation.
 * Otherwise, seed defaults to 0 for all ruptures.
 */
int rupgen_get_rupture_seed(int src_id, int rup_id) {
	return src_id*10000 + rup_id;
}

/* This function returns the seed used for a specific (source_id, rupture_id, slip, hypo) tuple.
 * This is useful for hb_high or other applications which don't call genslip but need to 
 * use the same seed.
 */

int rupgen_get_variation_seed(char* rup_geom_file, rg_stats_t *stats, int hypo_dist, int src_id, int rup_id, int slip, int hypo) {
    /* Create pseudo-program args */
    int MAX_ARGS = 30;
    int index = 0;
	int rgargc, j;
	char** rgargv;
    rgargc = MAX_ARGS;
    rgargv = malloc(rgargc * sizeof(char*));
    for (j = 0; j < MAX_ARGS; j++) {
        rgargv[j] = malloc(MAX_STR_ARGV);
        memset(rgargv[j], 0, MAX_STR_ARGV);
    }
    sprintf(rgargv[0], "%s", "rupgen");
    sprintf(rgargv[1], "infile=%s", rup_geom_file);
    if (hypo_dist==RUPGEN_RANDOM_HYPO) {
        sprintf(rgargv[2], "uniformgrid_hypo=0");
        sprintf(rgargv[3], "random_hypo=1");
    } else if (hypo_dist==RUPGEN_UNIFORM_HYPO) {
        sprintf(rgargv[2], "uniformgrid_hypo=1");
        sprintf(rgargv[3], "random_hypo=0");
    } else {
    	fprintf(stderr, "Error, did not specify a valid hypocenter location distribution, aborting.\n");
    	exit(2);
    }
    sprintf(rgargv[4], "doslip=%d", slip);
    sprintf(rgargv[5], "dohypo=%d", hypo);
	//Set initial seed
	int initial_seed = rupgen_get_rupture_seed(src_id, rup_id);
	sprintf(rgargv[6], "seed=%d", initial_seed);

	/* Run rupture generator */
    struct standrupformat srf;
	int seed;
#ifdef _USE_MEMCACHED
    seed = mc_genslip(rgargc, rgargv, stats, &srf, GET_SEED, mc_server);
#else
    seed = genslip(rgargc, rgargv, stats, &srf, GET_SEED);
#endif
    /* Free pseudo-program args */
    for (j = 0; j < MAX_ARGS; j++) {
  		free(rgargv[j]);
    }
    free(rgargv);

	return seed;	
}
