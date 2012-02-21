#include "include.h"
#include "structure.h"
#include "defs.h"
#include "function.h"

int rupgen_get_num_rv(char* rup_geom_file, rg_stats_t *stats) {
	int writeout = 0;

	int j, rgargc;
	char **rgargv = NULL;

	/* Create pseudo-program args */
	rgargc = NUM_GENSLIP_ARGS;
	rgargv = malloc(rgargc * sizeof(char*));
	for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
		rgargv[j] = malloc(MAX_STR_ARGV);
	}
	sprintf(rgargv[0], "%s", "rupgen");
	sprintf(rgargv[1], "infile=%s", rup_geom_file);
	sprintf(rgargv[2], "writeout=%d", writeout);

	/* Run rupture generator */
	struct standrupformat srf;
	genslip(rgargc, rgargv, stats, &srf, GET_STATS);
  
	/* Free pseudo-program args */
	for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
	    free(rgargv[j]);
	}

  	free(rgargv);
	return(0);

}

int rupgen_genslip(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf) {
        int writeout = 0;

        int i, j, rgargc;
        char **rgargv = NULL;

        /* Create pseudo-program args */
        rgargc = NUM_GENSLIP_ARGS;
        rgargv = malloc(rgargc * sizeof(char*));
        for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
                rgargv[j] = malloc(MAX_STR_ARGV);
        }
        sprintf(rgargv[0], "%s", "rupgen");
        sprintf(rgargv[1], "infile=%s", rup_geom_file);
        sprintf(rgargv[2], "writeout=%d", writeout);
	sprintf(rgargv[3], "doslip=%d", slip);
	sprintf(rgargv[4], "dohypo=%d", hypo);

        /* Run rupture generator */
        genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP);

	//Round strike, dip, rake to integers
        //Round flen, fwid, dtop, shyp, dhyp, slip1, slip2, slip3 to 2 places
	for (i=0; i<srf->srf_prect.nseg; i++) {
		srf->srf_prect.prectseg[i].stk = floor(srf->srf_prect.prectseg[i].stk+0.5);
                srf->srf_prect.prectseg[i].dip = floor(srf->srf_prect.prectseg[i].dip+0.5);
		srf->srf_prect.prectseg[i].flen = floor(srf->srf_prect.prectseg[i].flen*100+0.5)/100.0;
		srf->srf_prect.prectseg[i].fwid = floor(srf->srf_prect.prectseg[i].fwid*100+0.5)/100.0;
		srf->srf_prect.prectseg[i].dtop = floor(srf->srf_prect.prectseg[i].dtop*100+0.5)/100.0;
		srf->srf_prect.prectseg[i].shyp = floor(srf->srf_prect.prectseg[i].shyp*100+0.5)/100.0;
                srf->srf_prect.prectseg[i].dhyp = floor(srf->srf_prect.prectseg[i].dhyp*100+0.5)/100.0;
	}
	for (i=0; i<srf->srf_apnts.np; i++) {
		srf->srf_apnts.apntvals[i].stk = floor(srf->srf_apnts.apntvals[i].stk+0.5);
                srf->srf_apnts.apntvals[i].dip = floor(srf->srf_apnts.apntvals[i].dip+0.5);
                srf->srf_apnts.apntvals[i].rake = floor(srf->srf_apnts.apntvals[i].rake+0.5);
		srf->srf_apnts.apntvals[i].slip1 = floor(srf->srf_apnts.apntvals[i].slip1*100+0.5)/100.0;
		srf->srf_apnts.apntvals[i].slip2 = floor(srf->srf_apnts.apntvals[i].slip2*100+0.5)/100.0;
		srf->srf_apnts.apntvals[i].slip3 = floor(srf->srf_apnts.apntvals[i].slip3*100+0.5)/100.0;
	}

        /* Free pseudo-program args */
        for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
            free(rgargv[j]);
        }

        free(rgargv);
        return(0);
}
