#include "include.h"
#include "structure.h"
#include "defs.h"
#include "function.h"

char rupture_dir[MAX_STR_ARGV] = {0};

int rupgen_init(char* rup_dir) {
	strcpy(rupture_dir, rup_dir);
	return 0;
}

int rupgen_get_num_rv(int source, int rupture, rg_stats_t *stats) {
	char infile[1024];
	char logfile[1024];
	if (strlen(rupture_dir)==0 ) {
		fprintf(stderr, "Rupture directory was not initalized with rupgen_init()\n");
		exit(1);
	}

	sprintf(infile, "%s/%d/%d/%d_%d.txt", rupture_dir, source, rupture, source, rupture);

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
	sprintf(rgargv[1], "infile=%s", infile);
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

int rupgen_genslip(int source, int rupture, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf) {
        char infile[1024];
        char logfile[1024];

        if (strlen(rupture_dir)==0 ) {
                fprintf(stderr, "Rupture directory was not initalized with rupgen_init()\n");
                exit(1);
        }

        sprintf(infile, "%s/%d/%d/%d_%d.txt", rupture_dir, source, rupture, source, rupture);

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
        sprintf(rgargv[1], "infile=%s", infile);
        sprintf(rgargv[2], "writeout=%d", writeout);
	sprintf(rgargv[3], "doslip=%d", slip);
	sprintf(rgargv[4], "dohypo=%d", hypo);

        /* Run rupture generator */
        genslip(rgargc, rgargv, stats, srf, RUN_GENSLIP);

        /* Free pseudo-program args */
        for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
            free(rgargv[j]);
        }

        free(rgargv);
        return(0);
}
