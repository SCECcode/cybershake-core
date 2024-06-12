#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include "string.h"

#include "rupgen_api.h"

void get_rupture_root(char* rup_root) {
	//Determine executable location - Linux-only
	char exe_loc[256];
	readlink("/proc/self/exe", exe_loc, 256);
	//Determine directory
	char* tok, *last_tok;
	char directory[256];
	directory[0] = '\0';
	last_tok = strtok(exe_loc, "/");
	tok = strtok(NULL, "/");
	while (tok!=NULL) {
		sprintf(directory, "%s/%s", directory, last_tok);
		last_tok = tok;
		tok = strtok(NULL, "/");
	}
	//Go back 3 levels
	sprintf(directory, "%s/../../../", directory);
	//Config file
	char cs_config_file[512];
	sprintf(cs_config_file, "%s/cybershake.cfg", directory);
	FILE* fp_in = fopen(cs_config_file, "r");
	char key[256];
	char value[256];
	while (fscanf(fp_in, "%s = %s", key, value)>0) {
		if (strcmp(key, "RUPTURE_ROOT")==0) {
			fclose(fp_in);
			sprintf(rup_root, value);
			return;
		}
	}
	printf("Couldn't find key RUPTURE_ROOT in configuration file %s.\n", cs_config_file);
	exit(1);
	return;
}


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

	//Determine rupture root from cybershake.cfg
	char rupture_root[256];
	get_rupture_root(rupture_root);

    set_memcached_server("localhost");
    char* rupture_file = malloc(sizeof(char)*256);
    sprintf(rupture_file, "%s/Ruptures_erf%d/%d/%d/%d_%d.txt", rupture_root, erf_id, src, rup, src, rup);
    char outfile[256];
    rg_stats_t stats;
    int slip, hypo;
    rupgen_get_num_rv(rupture_file, &stats, RUPGEN_UNIFORM_HYPO);
    printf("%d slips, %d hypos.\n", stats.numslip, stats.numhypo);
    int i, j;
    struct standrupformat srf;
	struct timespec ts_start, ts_end;
	for (i=rv_start; i<rv_end; i++) {
    	for (j=0; j<stats.numhypo; j++) {
        	rupgen_genslip(rupture_file, i, j, &stats, &srf, RUPGEN_UNIFORM_HYPO, 0.05);
			sprintf(outfile, "%d_%d_%d_uniform.srf", src, rup, i);
			_write_srf(&srf, outfile, 0);
            free_srf_ptrs(&srf);
        }
	}
	free(rupture_file);
}

