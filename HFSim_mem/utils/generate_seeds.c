#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "rupgen_api.h"

int main(int argc, char** argv) {
	if (argc<3) {
		printf("Usage: %s <input rvs> <output file> [-r]\n", argv[0]);
		exit(1);
	}
	char* input_file = argv[1];
	char* output_file = argv[2];

	int restart = 0;
    if (argc>3) {
    	if (strcmp(argv[3], "-r")==0) {
        	restart = 1;
    	}
    }

    int src_restart, rup_restart, rv_restart;

	char* RUPTURE_ROOT = "/gpfs/alpine/proj-shared/geo112/CyberShake/ruptures/Ruptures_erf36";

	rg_stats_t stats;

	FILE* fp_in, *fp_out;
	int src_id, rup_id, rv_id;
	char rup_geom_file[256];
	int rup_var_seed;

        if (restart) {
                fp_in = fopen(output_file, "r");
                while (fscanf(fp_in, "%d %d %d %f", &src_id, &rup_id, &rv_id, &rup_var_seed)!=EOF);
                printf("Last entry: %d %d %d\n", src_id, rup_id, rv_id);
                src_restart = src_id;
                rup_restart = rup_id;
                rv_restart = rv_id;
                fclose(fp_in);
        }
    fp_in = fopen(input_file, "r");
    fp_out = fopen(output_file, "a");
	int i=0;
	while (fscanf(fp_in, "%d %d %d", &src_id, &rup_id, &rv_id)!=EOF) {
		if (restart) {
        	if (src_id<src_restart) {
            	continue;
            } else if (src_id==src_restart) {
            	if (rup_id<rup_restart) {
                	continue;
                } else if (rup_id==rup_restart) {
                	if (rv_id<=rv_restart) {
                		continue;
                    }
                }
            }
        }
		sprintf(rup_geom_file, "%s/%d/%d/%d_%d.txt", RUPTURE_ROOT, src_id, rup_id, src_id, rup_id);
		rup_var_seed = rupgen_get_variation_seed(rup_geom_file, &stats, RUPGEN_UNIFORM_HYPO, src_id, rup_id, rv_id, 0);
		fprintf(fp_out, "%d %d %d %d\n", src_id, rup_id, rv_id, rup_var_seed);
		i++;
	}
	fclose(fp_in);
	fflush(fp_out);
	fclose(fp_out);
	return 0;
}

