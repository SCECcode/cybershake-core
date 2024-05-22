#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "structure.h"

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <sgt_in> <raw_sgt_out> <head_out> [-r]\n", argv[0]);
		printf("-r: recalculate indx, converting from short form to 2016 longer form\n");
		exit(1);
	}

	FILE* fp_in = fopen(argv[1], "rb");
	FILE* sgt_out = fopen(argv[2], "wb");
	FILE* head_out = fopen(argv[3], "wb");
	int recompute_indx = 0;
	if (argc==5 && strcmp(argv[4], "-r")==0) {
		recompute_indx = 1;
		printf("Converting indx.\n");
	}

	struct sgtmaster mast;
	fread(&mast, sizeof(struct sgtmaster), 1, fp_in);
	fwrite(&mast, sizeof(struct sgtmaster), 1, head_out);
	int i, xsgt, ysgt, zsgt;
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex)*mast.globnp);
	fread(sgtindx, sizeof(struct sgtindex), mast.globnp, fp_in);
	if (recompute_indx==1) {
		//Old way: xsgt*100000000 + ysgt*10000 + zsgt
		//New way: xsgt*1000000000000 + ysgt*1000000 + zsgt
		for (i=0; i<mast.globnp; i++) {
			xsgt = sgtindx[i].indx/100000000;
			ysgt = (sgtindx[i].indx/10000) % 10000;
			zsgt = sgtindx[i].indx % 10000;
			sgtindx[i].indx = xsgt*(long long)1000000000000 + ysgt*(long long)1000000 + zsgt;
		}
	}
	fwrite(sgtindx, sizeof(struct sgtindex), mast.globnp, head_out);
	struct sgtheader sgthead;
	float* sgt_data = malloc(sizeof(float)*6*mast.nt);
	for (i=0; i<mast.globnp; i++) {
		fread(&sgthead, sizeof(struct sgtheader), 1, fp_in);
		fread(sgt_data, sizeof(float), 6*mast.nt, fp_in);
		fwrite(&sgthead, sizeof(struct sgtheader), 1, head_out);
		fwrite(sgt_data, sizeof(float), 6*mast.nt, sgt_out);
	}
	fclose(fp_in);
	fflush(sgt_out);
	fclose(sgt_out);
	fflush(head_out);
	fclose(head_out);
	return 0;
}
