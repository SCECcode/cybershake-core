#include <stdio.h>
#include <stdlib.h>
#include "structure.h"

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <sgt_in> <raw_sgt_out> <head_out>", argv[0]);
		exit(1);
	}

	FILE* fp_in = fopen(argv[1], "rb");
	FILE* sgt_out = fopen(argv[2], "wb");
	FILE* head_out = fopen(argv[3], "wb");
	struct sgtmaster mast;
	fread(&mast, sizeof(struct sgtmaster), 1, fp_in);
	fwrite(&mast, sizeof(struct sgtmaster), 1, head_out);
	int i;
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex)*mast.globnp);
	fread(sgtindx, sizeof(struct sgtindex), mast.globnp, fp_in);
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
