#include <stdio.h>
#include <stdlib.h>

#include "structure.h"

//This code inputs an sgt header file and outputs a list of the points SGTs are saved for.

int main(int argc, char** argv) {
	if (argc<3) {
		printf("Usage: %s <sgtheadfile> <output file>\n", argv[0]);
		exit(1);
	}

	char* sgtheadfile = argv[1];
	char* outputfile = argv[2];

	FILE* fp_in = fopen(sgtheadfile, "rb");
	FILE* fp_out = fopen(outputfile, "wb");

	struct sgtmaster sgtmast;
	fread(&sgtmast, sizeof(struct sgtmaster), 1, fp_in);
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex)*sgtmast.globnp);
	fread(sgtindx, sizeof(struct sgtindex), sgtmast.globnp, fp_in);
	struct sgtheader* sgthead = malloc(sizeof(struct sgtheader)*sgtmast.globnp);
	fread(sgthead, sizeof(struct sgtheader), sgtmast.globnp, fp_in);
	fclose(fp_in);

	int i;
	for (i=0; i<sgtmast.globnp; i++) {
		fprintf(fp_out, "%d %d %d\n", sgthead[i].xsgt, sgthead[i].ysgt, sgthead[i].zsgt);
	}
	fflush(fp_out);
	fclose(fp_out);
}
