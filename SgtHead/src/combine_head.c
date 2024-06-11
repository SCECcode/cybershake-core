#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "structure.h"

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <raw sgt in> <header in> <merged sgt out>\n", argv[0]);
		exit(1);
	}
	
	char* sgt_in = argv[1];
	char* header_in = argv[2];
	char* sgt_out = argv[3];

	FILE* sgt_fp_in = fopen(sgt_in, "rb");
	FILE* header_fp_in = fopen(header_in, "rb");
	FILE* sgt_fp_out = fopen(sgt_out, "wb");

	//Read in all the header data to parse SGTs appropriately
	struct sgtmaster mast;
	fread(&mast, sizeof(struct sgtmaster), 1, header_fp_in);
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex)*mast.globnp);
	fread(sgtindx, sizeof(struct sgtindex), mast.globnp, header_fp_in);
	struct sgtheader* sgtheads = malloc(sizeof(struct sgtheader)*mast.globnp);
	fread(sgtheads, sizeof(struct sgtheader), mast.globnp, header_fp_in);
	fclose(header_fp_in);

	fwrite(&mast, sizeof(struct sgtmaster), 1, sgt_fp_out);
	fwrite(sgtindx, sizeof(struct sgtindex), mast.globnp, sgt_fp_out);
	//Number of points to read at once
	int CHUNK_SIZE = 1000;
	int points_to_read;
	int point_index = 0;
	long long read_size;
	int i;
	float* buf = malloc(sizeof(float)*6*mast.nt*CHUNK_SIZE);
	while (point_index<mast.globnp) {
		printf("Point %d of %d.\n", point_index, mast.globnp);
		points_to_read = CHUNK_SIZE;
		if (point_index + points_to_read > mast.globnp) {
			points_to_read = mast.globnp - point_index;
		}
		read_size = 4*6*mast.nt*points_to_read;
		fread(buf, sizeof(float), 6*mast.nt*points_to_read, sgt_fp_in);
		for (i=0; i<points_to_read; i++) {
			fwrite(&(sgtheads[point_index+i]), sizeof(struct sgtheader), 1, sgt_fp_out);
			fwrite(&(buf[i*6*mast.nt]), sizeof(float), 6*mast.nt, sgt_fp_out);
		}
		point_index += points_to_read;
	}
	fclose(sgt_fp_in);
	fflush(sgt_fp_out);
	fclose(sgt_fp_out);

	free(sgtindx);
	free(sgtheads);	
	free(buf);
}

