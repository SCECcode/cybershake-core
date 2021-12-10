#include <stdio.h>
#include <stdlib.h>

#include "structure.h"

/*Read in an RWG-formatted SGT, output a version with fewer timesteps by a provided factor. */

int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <SGT header in> <SGT in> <decimation factor> <SGT header out> <SGT out>\n", argv[0]);
		exit(1);
	}	

	char* header_in = argv[1];
	char* sgt_in = argv[2];
	int decimate = atoi(argv[3]);
	char* header_out = argv[4];
	char* sgt_out = argv[5];

	FILE* fp_in = fopen(header_in, "rb");
	FILE* fp_out = fopen(header_out, "wb");
	struct sgtmaster mast_in, mast_out;
	fread(&mast_in, sizeof(struct sgtmaster), 1, fp_in);
	int new_nt = mast_in.nt/decimate;
	if (decimate*new_nt != mast_in.nt) {
		printf("Original nt of %d isn't evenly divisible by decimation factor of %d, aborting.\n", mast_in.nt, decimate);
		exit(2);
	}

	mast_out.geoproj = mast_in.geoproj;
	mast_out.modellon = mast_in.modellon;
	mast_out.modellat = mast_in.modellat;
	mast_out.modelrot = mast_in.modelrot;
	mast_out.xshift = mast_in.xshift;
	mast_out.yshift = mast_in.yshift;
	mast_out.globnp = mast_in.globnp;
	mast_out.localnp = mast_in.localnp;
	mast_out.nt = new_nt;
	fwrite(&mast_out, sizeof(struct sgtmaster), 1, fp_out);
	
	struct sgtindex* sgtindx = malloc(sizeof(struct sgtindex) * mast_in.globnp);
	fread(sgtindx, sizeof(struct sgtindex), mast_in.globnp, fp_in);
	fwrite(sgtindx, sizeof(struct sgtindex), mast_out.globnp, fp_out);

	struct sgtheader* sgthead = malloc(sizeof(struct sgtheader) * mast_in.globnp);
	fread(sgthead, sizeof(struct sgtheader), mast_in.globnp, fp_in);
	int i;
	for (i=0; i<mast_in.globnp; i++) {
		sgthead[i].nt = new_nt;
		sgthead[i].dt *= decimate;
	}
	fwrite(sgthead, sizeof(struct sgtheader), mast_in.globnp, fp_out);
	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);

	fp_in = fopen(sgt_in, "rb");
	fp_out = fopen(sgt_out, "wb");

	int BUFFER_POINTS = 2000;
	float* buffer = malloc(sizeof(float)*6*mast_in.nt*BUFFER_POINTS);
	float* buffer_out = malloc(sizeof(float)*6*mast_out.nt*BUFFER_POINTS);
	int pt_index = 0;
	int num_pts_to_read = 0;
	int j, k, out_offset, in_offset;
	while (pt_index<mast_in.globnp) {
		if (pt_index+BUFFER_POINTS<mast_in.globnp) {
			num_pts_to_read = BUFFER_POINTS;
		} else {
			num_pts_to_read = mast_in.globnp - pt_index;
		}
		fread(buffer, sizeof(float), 6*mast_in.nt*num_pts_to_read, fp_in);
		for (i=0; i<num_pts_to_read; i++) {
			for (j=0; j<6; j++) {
				in_offset = i*6*mast_in.nt + j*mast_in.nt;
				out_offset = i*6*mast_out.nt + j*mast_out.nt;
				for (k=0; k<mast_in.nt; k++) {
					//Decimate timeseries
					if (k % decimate == 0) {
						buffer_out[out_offset+k/decimate] = buffer[in_offset+k];
					}
				}
			}
		}
		fwrite(buffer_out, sizeof(float), num_pts_to_read*6*mast_out.nt, fp_out);
		pt_index += num_pts_to_read;
	}

	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);
	
	free(sgtindx);
	free(sgthead);
	free(buffer);
	free(buffer_out);
}
