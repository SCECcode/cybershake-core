#include <stdio.h>
#include <stdlib.h>
#include "string.h"

//Alters awp:  swaps XX and YY, XZ and YZ, multiplies XZ and YZ by -1
//Insert in order (RWG): XX YY ZZ XY XZ YZ

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <input AWP SGT> <timesteps> <output AWP SGT> <-z>", argv[0]);
		exit(1);
	}

	int z_comp = 0;
	if (argc==5 && strcmp(argv[4], "-z")==0) {
		z_comp = 1;
		printf("Will double all SGTs since they're the Z source\n");
	}
	FILE* fp_in;
	FILE* fp_out;
	float* sgt_data;
	size_t bytes_read;
	size_t tot_bytes_read = 0;
	int timesteps = atoi(argv[2]);
	int i;
	sgt_data = malloc(sizeof(float)*6*timesteps);

	int XX = 0;
	int YY = 1;
	int ZZ = 2;
	int XY = 3;
	int XZ = 4;
	int YZ = 5;
	fp_in = fopen(argv[1], "rb");
	fp_out = fopen(argv[3], "wb");

	bytes_read = fread(sgt_data, sizeof(float), 6*timesteps, fp_in);
	while (bytes_read==6*timesteps) {
		tot_bytes_read += bytes_read;
		if (z_comp==1) {
			//double everything and negate for flipped source
			for (i=0; i<timesteps*6; i++) {
				sgt_data[i] = -2.0*sgt_data[i];
			}
		}
		//YY -> XX
		fwrite(sgt_data+YY*timesteps, sizeof(float), timesteps, fp_out);
		//XX -> YY
		fwrite(sgt_data+XX*timesteps, sizeof(float), timesteps, fp_out);
		//ZZ -> ZZ
		fwrite(sgt_data+ZZ*timesteps, sizeof(float), timesteps, fp_out);
		//XY -> XY
		fwrite(sgt_data+XY*timesteps, sizeof(float), timesteps, fp_out);
		//apply the -1 to XZ and YZ
                for (i=0; i<timesteps; i++) {
                        sgt_data[XZ*timesteps+i] = -1.0*sgt_data[XZ*timesteps+i];
                        sgt_data[YZ*timesteps+i] = -1.0*sgt_data[YZ*timesteps+i];
                }
		//-YZ -> XZ
		fwrite(sgt_data+YZ*timesteps, sizeof(float), timesteps, fp_out);
		//-XZ -> YZ
                fwrite(sgt_data+XZ*timesteps, sizeof(float), timesteps, fp_out);
	        bytes_read = fread(sgt_data, sizeof(float), 6*timesteps, fp_in);
	}
	fclose(fp_in);
	fflush(fp_out);
	fclose(fp_out);
	return 0;
}
