#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <input mesh> <np> <output mesh prefix>", argv[0]);
		exit(1);
	}

	char* input_mesh = argv[1];
	long long np = atol(argv[2]);
	char* output_prefix = argv[3];

	long long i;
	int j;
	FILE* fp_in = fopen(input_mesh, "rb");
	char output_files[3][256];
	sprintf(output_files[0], "%s.vp", output_prefix);
	sprintf(output_files[1], "%s.vs", output_prefix);
	sprintf(output_files[2], "%s.rho", output_prefix);
	FILE* fp_outs[3];
	for (j=0; j<3; j++) {
		fp_outs[j] = fopen(output_files[j], "wb");
	}

	long long chunk_pts = 1000000000;
	float* buffer = malloc(sizeof(float)*3*chunk_pts);
	float* comp_buffers[3];
	for (i=0; i<3; i++) {
		comp_buffers[i] = malloc(sizeof(float)*chunk_pts);
	}
	long long index = 0;
	while (index+chunk_pts<np) {
		printf("Reading points %lld to %lld.\n", index, index+chunk_pts);
		fread(buffer, chunk_pts, 3*sizeof(float), fp_in);
		for (i=0; i<chunk_pts; i++) {
			for (j=0; j<3; j++) {
				comp_buffers[j][i] = buffer[3*i+j];
			}
		}
		for (j=0; j<3; j++) {
			fwrite(comp_buffers[j], chunk_pts, sizeof(float), fp_outs[j]);
		}
		index += chunk_pts;
	}
	if (index<np) {
		long long num_pts = np-index;
        printf("Reading points %lld to %lld.\n", index, index+num_pts);
		fread(buffer, num_pts, 3*sizeof(float), fp_in);
		for (i=0; i<num_pts; i++) {
            for (j=0; j<3; j++) {
                comp_buffers[j][i] = buffer[3*i+j];
            }
        }
        for (j=0; j<3; j++) {
            fwrite(comp_buffers[j], num_pts, sizeof(float), fp_outs[j]);
        }
	}	

	free(buffer);
	for (j=0; j<3; j++) {
		fflush(fp_outs[j]);
		fclose(fp_outs[j]);
		free(comp_buffers[j]);
	}
	fclose(fp_in);
}
