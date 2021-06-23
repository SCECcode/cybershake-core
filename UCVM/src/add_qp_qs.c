#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <input velocity mesh> <np> <qp> <qs> <output velocity mesh>\n", argv[0]);
		exit(1);
	}

	char* input_mesh = argv[1];
	int np = atoi(argv[2]);
	float qp = atof(argv[3]);
	float qs = atof(argv[4]);
	char* output_mesh = argv[5];

	FILE* fp_in = fopen(input_mesh, "rb");
	FILE* fp_out = fopen(output_mesh, "wb");

	int MAX_PTS_PER_READ = 150000000;
	float* buffer_in = malloc(sizeof(float) * MAX_PTS_PER_READ * 3);
	float* buffer_out = malloc(sizeof(float) * MAX_PTS_PER_READ * 5);
	int index = 0;
	int pts_to_read = 0;
	int i;
	while (index < np) {
		pts_to_read = MAX_PTS_PER_READ;
		if (index + pts_to_read > np) {
			pts_to_read = np - index;
		}
		fread(buffer_in, sizeof(float), pts_to_read * 3, fp_in);
		for (i=0; i<pts_to_read; i++) {
			buffer_out[5*i] = buffer_in[3*i];
			buffer_out[5*i+1] = buffer_in[3*i+1];
			buffer_out[5*i+2] = buffer_in[3*i+2];
			buffer_out[5*i+3] = qp;
			buffer_out[5*i+4] = qs;
		}
		fwrite(buffer_out, sizeof(float), pts_to_read * 5, fp_out);
		index += pts_to_read;
	}
	fclose(fp_in);
	fflush(fp_out);
	fclose(fp_out);
}
