#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<8) {
		printf("Usage: %s <nx> <ny> <nz> <vp> <vs> <rho> <output file>\n", argv[0]);
		exit(1);
	}
	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);
	float vp = atof(argv[4]);
	float vs = atof(argv[5]);
	float rho = atof(argv[6]);
	char* output_filename = argv[7];

	FILE* fp_out = fopen(output_filename, "wb");
	//Write 3*nx*ny elements a time
	int i;
	int buf_num = 3*nx*ny;
	float* buffer = malloc(buf_num*sizeof(float));
	for (i=0; i<nx*ny; i++) {
		buffer[3*i] = vp;
		buffer[3*i+1] = vs;
		buffer[3*i+2] = rho;
	}
	for (i=0; i<nz; i++) {
		fwrite(buffer, sizeof(float), buf_num, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
	return 0;
}
