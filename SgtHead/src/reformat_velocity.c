#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if (argc<6) {
		printf("Usage: %s <nx> <ny> <nz> <velocity file prefix> <awp velocity out>", argv[0]);
		exit(0);
	}
	
	int nx = atoi(argv[1]);
	int ny = atoi(argv[2]);
	int nz = atoi(argv[3]);

	char p_name[256];
	char s_name[256];
	char d_name[256];

	sprintf(p_name, "%s.p", argv[4]);
	sprintf(s_name, "%s.s", argv[4]);
	sprintf(d_name, "%s.d", argv[4]);

	FILE* p_fp = fopen(p_name, "rb");
	FILE* s_fp = fopen(s_name, "rb");
	FILE* d_fp = fopen(d_name, "rb");
	FILE* out_fp = fopen(argv[5], "wb");

	int i, j, k, in_offset, out_offset;
	float* p_data, *s_data, *d_data, *out_data;

	p_data = malloc(nx*ny*nz*sizeof(float));
        s_data = malloc(nx*ny*nz*sizeof(float));
        d_data = malloc(nx*ny*nz*sizeof(float));
	out_data = malloc(3*nx*ny*nz*sizeof(float));

	fread(p_data, sizeof(float), nx*ny*nz, p_fp);
        fread(s_data, sizeof(float), nx*ny*nz, s_fp);
        fread(d_data, sizeof(float), nx*ny*nz, d_fp);

	fclose(p_fp);
	fclose(s_fp);
	fclose(d_fp);

	for(i=0; i<nz; i++) {
		printf("%d of %d z-slices.\n", i, nz);
		for (j=0; j<nx; j++) {
			for (k=0; k<ny; k++) {
				in_offset = k*nx*nz + i*nx + j;
				out_offset = 3*(i*nx*ny + j*ny + k);
				out_data[out_offset] = p_data[in_offset]*1000.0;
				out_data[out_offset+1] = s_data[in_offset]*1000.0;
				out_data[out_offset+2] = d_data[in_offset]*1000.0;
			}
		}
	}

	fwrite(out_data, sizeof(float), 3*nx*ny*nz, out_fp);
	fflush(out_fp);
	fclose(out_fp);
	return 0;
}
