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

	int i, j, k, out_offset;
	float* p_data, *s_data, *d_data, *awp_buffer;


	//RWG is fast x, z, y (km)
	//AWP is fast y, x, z (actually fast x, y, z, but y and x are flipped)

	int chunk_size = (int)(175000000.0/(nx*nz));
	if (chunk_size>ny) {
		chunk_size = ny;
	}

	printf("Using chunk size %d.\n", chunk_size);	
	p_data = malloc(nx*chunk_size*nz*sizeof(float));
	s_data = malloc(nx*chunk_size*nz*sizeof(float));	
	d_data = malloc(nx*chunk_size*nz*sizeof(float));

	awp_buffer = malloc(3*chunk_size*sizeof(float));

	int y_index = 0;
	int in_arr_offset, max_y;
	long outfile_offset;
	while (y_index<ny) {
		max_y = chunk_size;
		if (y_index+chunk_size>ny) {
			max_y = ny-y_index;
		}

		fread(p_data, sizeof(float), nx*max_y*nz, p_fp);
		fread(s_data, sizeof(float), nx*max_y*nz, s_fp);
		fread(d_data, sizeof(float), nx*max_y*nz, d_fp);

		for(k=0; k<nz; k++) {
			//printf("%d-%d:  Z slice %d of %d.\n", y_index, y_index+max_y, k+1, nz);
			//fflush(stdout);
			for (i=0; i<nx; i++) {
				for (j=0; j<max_y; j++) {
					in_arr_offset = j*nx*nz + k*nx + i;
					//printf("In offset: %d\n", in_arr_offset);
					awp_buffer[3*j] = p_data[in_arr_offset]*1000.0;
					awp_buffer[3*j+1] = s_data[in_arr_offset]*1000.0;
					awp_buffer[3*j+2] = d_data[in_arr_offset]*1000.0;
				}
				outfile_offset = 3*sizeof(float)*(k*nx*ny + i*ny + y_index);
				fseek(out_fp, outfile_offset, SEEK_SET);
				fwrite(awp_buffer, sizeof(float), 3*max_y, out_fp);
			}
		}
		y_index += max_y;
	}

	fclose(p_fp);
	fclose(s_fp);
	fclose(d_fp);
	fflush(out_fp);
	fclose(out_fp);
	
	free(p_data);
	free(s_data);
	free(d_data);
	free(awp_buffer);
	return 0;
}
