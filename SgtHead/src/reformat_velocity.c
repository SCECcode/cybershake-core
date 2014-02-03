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

	//on Blue Waters, 64 GB per node -- 10 GB free for chunking (since you need in and out and 3 pieces of data)
	int chunk_size = (int)(((long)10*1024*1024*1024)/(long)(nx*nz));
	if (chunk_size>ny) {
		chunk_size = ny;
	}

	printf("Using chunk size %d.\n", chunk_size);	
	p_data = malloc((long)nx*(long)chunk_size*(long)nz*sizeof(float));
	s_data = malloc((long)nx*(long)chunk_size*(long)nz*sizeof(float));	
	d_data = malloc((long)nx*(long)chunk_size*(long)nz*sizeof(float));

	if (p_data==0 || s_data==0 || d_data==0) {
		fprintf(stderr, "Error allocating input buffers.\n");
		exit(1);
	}

	awp_buffer = malloc(3*chunk_size*sizeof(float));

	int y_index = 0;
	int max_y;
	long in_arr_offset, outfile_offset;
	while (y_index<ny) {
		max_y = chunk_size;
		if (y_index+chunk_size>ny) {
			max_y = ny-y_index;
		}

		fread(p_data, sizeof(float), (long)nx*(long)max_y*(long)nz, p_fp);
		fread(s_data, sizeof(float), (long)nx*(long)max_y*(long)nz, s_fp);
		fread(d_data, sizeof(float), (long)nx*(long)max_y*(long)nz, d_fp);

		for(k=0; k<nz; k++) {
			//printf("%d-%d:  Z slice %d of %d.\n", y_index, y_index+max_y, k+1, nz);
			//fflush(stdout);
			for (i=0; i<nx; i++) {
				for (j=0; j<max_y; j++) {
					in_arr_offset = (long)j*nx*nz + (long)k*nx + (long)i;
					//printf("In offset: %d\n", in_arr_offset);
					awp_buffer[3*j] = p_data[in_arr_offset]*1000.0;
					awp_buffer[3*j+1] = s_data[in_arr_offset]*1000.0;
					awp_buffer[3*j+2] = d_data[in_arr_offset]*1000.0;
				}
				outfile_offset = 3*sizeof(float)*((long)k*nx*ny + (long)i*ny + (long)y_index);
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
