#include <stdio.h>
#include <stdlib.h>
#include "structure.h"

int main(int argc, char** argv) {
	if (argc<7) {
		printf("Usage: %s <6-component sgt data> <nt> <dt> <high-pass freq> <low-pass freq> <output file>\n", argv[0]);
		exit(1);
	}

	char* sgt_file = argv[1];
	int nt = atoi(argv[2]);
	float dt = atof(argv[3]);
	float fhi = atof(argv[4]);
	float flo = atof(argv[5]);
	char* sgt_output_file = argv[6];
	int i;
	const int SGT_COMPS = 6;

	FILE* fp_in = fopen(sgt_file, "rb");
	FILE* fp_out = fopen(sgt_output_file, "wb");
	float* sgt_buffer[6];
	for (i=0; i<SGT_COMPS; i++) {
		sgt_buffer[i] = malloc(sizeof(float)*nt);
	}

	//Use wcc_tfilter to do the filtering
	char** param_string = malloc(sizeof(char*) * 5);
	for (i=0; i<5; i++) {
		param_string[i] = malloc(sizeof(char) * 128);
	}
	sprintf(param_string[1], "order=4");
	sprintf(param_string[2], "fhi=%f", fhi);
	sprintf(param_string[3], "flo=%f", flo);
	struct statdata shead1;
	shead1.nt = nt;
	shead1.dt = dt;
	printf("%s\n", param_string[0]);

	for (i=0; i<SGT_COMPS; i++) {
		fread(sgt_buffer[i], sizeof(float), nt, fp_in);
		wcc_tfilter(4, param_string, sgt_buffer[i], &shead1);
		fwrite(sgt_buffer[i], sizeof(float), nt, fp_out);
	}

	fclose(fp_in);
	fflush(fp_out);
	fclose(fp_out);
	for (i=0; i<5; i++) {
		free(param_string[i]);
	}
	free(param_string);

	for (i=0; i<SGT_COMPS; i++) {
		free(sgt_buffer[i]);
	}
}

