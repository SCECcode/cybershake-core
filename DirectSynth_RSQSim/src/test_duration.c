#include "include.h"
#include "defs.h"
#include "structure.h"
#include "duration.h"
#include "functions.h"

int my_global_id = 0;

int main(int argc, char** argv) {
	if (argc<4) {
		printf("Usage: %s <num rvs> <input grm> <output duration file>", argv[0]);
		exit(1);
	}
	int num_rvs = atoi(argv[1]);
	int num_comps = 2;
	char* input_file = argv[2];
	char* output_file = argv[3];
	int i;

	FILE* fp_in = fopen(input_file, "rb");
	FILE* fp_out = fopen(output_file, "wb");
	struct duration_record* duration_data = check_malloc(sizeof(int)+sizeof(struct duration_record)*num_comps*NUM_DURATION_MEASURES);

	struct seisheader* headers = check_malloc(sizeof(struct seisheader)*num_rvs);
	float** seis = check_malloc(sizeof(float*) * num_rvs);
	for (i=0; i<num_rvs; i++) {
		//Read seismogram
		fread(&headers[i], sizeof(struct seisheader), 1, fp_in);
		seis[i] = check_malloc(sizeof(float)*headers[i].nt*2);
		fread(seis[i], sizeof(float), headers[i].nt*2, fp_in);
		int num_dur_measures = NUM_DURATION_MEASURES;
		memcpy(duration_data, &num_dur_measures, sizeof(int));
		char* start_ptr = ((char*)(duration_data))+sizeof(int);
		duration(&headers[i], seis[i], (struct duration_record*)start_ptr);
		fwrite(&headers[i], sizeof(struct seisheader), 1, fp_out);
		fwrite(duration_data, sizeof(int) + sizeof(struct duration_record)*num_comps*NUM_DURATION_MEASURES, 1, fp_out);
	}
	fflush(fp_out);
	fclose(fp_out);
	fclose(fp_in);
	for (i=0; i<num_rvs; i++) {
		free(seis[i]);
	}
	free(seis);
	free(duration_data);
}
