#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "function.h"

extern void spectrad_(int *nx, int *ny, int *npts, float *dt, char* seis_units, int seis_units_len, char* output_units, int output_units_len, char* output_type, int output_type_len, char* period, int period_len, float *filter_high_hz, char* byteswap, int byteswap_len, char* input_file, int input_file_len, char* output_file, int output_file_len, float* seis);

int main(int argc, char** argv) {
	//command-line arguments
	char stat[64];
	float slon, slat;
	char sgt_xfilename[256];
	char sgt_yfilename[256];
	char rupmodfile[256];
	int ntout;
	int outputBinary = 1;
	int mergeOutput = 1;
	char seis_filename[256];	
	int i;

	memset(seis_filename, ' ', 256);

	setpar(argc, argv);

	mstpar("stat","s",stat);
	mstpar("slat","f",&slat);
	mstpar("slon","f",&slon);
	printf("slat %f, slon %f\n", slat, slon);
	mstpar("rupmodfile","s",rupmodfile);
	
	mstpar("sgt_xfile","s",sgt_xfilename);
	mstpar("sgt_yfile","s",sgt_yfilename);
	
	mstpar("ntout","d",&ntout);
	mstpar("seis_file","s",seis_filename);

	printf("seis_file: %s\n", seis_filename);

	//synthesis
	float* seis = jbsim3d_synth(stat, slon, slat, rupmodfile, sgt_xfilename, sgt_yfilename, outputBinary, mergeOutput, ntout, seis_filename);

	//psa
	int nx, ny, npts;
	float dt;
	char seis_units[120];
	char output_units[120];
	char output_type[120];
	char period[120];
	float filter_high_hz;
	char byteswap[120];
	char input_file[256];
	char output_file[256];
	int seis_units_len, output_units_len, output_type_len, period_len, byteswap_len, input_file_len, output_file_len;

	/*for (i=0; i<119; i++) {
		strcat(seis_units, " ");
		strcat(output_units, " ");
		strcat(output_type, " ");
		strcat(period, " ");
		strcat(byteswap, " ");
	}
	for (i=0; i<255; i++) {
		strcat(input_file, " ");
		strcat(output_file, " ");
	}*/

	memset(seis_units, ' ', 120);
	memset(output_units, ' ', 120);
	memset(output_type, ' ', 120);
	memset(period, ' ', 120);
	memset(byteswap, ' ', 120);
	memset(input_file, ' ', 256);
	memset(output_file, ' ', 256);
	
	mstpar("simulation_out_pointsX","d",&nx);
	mstpar("simulation_out_pointsY","d",&ny);
	mstpar("simulation_out_timesamples","d",&npts);
	mstpar("simulation_out_timeskip","f",&dt);
	mstpar("surfseis_rspectra_seismogram_units","s",seis_units);
	mstpar("surfseis_rspectra_output_units","s",output_units);
	mstpar("surfseis_rspectra_output_type","s",output_type);
	mstpar("surfseis_rspectra_period","s",period);
	mstpar("surfseis_rspectra_apply_filter_highHZ","f",&filter_high_hz);
	mstpar("surfseis_rspectra_apply_byteswap","s",byteswap);

	strcpy(input_file, seis_filename);
	mstpar("out","s",&output_file);

	//replace null terminator with space
	//fortran expects all the strings to be space-terminated
	seis_units_len = strlen(seis_units);
	seis_units[seis_units_len] = ' ';
	output_units_len = strlen(output_units);
	output_units[output_units_len] = ' ';
	output_type_len = strlen(output_type);
	output_type[output_type_len] = ' ';
	period_len = strlen(period);
	period[period_len] = ' ';
	byteswap_len = strlen(byteswap);
	byteswap[byteswap_len] = ' ';
	input_file_len = strlen(input_file);
	input_file[input_file_len] = ' ';
	output_file_len = strlen(output_file);
	output_file[output_file_len] = ' ';

	for (i=0; i<120; i++) {
		printf("%c", seis_units[i]);
	}
	printf("Entering spectrad.\n");
        //printf("locations: %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld, %ld\n", &nx, &ny, &npts, &dt, seis_units, output_units, output_type, period, &filter_high_hz, byteswap, input_file, output_file);
	fflush(stdout);
	spectrad_(&nx, &ny, &npts, &dt, seis_units, seis_units_len, output_units, output_units_len, output_type, output_type_len, period, period_len, &filter_high_hz, byteswap, byteswap_len, input_file, input_file_len, output_file, output_file_len, seis);

	free(seis);
}
