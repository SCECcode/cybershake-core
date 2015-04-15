#include "../include.h"
#include "../structure.h"
#include "../functions.h"
#include "../defs.h"

#define NUM_INTERP 3
#define MAX_PERIODS 200
#define MAXPTS 2000000

void write_result(FILE* fp_out, struct seisheader* seis_head, float* rotD100, int* rotD100ang, float* rotD50, int num_periods, float periods[], int inter_flag);
void calc_acc(float* acc, float* seis_data, int num_comps, int nt, float dt);
void copy_result(struct rotD_record* output, float* rotD100, int* rotD100ang, float* rotD50, int num_periods, float periods[], int inter_flag);

extern void calc_rotd_(int* inter_flag, int* npairs, int* npts, float* dt, float* acc1, float* acc2, float* rotD100, int* rD100ang, float* rotD50, int* nFreq, float* RSP_Period);

int rotd(struct seisheader* header, float* seis_data, struct rotD_record* rotD_records) {
	int i;
	//Constants for rotDcalc
	int inter_flag = 1;
	int npairs = 1;

	//Allocate buffers
	int num_comps;
	if (header->comps == (X_COMP_FLAG | Y_COMP_FLAG | Z_COMP_FLAG)) {
		num_comps = 3;
	} else if (header->comps == (X_COMP_FLAG | Y_COMP_FLAG )) {
		num_comps = 2;
	} else if (header->comps == X_COMP_FLAG || header->comps == Y_COMP_FLAG || header->comps == Z_COMP_FLAG) {
		num_comps = 1;
	} else {
		fprintf(stderr, "Could not determine which components are included in this seismogram, aborting rotd.\n");
		num_comps = 0;
		return -1;
	}

	if (num_comps!=2) {
		fprintf(stderr, "This RotD code must be used with two-component seismograms.\n");
		return 2;
	}

	float* acc = check_malloc(sizeof(float) * header->nt * num_comps);
	float* rotD100 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	float* rotD50 = check_malloc(sizeof(float*) * NUM_INTERP * MAX_PERIODS);
	int* rD100ang = check_malloc(sizeof(float*) * NUM_INTERP * MAX_PERIODS);

	//Added periods for 1 Hz	
	int num_periods = 22;
	float periods[] = {1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10.0};

	calc_acc(acc, seis_data, num_comps, header->nt, header->dt);
	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), acc, acc+header->nt, rotD100, rD100ang, rotD50, &num_periods, periods);
//	write_result(fp_out, header, rotD100, rD100ang, rotD50, num_periods, periods, inter_flag);
	copy_result(rotD_records, rotD100, rD100ang, rotD50, num_periods, periods, inter_flag);

//	fflush(stderr);
//	fflush(stdout);
	free(rotD100);
	free(rotD50);
	free(rD100ang);
	free(acc);

	return 0;
}

void calc_acc(float* acc, float* seis_data, int num_comps, int nt, float dt) {
	int i, j;
	const float g2cm = 980.665;
	for (i=0; i<num_comps; i++) {
		int comp_offset = i*nt;
		acc[comp_offset] = 0.0;
		for (j=1; j<nt; j++) {
			acc[comp_offset+j] = (seis_data[comp_offset+j]-seis_data[comp_offset+j-1])/(g2cm * dt);
		}
	}
}

void copy_result(struct rotD_record* output, float* rotD100, int* rotD100ang, float* rotD50, int num_periods, float periods[], int inter_flag) {
	int i;
	for (i=0; i<num_periods; i++) {
		output[i].period = periods[i];
		output[i].rotD100 = rotD100[i*NUM_INTERP + (inter_flag-1)];
		output[i].rotD100_angle = rotD100ang[i*NUM_INTERP + (inter_flag-1)];
		output[i].rotD50 = rotD50[i*NUM_INTERP + (inter_flag-1)];
	}
}

void write_result(FILE* fp_out, struct seisheader* seis_head, float* rotD100, int* rotD100ang, float* rotD50, int num_periods, float periods[], int inter_flag) {
	fwrite(seis_head, sizeof(struct seisheader), 1, fp_out);
	//Write number of periods next
	fwrite(&num_periods, sizeof(int), 1, fp_out);
	int i;
	struct rotD_record rotDrec;
	for (i=0; i<num_periods; i++) {
		rotDrec.period = periods[i];
		rotDrec.rotD100 = rotD100[i*NUM_INTERP + (inter_flag-1)];
		rotDrec.rotD100_angle = rotD100ang[i*NUM_INTERP + (inter_flag-1)];
		rotDrec.rotD50 = rotD50[i*NUM_INTERP + (inter_flag-1)];
		fwrite(&rotDrec, sizeof(struct rotD_record), 1, fp_out);
	}
}
