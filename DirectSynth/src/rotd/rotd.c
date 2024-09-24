#include "../include.h"
#include "../defs.h"
#include "../structure.h"
#include "../functions.h"

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

	if (num_comps<2) {
		fprintf(stderr, "This RotD code must be used with two- or three-component seismograms.\n");
		return 2;
	}

	float* acc = check_malloc(sizeof(float) * header->nt * num_comps);
	float* rotD100 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	float* rotD50 = check_malloc(sizeof(float*) * NUM_INTERP * MAX_PERIODS);
	int* rD100ang = check_malloc(sizeof(float*) * NUM_INTERP * MAX_PERIODS);
	float* acc_for_calc;

	//Added periods for 1 Hz	
	int num_periods = 27;
	float periods[] = {1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10.0, 12.0, 13.0, 15.0, 17.0, 20.0};
    //For PGV
    int num_pgv_periods = 1;
    float pgv_period[] = {1e-5};

	calc_acc(acc, seis_data, num_comps, header->nt, header->dt);
	//If there are 3 components, the vertical component is the first; want to pass just horizontal components to RotD calculation
	if (num_comps==3) {
		acc_for_calc = acc+header->nt;
	} else {
		acc_for_calc = acc;
	}
	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), acc_for_calc, acc_for_calc+header->nt, rotD100, rD100ang, rotD50, &num_periods, periods);
//	write_result(fp_out, header, rotD100, rD100ang, rotD50, num_periods, periods, inter_flag);
	copy_result(rotD_records, rotD100, rD100ang, rotD50, num_periods, periods, inter_flag);

    //Calculate PGV using velocity seismogram
    //Convert to g first
    const float cm2g = 1.0/980.665;
    float* vel_data = check_malloc(sizeof(float)*header->nt*num_comps);
    for (i=0; i<num_comps*header->nt; i++) {
        vel_data[i] = seis_data[i]*cm2g;
    }
	//Also need to only pass horizontal components for velocity calc
	float* vel_for_calc;
	if (num_comps==3) {
		vel_for_calc = vel_data+header->nt;
	} else {
		vel_for_calc = vel_data;
	}
    calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), vel_for_calc, vel_for_calc+header->nt, rotD100+NUM_INTERP*num_periods, rD100ang+NUM_INTERP*num_periods, rotD50+NUM_INTERP*num_periods, &num_pgv_periods, pgv_period);
    //printf("PGV:%f\n", rotD50[NUM_INTERP*num_periods]);

    //Add PGV to output data structure
    rotD_records[num_periods].period = pgv_period[0];
    rotD_records[num_periods].rotD100 = rotD100[num_periods*NUM_INTERP + (inter_flag-1)];
    rotD_records[num_periods].rotD100_angle = rD100ang[num_periods*NUM_INTERP + (inter_flag-1)];
    rotD_records[num_periods].rotD50 = rotD50[num_periods*NUM_INTERP + (inter_flag-1)];

//	fflush(stderr);
//	fflush(stdout);
	free(rotD100);
	free(rotD50);
	free(rD100ang);
	free(acc);
	free(vel_data);

	return 0;
}

int vert_rsp(struct seisheader* header, float* seis_data, struct vertical_rsp_record* rsp_records) {
    int i;
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

	if (num_comps!=3) {
		fprintf(stderr, "The vertical response code must be used with 3-component seismograms.\n");
		return 2;
	}

	float* acc = check_malloc(sizeof(float) * header->nt);
	float* rsp = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	float* rotD100_placeholder = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	int* angle_placeholder = check_malloc(sizeof(int) * NUM_INTERP * MAX_PERIODS);
	float* acc_for_calc;

	//Periods for 1 Hz
	int num_periods = 27;
    float periods[] = {1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10.0, 12.0, 13.0, 15.0, 17.0, 20.0};

	//Just Z-component
	calc_acc(acc, seis_data, 1, header->nt, header->dt);	

	acc_for_calc = acc;

	//Pass Z component twice
	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), acc_for_calc, acc_for_calc, rotD100_placeholder, angle_placeholder, rsp, &num_periods, periods);
	for (i=0; i<num_periods; i++) {
		rsp_records[i].period = periods[i];
		rsp_records[i].rsp = rsp[i*NUM_INTERP + (inter_flag-1)];
	}

	free(acc);
	free(rsp);
	free(rotD100_placeholder);
	free(angle_placeholder);

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
