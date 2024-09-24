#include "../include.h"
#include "../structure.h"
#include "../function.h"

#define NUM_INTERP 3
#define MAX_PERIODS 200
#define MAXPTS 2000000

#define PGA 1.0e-4
#define PGV 1.0e-5

void write_result(FILE* fp_out, struct seisheader* seis_head, float* rotD100, int* rotD100ang, float* rotD50, int num_periods, float periods[], int inter_flag);
void calc_acc(float* acc, float* seis_data, int num_comps, int nt, float dt);

extern void calc_rotd_(int* inter_flag, int* npairs, int* npts, float* dt, float* acc1, float* acc2, float* rotD100, int* rD100ang, float* rotD50, int* nFreq, float* RSP_Period);

extern const int X_COMP_FLAG;
extern const int Y_COMP_FLAG;
extern const int Z_COMP_FLAG;

int rotd(struct seisheader* header, float* seis_data, FILE* fp_out) {
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
	float* rotD50 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	int* rD100ang = check_malloc(sizeof(int) * NUM_INTERP * MAX_PERIODS);
	
	//Period list taken from BBP
	int num_periods = 69;
	float periods[] = {PGA, 0.01, 0.011, 0.012, 0.013, 0.015, 0.017, 0.02, 0.022, 0.025, 0.029, 0.032, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.075, 0.085, 0.1, 0.11, 0.12, 0.13, 0.15, 0.17, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.85, 1, 1.1, 1.2, 1.3, 1.5, 1.7, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.4, 5, 5.5, 6, 6.5, 7.5, 8.5, 10, 12, 13, 15, 17, 20, PGV};
	//int num_periods = 30;
	//float periods[] = {0.1, 0.125, 0.1666667, 0.2, 0.25, 0.3333333, 0.5, 0.6666667, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10.0};

	calc_acc(acc, seis_data, num_comps, header->nt, header->dt);
	float* acc_for_calc;
    //If there are 3 components, the vertical component is the first; want to pass just horizontal components to RotD calculation
    if (num_comps==3) {
    	acc_for_calc = acc+header->nt;
    } else {
        acc_for_calc = acc;
    }

	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), acc_for_calc, acc_for_calc+header->nt, rotD100, rD100ang, rotD50, &num_periods, periods);

	//Also calculate PGV by using velocity seismogram
	int num_v_periods = 1;
	float v_periods[] = {PGV};
	float* v_rotD100 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
    float* v_rotD50 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
	int* v_rD100ang = check_malloc(sizeof(int) * NUM_INTERP * MAX_PERIODS);
	//Convert to g
	const float cm2g = 1.0/980.665;
    float* vel_data = check_malloc(sizeof(float)*header->nt*num_comps);
	for (i=0; i<num_comps*header->nt; i++) {
        vel_data[i] = seis_data[i]*cm2g;
		//printf("%e\n", vel_data[i]);
    }

    //Also need to only pass horizontal components for velocity calc
    float* vel_for_calc;
    if (num_comps==3) {
        vel_for_calc = vel_data+header->nt;
    } else {
        vel_for_calc = vel_data;
    }

	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), vel_for_calc, vel_for_calc+header->nt, v_rotD100, v_rD100ang, v_rotD50, &num_v_periods, v_periods);
	//Append PGV value to acceleration arrays

	rotD100[num_periods*NUM_INTERP + (inter_flag-1)] = v_rotD100[0];
	rD100ang[num_periods*NUM_INTERP + (inter_flag-1)] = v_rD100ang[0];
	rotD50[num_periods*NUM_INTERP + (inter_flag-1)] = v_rotD50[0];
	num_periods++;

	write_result(fp_out, header, rotD100, rD100ang, rotD50, num_periods, periods, inter_flag);

	//fflush(stderr);
	//fflush(stdout);
	free(rotD100);
	free(rotD50);
	free(rD100ang);
	free(v_rotD100);
	free(v_rotD50);
	free(v_rD100ang);
	free(acc);
	free(vel_data);

	return 0;
}

int vert_rsp(struct seisheader* header, float* seis_data, FILE* fp_out) {
    int i;
    int inter_flag = 1;
    int npairs = 1;

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

	int num_periods = 69; //Not including PGV
    float periods[] = {PGA, 0.01, 0.011, 0.012, 0.013, 0.015, 0.017, 0.02, 0.022, 0.025, 0.029, 0.032, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.075, 0.085, 0.1, 0.11, 0.12, 0.13, 0.15, 0.17, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.85, 1, 1.1, 1.2, 1.3, 1.5, 1.7, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.4, 5, 5.5, 6, 6.5, 7.5, 8.5, 10, 12, 13, 15, 17, 20, PGV};

	calc_acc(acc, seis_data, 1, header->nt, header->dt);
	acc_for_calc = acc;

	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), acc_for_calc, acc_for_calc, rotD100_placeholder, angle_placeholder, rsp, &num_periods, periods);
	
    int num_v_periods = 1;
    float v_periods[] = {PGV};
    float* v_rotD100 = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
    float* v_rsp = check_malloc(sizeof(float) * NUM_INTERP * MAX_PERIODS);
    int* v_rD100ang = check_malloc(sizeof(int) * NUM_INTERP * MAX_PERIODS);

    const float cm2g = 1.0/980.665;
    float* vel_data = check_malloc(sizeof(float)*header->nt*num_comps);
    for (i=0; i<num_comps*header->nt; i++) {
        vel_data[i] = seis_data[i]*cm2g;
    }

	float* vel_for_calc = vel_data;

	calc_rotd_(&inter_flag, &npairs, &(header->nt), &(header->dt), vel_for_calc, vel_for_calc+header->nt, v_rotD100, v_rD100ang, v_rsp, &num_v_periods, v_periods);

	//Append PGV value to acceleration arrays
	rsp[num_periods] = v_rsp[0];
	num_periods++;
	
	fwrite(header, sizeof(struct seisheader), 1, fp_out);
	fwrite(&num_periods, sizeof(int), 1, fp_out);
	struct vertical_rsp_record rsp_record;
	for (i=0; i<num_periods; i++) {
		rsp_record.period = periods[i];
		rsp_record.rsp = rsp[i*NUM_INTERP + (inter_flag-1)];
		fwrite(&rsp_record, sizeof(struct vertical_rsp_record), 1, fp_out);
	}

    free(rotD100_placeholder);
    free(rsp);
    free(angle_placeholder);
    free(v_rotD100);
    free(v_rsp);
    free(v_rD100ang);
    free(acc);
	free(vel_data);

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
