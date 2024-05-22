#include "include.h"
#include "structure.h"
#include "function.h"

void parse_seis(char* seis_filename, float*** seis, struct seisheader* header) {
	FILE* fp_in = fopen(seis_filename, "rb");
	fread(header, sizeof(struct seisheader), 1, fp_in);
	*seis = malloc(sizeof(float*)*2);
	int i;
	for (i=0; i<2; i++) {
		(*seis)[i] = malloc(sizeof(float)*header->nt);
		fread((*seis)[i], sizeof(float), header->nt, fp_in);
	}
	fclose(fp_in);
}

float intensity(float* seis, int nt, float dt, int acc_flag) {
	//Calculates Arias intensity with acceleration or energy integral with velocity
	int i;
	float* cumulative_integral = malloc(sizeof(float)*nt);
        cumulative_integral[0] = 0.0;
        for (i=1; i<nt; i++) {
                cumulative_integral[i] = cumulative_integral[i-1] + 0.5*(seis[i]*seis[i] + seis[i-1]*seis[i-1])*dt;
        }
	if (acc_flag) {
		//Ia
		return cumulative_integral[nt-1]*3.14/(2*981.0);
	}
	//Ie
	return cumulative_integral[nt-1];
}


float significant_duration(float* seis, float low, float high, int nt, float dt, int acc_flag) {
	//Calculate v(t)^2 using trapezoid rule
	int i;
	//Start with just X
	float* cumulative_integral = malloc(sizeof(float)*nt);
	cumulative_integral[0] = 0.0;
	for (i=1; i<nt; i++) {
		cumulative_integral[i] = cumulative_integral[i-1] + 0.5*(seis[i]*seis[i] + seis[i-1]*seis[i-1])*dt;
	}
	if (!acc_flag) {
		printf("Ie = %f\n", cumulative_integral[nt-1]);
	} else {
		printf("Ia = %f\n", cumulative_integral[nt-1]*3.14/(2*981.0));
	}
	float low_bound = low*cumulative_integral[nt-1];
	float high_bound = high*cumulative_integral[nt-1];
	//printf("Searching for low_bound %f, high_bound %f\n", low_bound, high_bound);
	int low_bound_ts, high_bound_ts;
	for (i=0; i<nt; i++) {
		if (cumulative_integral[i]>low_bound) {
			low_bound_ts = i;
			break;
		}
	}
	for (i=0; i<nt; i++) {
		if (cumulative_integral[i]>high_bound) {
			high_bound_ts = i;
			break;
		}
	}
	float duration = (high_bound_ts - low_bound_ts)*dt;
	return duration;
}

float calc_cav(float* seis, int nt, float dt) {
        int i;
        float integral = 0.0;
        for (i=1; i<nt; i++) {
                integral += 0.5*(fabs(seis[i])+fabs(seis[i-1]))*dt;
        }
        return integral;
}

int duration(struct seisheader header, float** full_seis, FILE* fp_out) {
	float* seis;
	struct duration_record entry;	
	int i;
	fwrite(&header, sizeof(struct seisheader), 1, fp_out);
	//number of measures per rupture variation, per component
	int NUMBER_OF_MEASURES = 9;
	fwrite(&NUMBER_OF_MEASURES, sizeof(int), 1, fp_out);
	int nt = header.nt;
	float dt = header.dt;
	for (i=0; i<2; i++) {
		entry.component = i;
		seis = full_seis[i];

		float energy_integral = intensity(seis, nt, dt, 0);
		entry.type = ENERGY_INTEGRAL;
		entry.type_value = -1;
		entry.value = energy_integral;
		fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		float dv_5_75 = significant_duration(seis, 0.05, 0.75, nt, dt, 0);
		entry.type = DV;
		entry.type_value = D5_75;
		entry.value = dv_5_75;
		fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		float dv_5_95 = significant_duration(seis, 0.05, 0.95, nt, dt, 0);
        	entry.type = DV;
        	entry.type_value = D5_95;
        	entry.value = dv_5_95;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

        	float dv_20_80 = significant_duration(seis, 0.20, 0.80, nt, dt, 0);
        	entry.type = DV;
       		entry.type_value = D20_80;
        	entry.value = dv_20_80;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		integ_diff(0, seis, nt, dt);

		float arias_intensity = intensity(seis, nt, dt, 1);

        	entry.type = ARIAS_INTENSITY;
        	entry.type_value = -1;
        	entry.value = arias_intensity;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

        	float da_5_75 = significant_duration(seis, 0.05, 0.75, nt, dt, 1);
        	entry.type = DA;
        	entry.type_value = D5_75;
        	entry.value = da_5_75;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float da_5_95 = significant_duration(seis, 0.05, 0.95, nt, dt, 1);
        	entry.type = DA;
        	entry.type_value = D5_95;
        	entry.value = da_5_95;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float da_20_80 = significant_duration(seis, 0.20, 0.80, nt, dt, 1);
        	entry.type = DA;
        	entry.type_value = D20_80;
        	entry.value = da_20_80;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float cav = calc_cav(seis, nt, dt);
		entry.type = CAV;
		entry.type_value = -1;
		entry.value = cav;
        	fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	}
	fflush(fp_out);
	return 0;
}
