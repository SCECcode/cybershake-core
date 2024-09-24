#include "../include.h"
#include "../defs.h"
#include "duration.h"
#include "../structure.h"
#include "../functions.h"
#include <Python.h>

float duration_periods[27] = {1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.4, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10.0, 12.0, 13.0, 15.0, 17.0, 20.0};
char* duration_strings[3] = {"5-75", "5-95", "20-80"};
int duration_type_vals[3] = {D5_75, D5_95, D20_80};
//To prevent repeated initializations of the interpreter
int python_initialized = 0;
PyObject* construct_rec_expression;
PyObject* calc_dur_expression;
FILE* per_dur_fp;

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
		float ia = cumulative_integral[nt-1]*3.14/(2*981.0);
		free(cumulative_integral);
		return ia;
	}
	//Ie
	float ie = cumulative_integral[nt-1];
	free(cumulative_integral);
	return ie;
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
	free(cumulative_integral);
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

int duration(struct seisheader* header, float* full_seis, struct duration_record* out) {
	float* seis;
	struct duration_record* entry;	
	int i;
	//fwrite(&header, sizeof(struct seisheader), 1, fp_out);
	//number of measures per rupture variation, per component
	//fwrite(&NUM_DURATION_MEASURES, sizeof(int), 1, fp_out);
	int nt = header->nt;
	float dt = header->dt;
	//Determine number of components
	int comp_int = header->comps;
	int num_comps = 0;
	if (( comp_int & 1)!=0) {
		num_comps++;
	}
	if (( comp_int & 2)!=0) {
		num_comps++;
	}
	if (( comp_int & 4)!=0) {
		num_comps++;
	}

	//printf("nt = %d\n", header->nt);
	//printf("out addr: %p\n", out);
	//Save accel seismogram to debug
	//char acc_filename[256];
	//sprintf(acc_filename, "ACC_Seis_%d_%d_%d.grm", header->source_id, header->rupture_id, header->rup_var_id);
	//FILE* acc_out = fopen(acc_filename, "wb");
	for (i=0; i<num_comps; i++) {
		entry = out+i*NUM_DURATION_MEASURES;
                //printf("starting entry addr: %p\n", entry);
        if (i==X_COMP) {
			if (num_comps==2) {
				seis = full_seis;
			} else { //3 comps
				//Z is first
				seis = full_seis + nt;
			}
		} else if (i==Y_COMP) {
			if (num_comps==2) {
				seis = full_seis + nt;
			} else {
				seis = full_seis + 2*nt;
			}
		} else if (i==Z_COMP) {
			seis = full_seis;
		}

		float energy_integral = intensity(seis, nt, dt, 0);
		entry->component = i;
		entry->type = ENERGY_INTEGRAL;
		entry->type_value = -1;
		entry->value = energy_integral;
		//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		float dv_5_75 = significant_duration(seis, 0.05, 0.75, nt, dt, 0);
		entry += 1;
		entry->component = i;
		entry->type = DV;
		entry->type_value = D5_75;
		entry->value = dv_5_75;
		//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		float dv_5_95 = significant_duration(seis, 0.05, 0.95, nt, dt, 0);
		entry += 1;
		entry->component = i;
        	entry->type = DV;
        	entry->type_value = D5_95;
        	entry->value = dv_5_95;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

        	float dv_20_80 = significant_duration(seis, 0.20, 0.80, nt, dt, 0);
		entry += 1;
		entry->component = i;
        	entry->type = DV;
       		entry->type_value = D20_80;
        	entry->value = dv_20_80;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

		integ_diff(0, seis, nt, dt);
		//fwrite(seis, sizeof(float), nt, acc_out);	

		float arias_intensity = intensity(seis, nt, dt, 1);
		entry += 1;
		entry->component = i;
        	entry->type = ARIAS_INTENSITY;
        	entry->type_value = -1;
        	entry->value = arias_intensity;
		//printf("arias intensity value in entry: %f at addr %p\n", entry->value, entry);
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);

        	float da_5_75 = significant_duration(seis, 0.05, 0.75, nt, dt, 1);
		entry += 1;
		entry->component = i;
        	entry->type = DA;
        	entry->type_value = D5_75;
        	entry->value = da_5_75;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float da_5_95 = significant_duration(seis, 0.05, 0.95, nt, dt, 1);
		entry += 1;
		entry->component = i;
        	entry->type = DA;
        	entry->type_value = D5_95;
        	entry->value = da_5_95;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float da_20_80 = significant_duration(seis, 0.20, 0.80, nt, dt, 1);
		entry += 1;
		entry->component = i;
        	entry->type = DA;
        	entry->type_value = D20_80;
        	entry->value = da_20_80;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	
        	float cav = calc_cav(seis, nt, dt);
		entry += 1;
		entry->component = i;
		entry->type = CAV;
		entry->type_value = -1;
		entry->value = cav;
        	//fwrite(&entry, sizeof(struct duration_record), 1, fp_out);
	}
	//fflush(fp_out);
	//fflush(acc_out);
	//fclose(acc_out);
	return 0;
}

//This function wraps period_duration.py to perform period-dependent duration calculations
int period_durations(struct seisheader* header, float* full_seis, struct period_duration_record* out) {
	float* seis;
	struct period_duration_record* entry;
    int nt = header->nt;
    float dt = header->dt;
    int comp_int = header->comps;
    int num_comps = 0;
    if (( comp_int & 1)!=0) {
        num_comps++;
    }
    if (( comp_int & 2)!=0) {
        num_comps++;
    }
    if (( comp_int & 4)!=0) {
        num_comps++;
    }
	int i;
	//Python setup
	if (python_initialized==0) {
		//If this is run inside a Condor-submitted job, the HOME environment variable is reset to a different directory, and this leads to errors on PyRun_SimpleFile.  In case of this, we reset HOME to the correct location.
		setenv("HOME", "/home1/00349/scottcal", 1);
		Py_Initialize();
		char* filename = "/work2/00349/scottcal/frontera/CyberShake/software/DirectSynth/src/duration/period_duration.py";
		per_dur_fp = fopen(filename, "r");
		int rc = PyRun_SimpleFile(per_dur_fp, filename);
		if (rc!=0) {
			printf("Error with SimpleFile().\n");
			return -1;
		}
		PyObject* moduleString = PyUnicode_FromString("period_duration");
		if (moduleString==NULL) {
			printf("Error with moduleString.\n");
			return 1;
		}
    	PyObject* module = PyImport_Import(moduleString);
		if (module==NULL) {
			printf("Error with module.\n");
			return 2;
		}
		//def construct_rec(seis, times):
		construct_rec_expression = PyObject_GetAttrString(module, "construct_rec");
		if (construct_rec_expression==NULL) {
			printf("Error finding construct_rec.\n");
			return 3;
		}
	    calc_dur_expression = PyObject_GetAttrString(module, "c_calc_rec_durations");
	    if (calc_dur_expression==NULL) {
	        printf("Error finding c_calc_rec_durations.\n");
	        return 4;
	    }
		python_initialized = 1;
	}
	PyObject* construct_rec_args = PyTuple_New(2);
	
	/*PyObject* calc_dur_expression = PyObject_GetAttrString(module, "c_calc_rec_durations");
	if (calc_dur_expression==NULL) {
		printf("Error finding c_calc_rec_durations.\n");
		return 4;
	}*/
	//c_calc_rec_durations(freq, dur, rec):
	PyObject* calc_dur_args = PyTuple_New(3);
	
	//times is constant for all calls
	PyObject* times_list = PyList_New(nt);
	for (i=0; i<nt; i++) {
		PyList_SetItem(times_list, i, PyFloat_FromDouble(i*dt));
	}
	PyTuple_SetItem(construct_rec_args, 1, times_list);

	//Only do period-dependent duration calcs for horizontal components
	//PyObject* seis_list = PyList_New(nt);
	for (i=0; i<2; i++) {
		//printf("Beginning calcs for component %d.\n", i);
		//fflush(stdout);
        entry = out+i*NUM_DURATION_PERIODS*NUM_PERIOD_DURATION_MEASURES;
		if (i==X_COMP) {
			if (num_comps==2) {
				seis = full_seis;
			} else if (num_comps==3) {
				//Z is first, X is 2nd
				seis = full_seis + nt;
			} else {
				printf("Can't perform period-dependent duration calculations on single-component seismograms.\n");
				return -1;
			}
		} else if (i==Y_COMP) {
			if (num_comps==2) {
				seis = full_seis + nt;
			} else if (num_comps==3) {
				//Z, X, Y
				seis = full_seis + 2*nt;
			}
		}
		//Create Python list from seis float array
		PyObject* seis_list = PyList_New(nt);
		int j;
		for (j=0; j<nt; j++) {
			PyList_SetItem(seis_list, j, PyFloat_FromDouble(seis[j]));
		}
		PyTuple_SetItem(construct_rec_args, 0, seis_list);
		//printf("Constructing rec object for src %d rup %d rv %d comp %d\n", header->source_id, header->rupture_id, header->rup_var_id, i);
		//fflush(stdout);
		PyObject* rec = PyObject_CallObject(construct_rec_expression, construct_rec_args);
		PyTuple_SetItem(calc_dur_args, 2, rec);
		//Perform 5-75, 5-95, 20-80 on velocity
		for (j=0; j<NUM_PERIOD_DURATION_MEASURES; j++) {
			PyTuple_SetItem(calc_dur_args, 1, PyUnicode_FromString(duration_strings[j]));
			int k;
			for (k=0; k<NUM_DURATION_PERIODS; k++) {
				PyObject* freq = PyFloat_FromDouble(1.0/duration_periods[k]);
				PyTuple_SetItem(calc_dur_args, 0, freq);
				//printf("Running duration for src %d rup %d rv %d comp %d dur %s period %f\n", header->source_id, header->rupture_id, header->rup_var_id, i, duration_strings[j], duration_periods[k]);
				//fflush(stdout);
				PyObject* dur_val_item = PyObject_CallObject(calc_dur_expression, calc_dur_args);
				float dur_val = (float)(PyFloat_AsDouble(dur_val_item));
				//printf("Copying to entry at %p\n", entry);
				//fflush(stdout);
				entry->type = DV;
				entry->type_value = duration_type_vals[j];
				entry->component = i;
				entry->period = duration_periods[k];
				entry->value = dur_val;
				entry += 1;
			}
		}	
	}
    Py_DECREF(construct_rec_args);
    Py_DECREF(calc_dur_args);
    Py_DECREF(times_list);
	//Py_DECREF(seis_list);
    Py_DECREF(PyImport_ImportModule("threading"));
    //fclose(per_dur_fp);

	//Py_Finalize();
	//printf("Completed period duration calculation for src %d rup %d rv %d\n", header->source_id, header->rupture_id, header->rup_var_id);
	return 0;
}

void python_finalize() {
	if (python_initialized==1) {
		fclose(per_dur_fp);
		Py_Finalize();
	}
}
