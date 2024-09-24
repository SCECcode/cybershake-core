#include "include.h"
#include "structure.h"
#include "function.h"
#include <Python.h>

float duration_periods[66] = {0.01, 0.011, 0.012, 0.013, 0.015, 0.017, 0.02, 0.022, 0.025, 0.029, 0.032, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.075, 0.085, 0.1, 0.11, 0.12, 0.13, 0.15, 0.17, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.85, 1, 1.1, 1.2, 1.3, 1.5, 1.7, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 4.4, 5, 5.5, 6, 6.5, 7.5, 8.5, 10, 12, 15, 20};
int num_duration_periods = 66;
char* duration_strings[3] = {"5-75", "5-95", "20-80"};
int duration_type_vals[3] = {D5_75, D5_95, D20_80};
int num_period_duration_types = 3;
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
	//Determine number of components
	int comp_int = header.comps;
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
	
	for (i=0; i<num_comps; i++) {
		entry.component = i;

	    if (i==X_COMP) {
            if (num_comps==2) {
                seis = full_seis[0];
			} else { //3 comps
                //Z is first
               seis = full_seis[1];
			}
        } else if (i==Y_COMP) {
            if (num_comps==2) {
                seis = full_seis[1];
            } else {
                seis = full_seis[2];
            }
        } else if (i==Z_COMP) {
            seis = full_seis[0];
        }

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

//This function wraps period_duration.py to perform period-dependent duration calculations
int period_durations(struct seisheader header, float** full_seis, FILE* fp_out) {
    float* seis;
    struct period_duration_record entry;
    int nt = header.nt;
    float dt = header.dt;
    int comp_int = header.comps;
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
	//Write header info
	fwrite(&header, sizeof(struct seisheader), 1, fp_out);
	//Number of entries is 2 components x 3 measures x all the periods
	int num_entries = 2*num_period_duration_types*num_duration_periods;
	fwrite(&num_entries, sizeof(int), 1, fp_out);

    //Python setup
    if (python_initialized==0) {
		//If this is run inside a Condor-submitted job, the HOME environment variable is reset to a different directory, and this leads to errors on PyRun_SimpleFile.  In case of this, we reset HOME to the correct location.
		setenv("HOME", "/home1/00349/scottcal", 1);
        Py_Initialize();
        char* filename = "/work2/00349/scottcal/frontera/CyberShake/software/MergeIM/src/duration/period_duration.py";
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

    PyObject* calc_dur_args = PyTuple_New(3);
    //times is constant for all calls
    PyObject* times_list = PyList_New(nt);
    for (i=0; i<nt; i++) {
        PyList_SetItem(times_list, i, PyFloat_FromDouble(i*dt));
    }
    PyTuple_SetItem(construct_rec_args, 1, times_list);

    for (i=0; i<2; i++) {
		entry.component = i;
        if (i==X_COMP) {
            if (num_comps==2) {
                seis = full_seis[0];
            } else if (num_comps==3) {
                //Z is first, X is 2nd
                seis = full_seis[1];
            } else {
                printf("Can't perform period-dependent duration calculations on single-component seismograms.\n");
                return -1;
            }
        } else if (i==Y_COMP) {
            if (num_comps==2) {
                seis = full_seis[1];
            } else if (num_comps==3) {
                //Z, X, Y
                seis = full_seis[2];
            }
        }
        //Create Python list from seis float array
        PyObject* seis_list = PyList_New(nt);
        int j;
        for (j=0; j<nt; j++) {
            PyList_SetItem(seis_list, j, PyFloat_FromDouble(seis[j]));
        }
        PyTuple_SetItem(construct_rec_args, 0, seis_list);
        PyObject* rec = PyObject_CallObject(construct_rec_expression, construct_rec_args);
        PyTuple_SetItem(calc_dur_args, 2, rec);
        for (j=0; j<num_period_duration_types; j++) {
            PyTuple_SetItem(calc_dur_args, 1, PyUnicode_FromString(duration_strings[j]));
            int k;
            for (k=0; k<num_duration_periods; k++) {
                PyObject* freq = PyFloat_FromDouble(1.0/duration_periods[k]);
                PyTuple_SetItem(calc_dur_args, 0, freq);
                PyObject* dur_val_item = PyObject_CallObject(calc_dur_expression, calc_dur_args);
                float dur_val = (float)(PyFloat_AsDouble(dur_val_item));

                entry.type = DV;
                entry.type_value = duration_type_vals[j];
                entry.component = i;
                entry.period = duration_periods[k];
                entry.value = dur_val;
				fwrite(&entry, sizeof(struct period_duration_record), 1, fp_out);
            }
        }
    }
    Py_DECREF(construct_rec_args);
    Py_DECREF(calc_dur_args);
    Py_DECREF(times_list);
    Py_DECREF(PyImport_ImportModule("threading"));

	fflush(fp_out);

	return 0;
}

void python_finalize() {
    if (python_initialized==1) {
        fclose(per_dur_fp);
        Py_Finalize();
    }
}

