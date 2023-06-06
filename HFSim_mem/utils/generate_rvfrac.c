#include <stdio.h>
#include <stdlib.h>
#include <Python.h>

#include "rupgen_api.h"

float get_rvfac(double mean_rvfac, double range_rvfac, int seed) {
        //char* filename = "/gpfs/alpine/proj-shared/geo112/CyberShake/software/bbp/bbp-19.8.0-python3/bbp/comps/hfsims_cfg.py";
        char* filename = "/gpfs/alpine/scratch/callag/geo112/bbp/bbp/comps/hfsims_cfg.py";
        FILE* hfsims_cfg_fp = fopen(filename, "r");
        Py_Initialize();
        PyRun_SimpleFile(hfsims_cfg_fp, filename);
        PyObject* main_module = PyImport_AddModule("__main__");
        PyObject* global_dict = PyModule_GetDict(main_module);
        PyObject* expression = PyDict_GetItemString(global_dict, "calculate_rvfac");
        PyObject* args = PyTuple_New(3);
        PyTuple_SetItem(args, 0, PyFloat_FromDouble(mean_rvfac));
        PyTuple_SetItem(args, 1, PyFloat_FromDouble(range_rvfac));
        PyTuple_SetItem(args, 2, PyLong_FromLong((long)seed));
        PyObject* py_rvfac = PyObject_CallObject(expression, args);

        double rvfac = PyFloat_AsDouble(py_rvfac);
        Py_DECREF(args);
        Py_Finalize();

        fclose(hfsims_cfg_fp);
        return (float)rvfac;
}


int main(int argc, char** argv) {
	if (argc<3) {
		printf("Usage: %s <rv input file> <rvfac output file> [-r]\n", argv[0]);
		return -1;
	}

	char* RUPTURE_ROOT = "/gpfs/alpine/proj-shared/geo112/CyberShake/ruptures/Ruptures_erf36";

	char* input_file = argv[1];
	char* output_file = argv[2];
	int restart = 0;
	if (argc>3) {
		if (strcmp(argv[3], "-r")==0) {
			restart = 1;
		}
	}

	int src_id, rup_id, rv_id;
	int src_restart, rup_restart, rv_restart;
    double mean_rvfac = 0.775;
    double range_rvfac = 0.1;

	char rup_geom_file[256];
	int variation_seed;
	float rvfrac;
	rg_stats_t stats;

	FILE* fp_in, *fp_out;
	if (restart) {
		fp_in = fopen(output_file, "r");
		while (fscanf(fp_in, "%d %d %d %f", &src_id, &rup_id, &rv_id, &rvfrac)!=EOF);
		printf("Last entry: %d %d %d\n", src_id, rup_id, rv_id);
		src_restart = src_id;
		rup_restart = rup_id;
		rv_restart = rv_id;
		fclose(fp_in);
	}

	fp_in = fopen(input_file, "r");
	fp_out = fopen(output_file, "a");
	int i=0;
	while (fscanf(fp_in, "%d %d %d", &src_id, &rup_id, &rv_id)!=EOF) {
		if (restart) {
			if (src_id<src_restart) {
				continue;
			} else if (src_id==src_restart) {
				if (rup_id<rup_restart) {
					continue;
				} else if (rup_id==rup_restart) {
					if (rv_id<=rv_restart) {
						continue;
					}
				}
			}
		}
		sprintf(rup_geom_file, "%s/%d/%d/%d_%d.txt", RUPTURE_ROOT, src_id, rup_id, src_id, rup_id);
		variation_seed = rupgen_get_variation_seed(rup_geom_file, &stats, RUPGEN_UNIFORM_HYPO, src_id, rup_id, rv_id, 0);
		rvfrac = get_rvfac(mean_rvfac,range_rvfac,variation_seed);
		fprintf(fp_out, "%d %d %d %f\n", src_id, rup_id, rv_id, rvfrac);
		fflush(fp_out);
		i++;
	}
	fclose(fp_in);
	fflush(fp_out);
	fclose(fp_out);

	return 0;
}
