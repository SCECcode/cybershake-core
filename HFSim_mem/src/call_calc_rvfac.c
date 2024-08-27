#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>

int main(int argc, char** argv) {
	Py_Initialize();
	/*PyRun_SimpleString("foo=1");
	PyRun_SimpleString("bar=2");
	PyRun_SimpleString("foo += bar");
	PyRun_SimpleString("print(foo)");
	Py_Finalize();*/
	char* filename = "/work2/00349/scottcal/frontera/CyberShake/software/bbp/bbp-22.4.0/bbp/comps/hfsims_cfg.py";
	//char* filename = "test_hfsims_cfg.py";
	FILE* hfsims_cfg_fp = fopen(filename, "r");
	// = fopen("/gpfs/alpine/proj-shared/geo112/CyberShake/software/bbp/bbp-19.8.0-python3/bbp/comps/hfsims_cfg.py", "r")
	PyRun_SimpleFile(hfsims_cfg_fp, filename);
	char call[128];
	double mean_rvfac = 0.1;
	double range_rvfac = 0.2;
	long seed = 2379646;
	PyObject* main_module = PyImport_AddModule("__main__");
	if (main_module==NULL) {
		printf("Error getting main module, aborting.\n");
		exit(2);
	}
	PyObject* global_dict = PyModule_GetDict(main_module);

	PyObject* expression = PyDict_GetItemString(global_dict, "calculate_rvfac");

	PyObject* args = PyTuple_New(3);
	if (PyTuple_SetItem(args, 0, PyFloat_FromDouble(mean_rvfac))!=0) {
		printf("Error adding mean_rvfac, aborting.\n");
		exit(1);
	}
	PyTuple_SetItem(args, 1, PyFloat_FromDouble(range_rvfac));
	PyTuple_SetItem(args, 2, PyLong_FromLong(seed));

	printf("foo\n");
	fflush(stdout);
	PyObject* py_rvfac = PyObject_CallObject(expression, args);
	//PyObject rvfac = PyObject_CallFunction(expression, "ffi", mean_rvfac, range_rvfac, seed);
	double rvfac = PyFloat_AsDouble(py_rvfac);
	printf("rvfac = %lf\n", rvfac);
	Py_Finalize();
	return 0;
}
