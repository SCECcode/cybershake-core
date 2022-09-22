#include "structure.h"
#include "StandRupFormat/function.h"
#include "WccFormat/Progs/function.h"

float ucvm_vs30(float lon, float lat, char* model);

void hfsim(float** seisC, char* stat, float slon, float slat, char* local_vmod, FILE* output_fp, FILE* raw_output_fp, float vs30, struct seisheader* header, float modelrot, struct slipfile* sfile, int num_comps, int do_site_response, float vref, float vpga, int seed, float rvfrac, int debug);

int add_param_reset(char** param_string, char* param, int reset);
int add_param(char** param_string, char* param);

