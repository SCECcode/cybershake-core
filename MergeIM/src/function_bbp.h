#include "WccFormat/Progs/function.h"

void merge_data(float*** lf_seis, struct seisheader lf_header, float** hf_seis, struct seisheader hf_header, float match_freq, int num_comps, float** merged_seis, struct seisheader* merged_header, int debug);

int rotd(struct seisheader* header, float* seis_data, FILE* fp_out);

int duration(struct seisheader header, float** full_seis, FILE* fp_out);

int add_param_reset(char** param_string, char* param, int reset);
int add_param(char** param_string, char* param);
