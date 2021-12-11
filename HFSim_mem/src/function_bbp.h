#include "structure.h"

void *check_malloc(size_t);
void *check_realloc(void *,size_t);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);
float *read_wccseis(char *ifile,struct statdata *shead,float *s,int bflag);
void getheader(char *str,struct statdata *hd);

double frand(void);
double sfrand(long *);
double gaussian_rand(float *,float *,long *);

int gen_esg2006_stf(float *,float *,float *,int,float *,float *);
int gen_2tri_stf(float *,float *,float *,int,float *,float *);
int gen_ucsb_stf(float *,float *,float *,int,float *,float *);
int gen_brune_stf(float *,float *,float *,int,float *,float *);
int gen_cos_stf(float *,float *,float *,int,float *,float *);
int gen_seki_stf(float *,float *,float *,int,float *,float *);

void set_ne(float *,float *,float *,float *,float *,float *);
void set_ll(float *,float *,float *,float *,float *,float *);
void swap_in_place(int,char *);

void init_plane_srf(struct standrupformat *,float *,float *,int,int,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void load_slip_srf(struct standrupformat *,struct stfpar *,struct pointsource *);
void load_rupt_srf(struct standrupformat *,struct pointsource *,float *,float *);

void read_srf(struct standrupformat *,char *,int);
void write_srf(struct standrupformat *,char *,int);

void free_srf_stf(struct standrupformat *);

int write_xyz(char *,struct standrupformat *,char *,int,int,float *,float *,int,int,int,float *,float *,int,int);
void write_maxsvf(char *,struct standrupformat *,char *,int,float *);
void get_vmax2slip(char *,struct standrupformat *,char *,int);
int write_lld(char *,struct standrupformat *,int,float *,float *);

void get_moment(struct standrupformat *,struct velmodel *);
void read_velmodel(char *,struct velmodel *); 

void sum_srf(struct standrupformat *,struct standrupformat *,struct standrupformat *,float *);
void join_srf(struct standrupformat *,struct standrupformat *,struct standrupformat *);

float ucvm_vs30(float lon, float lat, char* model);

int wcc_siteamp09(float* seis, int nt, float dt, float pga, float vs30);

void hfsim(float** seisC, char* stat, float slon, float slat, char* local_vmod, FILE* output_fp, FILE* raw_output_fp, float vs30, struct seisheader* header, float modelrot, struct slipfile* sfile, int num_comps, int do_site_response, float vref, float vpga, int seed, int debug);

int add_param_reset(char** param_string, char* param, int reset);
int add_param(char** param_string, char* param);
