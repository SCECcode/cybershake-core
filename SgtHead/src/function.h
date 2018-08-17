/* Started 07/07/04- attempt to go towards ANSI C */


void swap_in_place(int,char *);
void swap_in_place_d(int,char *);
void set_tmpdir(char *,int,char *,char *);

void check_procname(char *name,int *len);
void get_ny1ny2(int,int,int,int,int,int *,int *,int *);

void init_field_val(float,float *,int);
void floor_c(float *,float *,int);

void *check_malloc(size_t);
void *check_realloc(void *,size_t);

void resample(float *,int,float *,int,int,int,float *,float *,int,float *,float *);
void taper_norm(float *,float *,int,float *);
void norm(float *,float *,int);
double nt_tol(float,int);

void get_dc_par(struct pntsrcs *,float *,float *,int);
void get_pointmt_par(struct pntsrcs *,float *,float *,int);
void get_ff_par(struct pntsrcs *,float *,float *,int,float *);
void get_srf_par(struct pntsrcs *,struct runparams *,int,float *,float *,float *);

void init_outp(struct outputfields *,struct runparams *,char *,struct pntsrcs *,struct restartinfo *,struct dump_output_info *);
void check_rpars(char *,struct runparams *,struct runparams *);
void dump_restart(struct restartinfo *,char *,int,struct runparams *,struct modelstorage *,float *,int,int,int,int);
void get_restart(struct restartinfo *,struct modelstorage *,float *,int,int,int);
void print_rpars(struct runparams *,struct runparams *);
void dump_files(struct outputfields *,int,struct dump_output_info *);
void rm_dump_files(struct outputfields *,int,struct dump_output_info *);

float *reed_model(struct modelstorage *,float *,int);

int eff_media(struct media_input *,float *,int,int,int,int,int,struct modelstorage *,int);
void init_media2(struct media_input *,float *,int,int,int,int,int,int,float *,struct qvalues *,struct modelstorage *,int,int);
void init_media(struct media_input *,float *,struct runparams *,struct qvalues *,struct modelstorage *,int);

void reedmedslice(int,int,int,float *,int,int,int,int,struct qvalues *,int,int,float *,float *,int);

size_t chunk_reed(int fd, void* pntr, size_t length);
size_t reed(int,void *,size_t);
size_t rite(int,void *,size_t);
int cropfile_rw(char *);
int fcp_read_write(char *,char *);
int opfile(char*);
int opfile_ro(char*);

void write_mt(struct pntsrcs *,float *,float *,int,int,int);

void init_modelc(struct runparams *);

void gcproj(float *,float *,float *,float *,float *,double *,double *,double *,double *,int);
void gen_matrices(double *,double *,float *,float *,float *);

void latlon2km(float *,float *,float *,float *,float *);
void set_g2(float *,float *);
void geocen(float *,double);

void check_nan(float *p,float *m,int nx,int nz,int iy,int it,int j);

void init_dc(struct pntsrcs *,float *);
void init_pointmt(struct pntsrcs *,float *,float *);
void init_mt(struct pntsrcs *,float *,struct modelstorage *,float *,int,int,int,int,int,char *,float *,char *,char *,int);

int tshift_st(float* s,float* dt,int nt,float* flo,float* tsh);

FILE *fopfile(char* name,char* mode);

