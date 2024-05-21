void *check_malloc(size_t);
void *check_realloc(void *, size_t);
float *read_wccseis(char *, struct statdata *, float *, int);
void write_wccseis(char *, struct statdata *, float *, int);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);
void getheader(char *,struct statdata *);

void swap_in_place(int,char *);

void gcproj(float *,float *,float *,float *,float *,double *,double *,double *,double *,int);
void gen_matrices(double *,double *,float *,float *,float *);

void latlon2km(float *,float *,float *,float *,float *);
void set_g2(float *,float *);
void geocen(float *,double);

void cfft_r(struct complex* x,int n,int isign);
