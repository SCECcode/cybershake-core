#ifndef RUPGEN_FUNC_H
#define RUPGEN_FUNC_H

#include "structure.h"

void *check_malloc(size_t);
void *check_realloc(void *,size_t);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);

void fft2d(struct complex *, int, int, int,float *,float *);
void kfilt(struct complex *,int,int,float *,float *,float *,float *,long *,int);
void init_slip(struct complex *,int,int,float *,float *);
void init_slip_IO(struct complex *,int,int,float *,float *,int,char *);

void taper_slip(struct complex *,int,int,float *,float *);
void taper_slip_all(struct complex *,int,int,float *,float *,float *);

void scale_slip(struct pointsource *,struct complex *,int,int,int,float *,float *,float *,float *,float *,struct velmodel *,float *,float *);

void write_field(char *,struct pointsource *,char *,int,int,float *,float *);
void write_spec(char *,float *,struct complex *,int,int,float *,float *,float *,float *,float *,float *,int);
void write_avgspec(char *,float *,int,int,int,float *,float *);

void default_velmodel(struct velmodel *);
void read_velmodel(char *,struct velmodel *);
void conv2vrup(struct velmodel *,struct velmodel *,float *,float *,float *,float *,float *);
void conv2vrup_dd(struct velmodel *,struct velmodel *,float *,float *,float *,float *,float *,float *,float *);

double frand(void);
double sfrand(long *);
double gaus_rand(float *,float *,long *);

int gen_brune_stf(float *,float *,float *,int,float *);
int gen_ucsb_stf(float *,float *,float *,int,float *);
int gen_tri_stf(float *,float *,float *,int,float *);
int gen_2tri_stf(float *,float *,float *,int,float *,float *);

void set_ll(float *,float *,float *,float *,float *,float *);
void swap_in_place(int,char *);

struct pointsource *read_ruppars(char *,struct pointsource *,float *,int *,int *,float *,float *,float *,float *,float *,float *,float *);
struct pointsource *read_gsfpars(char *,struct pointsource *,struct generic_slip *,float *,float *,float *,float *);
struct pointsource *set_ruppars(struct pointsource *,float *,int *,int *,float *,float *,float *,float *,float *,float *,float *,float *);

void init_plane_srf(struct standrupformat *,struct generic_slip *,float *,float *,int,int,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void load_slip_srf(struct standrupformat *,struct stfpar *,struct pointsource *);
void load_rupt_srf(struct standrupformat *,struct pointsource *,float *,float *);
void write_srf(struct standrupformat *,char *,int);
void write2gsf(struct generic_slip *,struct pointsource *,char *,char *);

void free_srf_stf(struct standrupformat *);

/* Copy file from source to dest */
int cp(const char *to, const char *from);
/* Check if file exists */
int file_exists(const char *file);
/* Check that file exists and is a file */
int rg_is_file(const char *path);
/* Safe string copy */
int rg_strcpy(char *str1, const char *str2, int str1len);

//int fwrite_buffered(FILE *fd, char *pntr, int length);
//int fwrite_flush(FILE *fd);

#endif


