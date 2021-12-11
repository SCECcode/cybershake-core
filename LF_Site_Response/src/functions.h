#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "include.h"
#include "structure.h"

void *check_malloc(size_t len);
void *check_realloc(void *ptr,size_t len);

float ucvm_vs30(float lon, float lat, char* model);
int wcc_siteamp14(float* seis, int nt, float dt, float pga, float vs30, float vref, char model[128]);
void integ_diff(int integ_notdiff, float* seis, int nt, float dt);
float wcc_getpeak(float* seis, int nt, float dt);

int getnt_p2(int nt);

void forfft(register struct complex* x,int n,int isign);
void invfft(register struct complex* x,int n,int isign);

#endif
