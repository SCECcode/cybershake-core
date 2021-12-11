#ifndef RUPGEN_GENSLIP_H
#define RUPGEN_GENSLIP_H

#include "structure.h"

#ifdef _USE_MEMCACHED
int mc_genslip(int ac,char **av, rg_stats_t *stats, struct standrupformat* srf, int state, char* mc_server);
struct pointsource* _mc_read_ruppars(char *file,struct pointsource *psrc,float *mag,int *nx,int *ny,float *dx,float *dy,float *dtop,float *
stk,float *dip,float *elon,float *elat, char* mc_server);
#else
int genslip(int ac,char **av, rg_stats_t *stats, struct standrupformat* srf, int state);
#endif

void write2srf(struct standrupformat *srf,char *str,int outbin);

#endif
