#ifndef RUPGEN_GENSLIP_H
#define RUPGEN_GENSLIP_H

#include "structure.h"

void write2srf(struct standrupformat *,char *,int);
int genslip(int ac,char **av, rg_stats_t *stats);

#endif
