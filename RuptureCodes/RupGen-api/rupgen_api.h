#ifndef RUPGEN_API_H
#define RUPGEN_API_H

/* Rupture generation statistics */
typedef struct rg_stats_t {
  int numslip;
  int numhypo;
} rg_stats_t;

int rupgen_get_num_rv(char* rup_geom_file, rg_stats_t *stats);
int rupgen_genslip(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf);

#endif
