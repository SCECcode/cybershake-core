#ifndef RUPGEN_API_H
#define RUPGEN_API_H

#include "structure.h"

#ifdef _USE_MEMCACHED
int set_memcached_server(char* memcached_server);
#endif

#define RUPGEN_RANDOM_HYPO 0
#define RUPGEN_UNIFORM_HYPO 1

/* Rupture generation statistics */
typedef struct rg_stats_t {
  int numslip;
  int numhypo;
} rg_stats_t;

int rupgen_get_num_rv(char* rup_geom_file, rg_stats_t *stats, int hypo_dist);
int rupgen_genslip_with_params(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf, int hypo_dist, float dt, char* params);
int rupgen_genslip(char* rup_geom_file, int slip, int hypo, rg_stats_t *stats, struct standrupformat* srf, int hypo_dist, float dt);
void free_srf_ptrs(struct standrupformat* srf);

#endif
