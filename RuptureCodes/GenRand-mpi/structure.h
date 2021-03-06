#ifndef RUPGEN_STRUCT_H
#define RUPGEN_STRUCT_H

#include "defs.h"


struct complex
   {
   float re;
   float im;
   };

struct velmodel
   {
   int nlay;
   float *vp;
   double *vs;    /* need double for ray tracing to get ruptime */
   float *den;
   float *th;
   float *dep;
   float *mu;     /* in CMS units */
   double *invb2; /* need double for ray tracing to get ruptime */
   };

struct stfpar
   {
   int nt;
   float dt;
   float trise;
   float risetimedep;
   float risetimedep_range;
   float risetimefac;
   float rt_scalefac;
   char stype[32];
   };

struct pointsource
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float rak;
   float area;
   float slip;
   float rupt;
   };

struct srf_apointvalues
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float area;
   float tinit;
   float dt;
   float rake;
   float slip1;
   int nt1;
   float slip2;
   int nt2;
   float slip3;
   int nt3;
   float *stf1;
   float *stf2;
   float *stf3;
   };

struct srf_allpoints
   {
   int np;
   struct srf_apointvalues *apntvals;
   };

struct srf_prectsegments
   {
   float elon;
   float elat;
   int nstk;
   int ndip;
   float flen;
   float fwid;
   float dlen;
   float dwid;
   float stk;
   float dip;
   float dtop;
   float shyp;
   float dhyp;
   };

struct srf_planerectangle
   {
   int nseg;
   struct srf_prectsegments *prectseg;
   };

struct standrupformat
   {
   char version[32];
   char type[32];
   struct srf_planerectangle srf_prect;
   struct srf_allpoints srf_apnts;
   };

struct slippars
   {
   float lon;
   float lat;
   float dep;
   float ds;
   float dw;
   float stk;
   float dip;
   float rake;
   float slip;
   float tinit;
   int segno;
   };

struct generic_slip
   {
   int np;
   struct slippars *spar;
   };


/* Rupture generator config */
typedef struct rg_conf_t {
  char rupdir[MAX_FILENAME];
  char logdir[MAX_FILENAME];
  char erfid[MAX_ID_LEN];
  char rvid[MAX_ID_LEN];
  char gridftp[MAX_FILENAME];
  char pool[MAX_FILENAME];
} rg_conf_t;


/* Rupture generation statistics */
typedef struct rg_stats_t {
  int numslip;
  int numhypo;
} rg_stats_t;


/* Rupture file info */
typedef struct rg_rfile_t {
  char filename[MAX_FILENAME];
  int index;
  float mag;
  int src;
  int rup;
  rg_stats_t stats;
  int request;
} rg_rfile_t;

#endif

