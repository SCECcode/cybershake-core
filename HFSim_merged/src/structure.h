#ifndef STRUCTURE_H
#define STRUCTURE_H

#define FD_STATCHAR 8
#define STATCHAR 12
#define COMPCHAR 4
#define TITLCHAR 64

//for compatibility with params.h
#define mmv 30000
#define LV 1000
#define NP 50
#define NQ 600

struct complex
   {
   float re;
   float im;
   };

struct statdata
   {
   char stat[STATCHAR];
   char comp[COMPCHAR];
   char stitle[TITLCHAR];
   int nt;
   float dt;
   int hr;
   int min;
   float sec;
   float edist;
   float az;
   float baz;
   };


struct stfpar
   {
   int nt;
   float dt;
   float trise;
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

struct slipfile {
   int nseg;
   float elon[LV];
   float elat[LV];
   int nx[LV];
   int ny[LV];
   float dx[LV];
   float dy[LV];
   float strike[LV];
   float dip[LV];
   float ravg[LV];
   float dtop[LV];
   float dhypo[LV];
   float shypo[LV];
//   float sp[NP][NQ][LV];
//   float tr[NP][NQ][LV];
//   float ti[NP][NQ][LV];
   float* sp;
   float* tr;
   float* ti;
};

#endif
