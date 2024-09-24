#ifndef STRUCTURE_H
#define STRUCTURE_H

#define FD_STATCHAR 8
#define STATCHAR 12
#define COMPCHAR 4
#define TITLCHAR 64

#define mmv 30000

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

//All the information needed to regenerate a seismogram
//For version 12.10, try avoiding CyberShake-specific data
struct seisheader {
  char version[8];
  char site_name[8];
  //in case we think of something later
  char padding[8];
  int source_id;
  int rupture_id;
  int rup_var_id;
  float dt;
  int nt;
  int comps;
  float det_max_freq;
  float stoch_max_freq;
};

struct rotD_record {
  float period;
  float rotD100;
  int rotD100_angle;
  float rotD50;
};

struct vertical_rsp_record {
  float period;
  float rsp;
};

#define ARIAS_INTENSITY 0
#define ENERGY_INTEGRAL 1
#define CAV 2
#define DV 3
#define DA 4
#define D5_75 5
#define D5_95 6
#define D20_80 7

#define X_COMP 0
#define Y_COMP 1
#define Z_COMP 2

struct duration_record {
  int type;
  int type_value;
  int component;
  float value;
};

struct period_duration_record {
        int type;
        int type_value;
        int component;
        float period;
        float value;
};

#endif
