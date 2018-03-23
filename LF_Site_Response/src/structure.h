#ifndef STRUCTURE_H
#define STRUCTURE_H

#define FD_STATCHAR 8
#define STATCHAR 12
#define COMPCHAR 4
#define TITLCHAR 64

//for compatibility with params.h
#define mmv 60000
#define LV 1000
#define NP 100
#define NQ 600

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
   float qfexp;
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

struct rupture_variation {
  int rup_var_id;
  int slip_id;
  int hypo_id;
};

#endif
