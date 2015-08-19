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

#endif
