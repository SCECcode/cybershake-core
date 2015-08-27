#include "include.h"
#include "structure.h"
#include "function.h"

double frand(void);
double sfrand(long *);

int wcc_add(float** lf_seis, float** hf_seis, float** merged_seis, int ncomp, float dt, int nt) {
struct statdata shead1, shead2, shead3;
float *s1, *s2, *p;

float f1 = 1.0;
float f2 = 1.0;

float t1 = 0.0;
float t2 = 0.0;

char infile1[128];
char infile2[128];
char outfile[128];

int inbin1 = 0;
int inbin2 = 0;
int outbin = 0;

int it;
float add_rand = 0.0;

/*
setpar(ac,av);
getpar("f1","f",&f1);
getpar("t1","f",&t1);
mstpar("infile1","s",infile1);
getpar("f2","f",&f2);
getpar("t2","f",&t2);
mstpar("infile2","s",infile2);
mstpar("outfile","s",outfile);
getpar("inbin1","d",&inbin1);
getpar("inbin2","d",&inbin2);
getpar("outbin","d",&outbin);
getpar("add_rand","f",&add_rand);
endpar();*/

shead1.nt = nt;
shead1.dt = dt;
shead2.nt = nt;
shead2.dt = dt;
shead2.hr = shead1.hr = 0;
shead2.min = shead1.min = 0;
shead2.sec = shead1.sec = 0.0;
shead2.edist = shead1.edist = 0.0;
shead2.az = shead1.az = 0.0;
shead2.baz = shead1.baz = 0.0;


s1 = NULL;
//s1 = read_wccseis(infile1,&shead1,s1,inbin1);
s2 = NULL;
//s2 = read_wccseis(infile2,&shead2,s2,inbin2);

int i;
for (i=0; i<ncomp; i++) {

s1 = lf_seis[i];
s2 = hf_seis[i];
p = merged_seis[i];

//p = (float *) check_malloc ((shead1.nt+shead2.nt)*size_float);

sum(s1,&shead1,&f1,&t1,s2,&shead2,&f2,&t2,p,&shead3);

//strcpy(shead3.stat,shead1.stat);
//strcpy(shead3.comp,shead1.comp);
//sprintf(shead3.stitle,"summed output");

if(add_rand > (float)(0.0))
   {
   for(it=0;it<shead3.nt;it++)
      p[it] = p[it] + add_rand*frand();
   }

} //end for
//write_wccseis(outfile,&shead3,p,outbin);
}

int sum(float* s1,struct statdata* h1,float* f1,float* t1,float* s2,struct statdata* h2,float* f2,float* t2,float* p,struct statdata* hd) {
float time1, time2;
int it, n;

hd->dt = h1->dt;
hd->edist = h1->edist;
hd->az = h1->az;
hd->baz = h1->baz;

/* set output start time to be earliest start time of two input files */

time1 = h1->hr*3600.0 + h1->min*60.0 + h1->sec + *t1;
time2 = h2->hr*3600.0 + h2->min*60.0 + h2->sec + *t2;

n = 0;
if(time1 < time2)
   {
   hd->hr = h1->hr;
   hd->min = h1->min;
   hd->sec = h1->sec + *t1;

   n = ((time2-time1)/hd->dt + 0.5);

   for(it=0;it<n;it++)
      p[it] = (*f1)*s1[it];

   h1->nt = h1->nt - n;
   s1 = s1 + n;
   p = p + n;
   }
else
   {
   hd->hr = h2->hr;
   hd->min = h2->min;
   hd->sec = h2->sec + *t2;

   n = ((time1-time2)/hd->dt + 0.5);

   for(it=0;it<n;it++)
      p[it] = (*f2)*s2[it];

   h2->nt = h2->nt - n;;
   s2 = s2 + n;
   p = p + n;
   }

if(h1->nt > h2->nt)
   {
   hd->nt = h1->nt + n;

   for(it=0;it<h2->nt;it++)
      p[it] = (*f1)*s1[it] + (*f2)*s2[it];

   for(it=h2->nt;it<h1->nt;it++)
      p[it] = (*f1)*s1[it];
   }
else
   {
   hd->nt = h2->nt + n;

   for(it=0;it<h1->nt;it++) {
      p[it] = (*f1)*s1[it] + (*f2)*s2[it];
   }

   for(it=h1->nt;it<h2->nt;it++)
      p[it] = (*f2)*s2[it];
   }
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
/*
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}*/

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
/*double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}*/
