#include "include.h"
#include "defs.h"
#include "structure.h"
#include "function.h"

void do_cnvlv(float *s,float *u,int nt,float *stf,int nstf)
{
float *sv, *sn, *se, *uv, *un, *ue;
int it, k, kend;

sv = s;
sn = s + nt;
se = s + 2*nt;

uv = u;
un = u + nt;
ue = u + 2*nt;

float acc[3];
for(it=0;it<nt;it++)
   {
   kend = it + 1;
   if(kend > nstf)
      kend = nstf;

   acc[0] = acc[1] = acc[2] = 0.0;
   for(k=0;k<kend;k++) {
	acc[0] += stf[k]*uv[it-k];
	acc[1] += stf[k]*un[it-k];
	acc[2] += stf[k]*ue[it-k];
   }
   sv[it] += acc[0];
   sn[it] += acc[1];
   se[it] += acc[2];
   }

   /*for(k=0;k<kend;k++)
      {
      sv[it] = sv[it] + stf[k]*uv[it-k];
      sn[it] = sn[it] + stf[k]*un[it-k];
      se[it] = se[it] + stf[k]*ue[it-k];
      }
   }*/
}

void sum_nostf(float *s,float *u,float *slip,int nt)
{
float *sv, *sn, *se, *uv, *un, *ue;
int it;

sv = s;
sn = s + nt;
se = s + 2*nt;

uv = u;
un = u + nt;
ue = u + 2*nt;

for(it=0;it<nt;it++)
   {
   sv[it] = sv[it] + (*slip)*uv[it];
   sn[it] = sn[it] + (*slip)*un[it];
   se[it] = se[it] + (*slip)*ue[it];
   }
}
