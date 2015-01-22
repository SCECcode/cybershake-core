#include "include.h"
#include "structure.h"
#include "function.h"

void get_srfpars(struct standrupformat *srf,int off, int ip,float *rt,float *vs,float *stk,float *dip,float *rak,struct mechparam *mpar)
{
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

mpar->nmech = 0;
mpar->flag[0] = 0;
mpar->flag[1] = 0;
mpar->flag[2] = 0;

if(apval_ptr[ip].nt1 > 0)
   {
   mpar->flag[mpar->nmech] = U1FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt2 > 0)
   {
   mpar->flag[mpar->nmech] = U2FLAG;
   mpar->nmech = mpar->nmech + 1;
   }
if(apval_ptr[ip].nt3 > 0)
   {
   mpar->flag[mpar->nmech] = U3FLAG;
   mpar->nmech = mpar->nmech + 1;
   }

*vs = sqrt(apval_ptr[ip].slip1*apval_ptr[ip].slip1
         + apval_ptr[ip].slip2*apval_ptr[ip].slip2
	 + apval_ptr[ip].slip3*apval_ptr[ip].slip3);
*stk = apval_ptr[ip].stk;
*dip = apval_ptr[ip].dip;
*rak = apval_ptr[ip].rake;
*rt = apval_ptr[ip].tinit;
}

void srf_stf(struct standrupformat *srf,int off,int ip,float *s,float *u,float *stf,int nt,float *dt,struct mechparam mp,float *space)
{
FILE *fpw;
int it, nstf, im;
float sum, *sptr, *uptr;

struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float fnt;
int resamp, ntpad, ntrsmp, gnt;

float tol = 1.0e-02;

float pratio_tol, mratio_tol;
float ratio_tol = 0.00001;

pratio_tol = 1.0 + ratio_tol;
mratio_tol = 1.0 - ratio_tol;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

zapit(stf,nt);

/*
if(apval_ptr[ip].nt1 == 0)
   return;
*/

/* for now, simply copy STF
   should add option to resample to dtout
   */

for(im=0;im<mp.nmech;im++)
   {
   if(mp.flag[im] == U1FLAG)
      {
      nstf = apval_ptr[ip].nt1;
      sptr = apval_ptr[ip].stf1;
      }
   else if(mp.flag[im] == U2FLAG)
      {
      nstf = apval_ptr[ip].nt2;
      sptr = apval_ptr[ip].stf2;
      }
   else if(mp.flag[im] == U3FLAG)
      {
      nstf = apval_ptr[ip].nt3;
      sptr = apval_ptr[ip].stf3;
      }

/* add factor of dt to prenormalize convolution */
   for(it=0;it<nstf;it++)
      stf[it] = (*dt)*sptr[it];

/* resample if needed */
   if((*dt)/apval_ptr[ip].dt > pratio_tol || (*dt)/apval_ptr[ip].dt < mratio_tol)
      {
      /*
      fprintf(stderr,"*** RESAMPLED diff= %13.5e ratio= %13.5e\n",(*dt)-apval_ptr[ip].dt,(*dt)/apval_ptr[ip].dt);
      */

      ntpad = 2*nstf;
      fnt = ntpad*apval_ptr[ip].dt/(*dt);
      gnt = (int)(fnt + 0.5);

      while(nt_tol(fnt,gnt) > tol)
         {
         ntpad++;
         fnt = ntpad*apval_ptr[ip].dt/(*dt);
         gnt = (int)(fnt + 0.5);
         }

      ntrsmp = (int)(fnt);
      if(ntrsmp > nt)
         {
         fprintf(stderr,"*** resampled nt > ntsum, exiting...\n");
         exit(-1);
         }

      if((*dt) < apval_ptr[ip].dt)
         resamp = 1;
      else
         resamp = -1;

      resample(stf,nstf,&apval_ptr[ip].dt,resamp,ntpad,ntrsmp,dt,space);

      nstf = ntrsmp;
      }

   uptr = u + 3*im*nt;
   do_cnvlv(s,uptr,nt,stf,nstf);
   }
}
