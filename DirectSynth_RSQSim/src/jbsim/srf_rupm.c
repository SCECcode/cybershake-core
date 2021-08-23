#include "include.h"
#include "defs.h"
#include "structure.h"
#include "srf_structure.h"
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
	 //printf("fnt=%f, gnt=%d, tol=%f.\n", fnt, gnt, tol);
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

      //fprintf(stderr, "Resampling with nstf=%d, ntpad=%d, ntrsmp=%d, dt=%f\n", nstf, ntpad, ntrsmp, *dt);
      //Moved resample checks here, since resample() doesn't have access to nt
      if (nt<2*ntpad) {
	fprintf(stderr, "Error in resample: nt=%d, but ntpad=%d, so resample would overflow stf array.  Aborting.", nt, ntpad);
	exit(2);
      }
      if (nt<2*ntrsmp) {
	fprintf(stderr, "Error in resample: nt=%d, but ntrsmp=%d, so resample would overflow stf array.  Aborting.", nt, ntrsmp);
	exit(2);
      }

      resample(stf,nstf,&apval_ptr[ip].dt,resamp,ntpad,ntrsmp,dt,space);

      nstf = ntrsmp;
      }

   uptr = u + 3*im*nt;
   do_cnvlv(s,uptr,nt,stf,nstf);
   }
}

void get_srfpars_v2(struct standrupformat *srf,int off, int ip,float *rt,float *vs,struct mechparam *mpar)
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

mpar->stk = apval_ptr[ip].stk;
mpar->dip = apval_ptr[ip].dip;
mpar->rak = apval_ptr[ip].rake;

*rt = apval_ptr[ip].tinit;
}


void srf_stf_v2(struct standrupformat *srf,int off,int ip,float *s,float *u,float *stf,int nt,float *dt,struct mechparam mp,float *unused)
{
FILE *fpw;
int it, nstf, im;
float sum, *sptr, *uptr;
float *space, *sptr2;
int i;

struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float fnt;
int resamp, ntpad, ntrsmp, gnt;

float slip, sfac, dt_stf, dt_tmp, ds_dt;
int nt_fac, nt_tmp, jj, kk;
int realloc_flag = 0;

float tol = 1.0e-02;

float pratio_tol, mratio_tol;
float ratio_tol = 0.001;

pratio_tol = 1.0 + ratio_tol;
mratio_tol = 1.0 - ratio_tol;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals + off;

zapit(stf,nt);

space = NULL;
sptr2 = NULL;

for(im=0;im<mp.nmech;im++)
   {
   if(mp.flag[im] == U1FLAG)
      {
      nstf = apval_ptr[ip].nt1;
	  //Make a copy of stf1 contents to avoid freeing here
      //sptr = apval_ptr[ip].stf1;
      sptr = check_malloc(nstf*sizeof(float));
	  memcpy(sptr, apval_ptr[ip].stf1, sizeof(float)*nstf);
      slip = apval_ptr[ip].slip1;
      }
   else if(mp.flag[im] == U2FLAG)
      {
      nstf = apval_ptr[ip].nt2;
      //sptr = apval_ptr[ip].stf2;
      sptr = check_malloc(nstf*sizeof(float));
      memcpy(sptr, apval_ptr[ip].stf2, sizeof(float)*nstf);
      slip = apval_ptr[ip].slip2;
      }
   else if(mp.flag[im] == U3FLAG)
      {
      nstf = apval_ptr[ip].nt3;
      //sptr = apval_ptr[ip].stf3;
      sptr = check_malloc(nstf*sizeof(float));
      memcpy(sptr, apval_ptr[ip].stf3, sizeof(float)*nstf);
      slip = apval_ptr[ip].slip3;
      }

   dt_stf = apval_ptr[ip].dt;

/* resample if needed */

   if((*dt)/dt_stf > pratio_tol || (*dt)/dt_stf < mratio_tol)
      {
      if(nstf == 1)   /* add zero points in front and back for this case */
         {
         nstf = nstf + 2;
         sptr = (float *)check_realloc((void *)sptr,nstf*sizeof(float));
         realloc_flag = 1;

         sptr[2] = 0.0;
         sptr[1] = sptr[0];
         sptr[0] = 0.0;
         }
      else   /* add zero point in back for this case */
         {
         nstf++;
         sptr = (float *)check_realloc((void *)sptr,nstf*sizeof(float));
         realloc_flag = 1;

         sptr[nstf-1] = 0.0;
         }

/* first try linear interpolation as this best mimics target slip-rate */

      nt_fac = (int)(dt_stf/(*dt) + 0.5);
      if(nt_fac >= 2)  /* otherwise, skip linear interpolation */
         {
         nt_tmp = nt_fac*(nstf-1) + 1;
         dt_tmp = dt_stf/nt_fac;
         sptr2 = (float *)check_realloc((void *)sptr2,nt_tmp*sizeof(float));
         for(it=0;it<nt_tmp;it++)
            sptr2[it] = 0.0;

         sptr2[0] = sptr[0];
         for(it=1;it<nstf;it++)
            {
            ds_dt = (sptr[it]-sptr[it-1])/dt_stf;

            for(jj=1;jj<=nt_fac;jj++)
               {
               kk = jj + (it-1)*nt_fac;
               sptr2[kk] = sptr[it-1] + jj*dt_tmp*ds_dt;
               }
            }

         sptr = (float *)check_realloc((void *)sptr,nt_tmp*sizeof(float));
         realloc_flag = 1;
         for(it=0;it<nt_tmp;it++)
            sptr[it] = sptr2[it];

         nstf = nt_tmp;
         dt_stf = dt_tmp;
         }

/* second, check to see if additional adjustment is still needed */

      if((*dt)/dt_stf > pratio_tol || (*dt)/dt_stf < mratio_tol)
         {
         ntpad = 2*nstf;
         fnt = ntpad*dt_stf/(*dt);
         gnt = (int)(fnt + 0.5);
         while(nt_tol(fnt,gnt) > tol)
            {
            ntpad++;
            fnt = ntpad*dt_stf/(*dt);
            gnt = (int)(fnt + 0.5);
            }

         ntrsmp = (int)(fnt);

         if((*dt) < dt_stf)
        {
        space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
        sptr2 = (float *)check_realloc((void *)sptr2,2*ntrsmp*sizeof(float));
            for(it=0;it<2*ntrsmp;it++)
               sptr2[it] = 0.0;

            resamp = 1;
        }
         else
        {
        space = (float *) check_realloc ((void *)space,2*ntpad*sizeof(float));
        sptr2 = (float *)check_realloc((void *)sptr2,2*ntpad*sizeof(float));
            for(it=0;it<2*ntpad;it++)
               sptr2[it] = 0.0;

            resamp = -1;
        }

         for(it=0;it<nstf;it++)
            sptr2[it] = sptr[it];

         resample(sptr2,nstf,&dt_stf,resamp,ntpad,ntrsmp,dt,space);

         sptr = sptr2;
         nstf = ntrsmp;
         }

/* double check to ensure slip-rate integrates to final slip */

      if(slip > 0.0)
         {
         sfac = 0.0;
         for(it=0;it<nstf;it++)
            sfac = sfac + (*dt)*sptr[it];
         sfac = slip/sfac;
         for(it=0;it<nstf;it++)
            sptr[it] = sfac*sptr[it];
         }
      }

   if(nstf > nt)
      nstf = nt;

/* add factor of dt to prenormalize convolution */
   for(it=0;it<nstf;it++)
      stf[it] = (*dt)*sptr[it];

   uptr = u + 3*im*nt;
   do_cnvlv(s,uptr,nt,stf,nstf);
   }

if(realloc_flag == 1) {
   free(sptr);
   sptr = NULL;
}

free(space);
free(sptr2);
}

