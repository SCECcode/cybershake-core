#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "misc.h"

void _init_plane_srf(struct standrupformat *srf,struct generic_slip *gslip,float *elon,float *elat,int nx,int ny,float *fl,float *fw,float *dx,float *dy,float *stk,float *dip,float *dtop,float *sh,float *dh)
{
struct srf_prectsegments *prseg_ptr;
float avglon, avglat, dmin, dold, se, sn;
int ig, iseg, nseg, nc, i, ip;
struct slippars *spar, *spar2;

float rperd = 0.017453293;

 dmin = 0.0;
sprintf(srf[0].type,"PLANE");

if(gslip->np > 0)
   {
   spar = gslip->spar;
   spar2 = (struct slippars *)_check_malloc(gslip->np*sizeof(struct slippars));

   nseg = 0;
   for(i=0;i<gslip->np;i++)
      {
      if(spar[i].segno > nseg)
         nseg = spar[i].segno;
      }
   nseg++;

   srf[0].srf_prect.nseg = nseg;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)_check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr = srf[0].srf_prect.prectseg;

   nc = 0;
   for(i=0;i<nseg;i++)
      {
      avglon = 0.0;
      avglat = 0.0;
      prseg_ptr[i].stk = 0.0;
      prseg_ptr[i].dip = 0.0;
      prseg_ptr[i].dlen = 0.0;
      prseg_ptr[i].dwid = 0.0;
      iseg = 0;
      dold = -1;
      prseg_ptr[i].nstk = 0;
      prseg_ptr[i].ndip = 0;
      for(ip=0;ip<gslip->np;ip++)
         {
         if(spar[ip].segno == i)
            {
            spar2[nc].lon = spar[ip].lon;
            spar2[nc].lat = spar[ip].lat;
            spar2[nc].dep = spar[ip].dep;
            spar2[nc].ds = spar[ip].ds;
            spar2[nc].dw = spar[ip].dw;
            spar2[nc].stk = spar[ip].stk;
            spar2[nc].dip = spar[ip].dip;
            spar2[nc].rake = spar[ip].rake;
            spar2[nc].slip = spar[ip].slip;
            spar2[nc].tinit = spar[ip].tinit;
            spar2[nc].segno = spar[ip].segno;

            if(iseg == 0)
               dmin = spar[ip].dep;

            if(spar[ip].dep == dmin)
               {
               avglon = avglon + spar[ip].lon;
               avglat = avglat + spar[ip].lat;
               prseg_ptr[i].nstk++;
               }

           if(spar[ip].dep != dold)
               {
               dold = spar[ip].dep;
               prseg_ptr[i].ndip++;
               }

            prseg_ptr[i].dlen = prseg_ptr[i].dlen + spar[ip].ds;
            prseg_ptr[i].dwid = prseg_ptr[i].dwid + spar[ip].dw;

            nc++;
            iseg++;

            if(iseg == 1)
               {
               prseg_ptr[i].stk = spar[ip].stk;
               prseg_ptr[i].dip = spar[ip].dip;
               }
            else
               {
               if(spar[ip].stk > prseg_ptr[i].stk/(float)(iseg-1.0) + 90.0)
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk - 180.0;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
                  }
               else if(spar[ip].stk < prseg_ptr[i].stk/(float)(iseg-1.0) - 90.0)
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk + 180.0;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
                  }
               else
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + spar[ip].dip;
                  }
               }
            }
         }
      prseg_ptr[i].stk = prseg_ptr[i].stk/(float)(iseg);
      prseg_ptr[i].dip = prseg_ptr[i].dip/(float)(iseg);

      if(prseg_ptr[i].dip > 90.0)
         {
         prseg_ptr[i].dip = 180.0 - prseg_ptr[i].dip;
         prseg_ptr[i].stk = prseg_ptr[i].stk - 180.0;
         }

      while(prseg_ptr[i].stk < 0.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk + 360.0;
      while(prseg_ptr[i].stk >= 360.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk - 360.0;

      prseg_ptr[i].dlen = prseg_ptr[i].dlen/(float)(iseg);
      prseg_ptr[i].dwid = prseg_ptr[i].dwid/(float)(iseg);

      prseg_ptr[i].flen = prseg_ptr[i].nstk*prseg_ptr[i].dlen;
      prseg_ptr[i].fwid = prseg_ptr[i].ndip*prseg_ptr[i].dwid;

      prseg_ptr[i].dtop = dmin - 0.5*prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip);
      if(prseg_ptr[i].dtop < 0.0)
         prseg_ptr[i].dtop = 0.0;

      avglon = avglon/(float)(prseg_ptr[i].nstk);
      avglat = avglat/(float)(prseg_ptr[i].nstk);
      se = -prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);

      _set_ll(&avglon,&avglat,&prseg_ptr[i].elon,&prseg_ptr[i].elat,&sn,&se);

      prseg_ptr[i].shyp = -999.9;
      prseg_ptr[i].dhyp = -999.9;

      fprintf(stderr,"%3d: %11.5f %11.5f %2d %2d %10.4f %10.4f %.1f %.1f %10.4f\n",i,prseg_ptr[i].elon,prseg_ptr[i].elat,prseg_ptr[i].nstk,prseg_ptr[i].ndip,prseg_ptr[i].flen,prseg_ptr[i].fwid,prseg_ptr[i].stk,prseg_ptr[i].dip,prseg_ptr[i].dtop);
      }

   free(spar2);
   }
else
   {
   srf[0].srf_prect.nseg = 1;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)_check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
   prseg_ptr = srf[0].srf_prect.prectseg;

   prseg_ptr[0].elon = *elon;
   prseg_ptr[0].elat = *elat;
   prseg_ptr[0].nstk = nx;
   prseg_ptr[0].ndip = ny;
   prseg_ptr[0].flen = *fl;
   prseg_ptr[0].fwid = *fw;
   prseg_ptr[0].dlen = *dx;
   prseg_ptr[0].dwid = *dy;
   prseg_ptr[0].stk = *stk;
   prseg_ptr[0].dip = *dip;
   prseg_ptr[0].dtop = *dtop;
   prseg_ptr[0].shyp = *sh;
   prseg_ptr[0].dhyp = *dh;
   }

srf[0].srf_apnts.np = 0;
for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
   srf[0].srf_apnts.np = srf[0].srf_apnts.np + (prseg_ptr[ig].nstk)*(prseg_ptr[ig].ndip);

srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)_check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));
}

void _write_srf(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ip, ig, ntot;

char pword[32];
int fdw;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = _croptrfile(file);

   _rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      _rite(fdw,srf->type,sizeof(srf->type));
      _rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      _rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   _rite(fdw,pword,sizeof(pword));
   _rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      _rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      _rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      _rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      _rite(fdw,&(apval_ptr[i].area),sizeof(float));
      _rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      _rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      _rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      _rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      _rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      _rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      _rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      _rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      _rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      _rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = _fopfile(file,"w");

   fprintf(fpw,"%s\n",srf->version);

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
        fprintf(fpw,"%10.4f %9.4f %5d %5d %8.2f %8.2f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %8.2f %8.2f %8.2f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   fprintf(fpw,"POINTS %d\n",apnts_ptr->np);
   for(i=0;i<apnts_ptr->np;i++)
      {
      fprintf(fpw,"%10.4f %9.4f %9.4f %4.0f %4.0f %12.5e %10.4f %12.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt);
      fprintf(fpw,"%4.0f %8.2f %6d %8.2f %6d %8.2f %6d\n",
                                              apval_ptr[i].rake,
                                              apval_ptr[i].slip1,
                                              apval_ptr[i].nt1,
                                              apval_ptr[i].slip2,
                                              apval_ptr[i].nt2,
                                              apval_ptr[i].slip3,
                                              apval_ptr[i].nt3);

      stf = apval_ptr[i].stf1;
      nt6 = (apval_ptr[i].nt1)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt1))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf2;
      nt6 = (apval_ptr[i].nt2)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt2))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf3;
      nt6 = (apval_ptr[i].nt3)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt3))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }
      }
   fclose(fpw);
   }
}

void _load_slip_srf(struct standrupformat *srf,struct stfpar *spar,struct pointsource *ps)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
float area;
float *stf;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot;

float dmin = 4.0;
float dmax = 6.0;
float rtfac, tzero;
float rtfac0, sabs;
float rtfac_min = 0.5;

dmin = spar->risetimedep - spar->risetimedep_range;
dmax = spar->risetimedep + spar->risetimedep_range;
rtfac0 = spar->risetimefac - 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
         ip = noff + i + j*prseg_ptr[iseg].nstk;
         ip0 = i + ioff + j*ntot;

         apval_ptr[ip].stf1 = (float *)_check_malloc(spar->nt*sizeof(float));
         stf = apval_ptr[ip].stf1;

         apval_ptr[ip].dt = spar->dt;

         sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
         if(sabs > MINSLIP)
            {
            if(ps[ip0].dep >= dmax)
               rtfac = 1.0;
            else if(ps[ip0].dep < dmax && ps[ip0].dep > dmin)
               rtfac = 1.0 + rtfac0*(dmax-(ps[ip0].dep))/(dmax-dmin);
            else
               rtfac = 1.0 + rtfac0;

            if(spar->rt_scalefac > 0.0)
               {
               rtfac = rtfac*exp(0.5*log(sabs))/spar->rt_scalefac;
               if(rtfac < rtfac_min)
                  rtfac = rtfac_min;
               }

            if(strcmp(spar->stype,"brune") == 0)
               {
               tzero = 0.1*exp(-1.0)*sqrt(ps[ip0].slip)/(1.2); /* assume slip in cm */    
               tzero = tzero*rtfac;

               apval_ptr[ip].nt1 = _gen_brune_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"urs") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = _gen_2tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt,&ps[ip0].dep);
               }
            else if(strcmp(spar->stype,"ucsb") == 0)
               {
               tzero = rtfac*spar->trise;
               apval_ptr[ip].nt1 = _gen_ucsb_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"tri") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = _gen_tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            }
         else
            apval_ptr[ip].nt1 = 0;

         if(apval_ptr[ip].nt1)
            apval_ptr[ip].stf1 = (float *)_check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
         else
            {
            free(apval_ptr[ip].stf1);
            apval_ptr[ip].stf1 = NULL;
            }

         apval_ptr[ip].lon = ps[ip0].lon;
         apval_ptr[ip].lat = ps[ip0].lat;
         apval_ptr[ip].dep = ps[ip0].dep;
         apval_ptr[ip].stk = ps[ip0].stk;
         apval_ptr[ip].dip = ps[ip0].dip;
         apval_ptr[ip].area = ps[ip0].area;
         apval_ptr[ip].rake = ps[ip0].rak;
         apval_ptr[ip].slip1 = ps[ip0].slip;

         apval_ptr[ip].slip2 = 0.0;
         apval_ptr[ip].nt2 = 0;
         apval_ptr[ip].stf2 = NULL;
         apval_ptr[ip].slip3 = 0.0;
         apval_ptr[ip].nt3 = 0;
         apval_ptr[ip].stf3 = NULL;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void _load_rupt_srf(struct standrupformat *srf,struct pointsource *ps,float *sh,float *dh)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot;

(srf->srf_prect).prectseg[0].shyp = *sh;
(srf->srf_prect).prectseg[0].dhyp = *dh;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
         ip = noff + i + j*prseg_ptr[iseg].nstk;
         ip0 = i + ioff + j*ntot;

         apval_ptr[ip].tinit = ps[ip0].rupt;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

int _gen_2tri_stf(float *slip,float *trise,float *stf,int nt,float *dt,float *z0)
{
int it, nstf;
int ip, it0, it1, it2;
float tr, amp, a0;
float sum;
float alpha = 0.1;      /* 1st triangle has pulse width = 2*alpha*trise */
float betadeep = 0.2;       /* 2nd triangle has amplitude = beta*A (z0>dmax)*/
float betashal = 0.5;       /* 2nd triangle has amplitude = beta*A (z0<dmin)*/
float beta, dbdd;

float dmin = 4.0;
float dmax = 6.0;

dbdd = (betadeep - betashal)/(dmax-dmin);

if((*z0) >= dmax)
   beta = betadeep;
else if((*z0) < dmax && (*z0) > dmin)
   beta = betadeep - (dmax-(*z0))*dbdd;
else
   beta = betashal;

_zapit(stf,nt);

tr = (*trise);
alpha = alpha*tr;

it0 = (int)((alpha)/(*dt) + 0.5);
if(it0 < 2)
   it0 = 2;
it1 = (int)((tr)/(*dt) + 0.5);
if(it1 < 4)
   it1 = 4;

it2 = (2 - beta)*it0;

a0 = 1.0;
amp = a0/(float)(it0);

for(it=0;it<it0;it++)
   stf[it] = it*amp;

for(it=it0;it<it2;it++)
   stf[it] = (2*it0-it)*amp;

amp = beta*a0/(float)(it1-it2);

for(it=it2;it<it1;it++)
   stf[it] = beta*a0 + (it2-it)*amp;

nstf = nt-1;
while(stf[nstf] == (float)(0.0) && nstf)
   nstf--;

if(nstf == 0)
   return(0);

if(nstf < nt-1)
   nstf = nstf + 2;;

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int _gen_tri_stf(float *slip,float *t0,float *stf,int nt,float *dt)
{
int it, nstf;
float tau, dadt, amp, tau2;
float sum, t, alpha;
float pi = 3.141592654;

_zapit(stf,nt);

tau = (*t0);

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

tau2 = 0.5*tau;
amp = (*slip)/tau2;
dadt = amp/tau2;
for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t <= tau2)
      alpha = t*dadt;
   else if(t <= tau)
      alpha = 2.0*amp - t*dadt;

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int _gen_ucsb_stf(float *slip,float *t0,float *stf,int nt,float *dt)
{
int it, nstf;
float tau, tau1, tau2, tau1x2, arg1, arg2;
float sum, t, alpha;
float pi = 3.141592654;

_zapit(stf,nt);

tau = (*t0);
tau1 = 0.13*tau;
tau2 = tau - tau1;
tau1x2 = 2.0*tau1;

nstf = (int)((tau)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

for(it=0;it<nstf;it++)
   {
   t = it*(*dt);

   alpha = 0.0;
   if(t < tau1)
      {
      arg1 = pi*t/tau1;
      arg2 = 0.5*arg1;
      alpha = 0.7 - 0.7*cos(arg1) + 0.6*sin(arg2);
      }
   else if(t < tau1x2)
      {
      arg1 = pi*t/tau1;
      arg2 = pi*(t - tau1)/tau2;
      alpha = 1.0 - 0.7*cos(arg1) + 0.3*cos(arg2);
      }
   else if(t < tau)
      {
      arg1 = pi*(t - tau1)/tau2;
      alpha = 0.3 + 0.3*cos(arg1);
      }

   stf[it] = alpha;
   }

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

/* scale STF by slip */
sum = (*slip)/sum;
for(it=0;it<nstf;it++)
   stf[it] = stf[it]*sum;

return(nstf);
}

int _gen_brune_stf(float *slip,float *t0,float *stf,int nt,float *dt)
{
int it, nstf;
float t95, tend, sfac, tfac;
float sum;

_zapit(stf,nt);

t95 = 1.745*exp(1.0)*(*t0);
tend = 3.0*t95;

nstf = (int)((tend)/(*dt) + 0.5);
if(nstf > nt)
   nstf = nt;

if(nstf == 0)
   return(0);

sfac = (*slip)/(*t0);
tfac = (*dt)/(*t0);
for(it=0;it<nstf;it++)
   stf[it] = sfac*(it*tfac)*exp(-it*tfac);

sum = 0.0;
for(it=0;it<nstf;it++)
   sum = sum + (*dt)*stf[it];

if(sum <= 0.0)
   return(0);

return(nstf);
}

