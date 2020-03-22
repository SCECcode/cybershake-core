#include "include.h"

add_src(stfp,nt,dt,pvptr,medptr,isrc,nx,nz,it,srcs,iflag)
struct pntsrcs *srcs;
register float *stfp, *dt, **pvptr, **medptr;
int isrc, nx, it, nz, iflag, nt;
{
struct momenttensor *mt;
struct adjointsource *adj;
float *vx, *vy, *vz, *medx, *medy, *medz;
float *vx0, *vy0, *vz0, *vx1, *vy1, *vz1, *vx2, *vy2, *vz2;
float *txx, *tyy, *tzz, *txy0, *txy1, *txz, *tyz0, *tyz1;
float f;
int itx, itp, itm, ip, nxm, nxp;
float *st;

float qrt = 0.25;

if(!it)
   return;

if(srcs->ffault != 2 && srcs->adjoint != 1)
   st = stfp + nt*srcs->stfindx[isrc];

ip = srcs->iz[isrc]*nx + srcs->ix[isrc];
nxm = nx - 1;
nxp = nx + 1;

itm = it - 1;
itp = it + 1;

if(iflag == FSRC)
   {
   vx  = pvptr[2];
   vy  = pvptr[2]  +   nx*nz;
   vz  = pvptr[2]  + 2*nx*nz;
   medx = medptr[2] + 5*nx*nz;  /* med = 1/rho */
   medy = medptr[2] + 6*nx*nz;  /* med = 1/rho */
   medz = medptr[2] + 7*nx*nz;  /* med = 1/rho */

   if(srcs->bforce)
      {
      vx[ip] = vx[ip] + srcs->fxsrc[isrc]*st[it]*medx[ip];
      vy[ip] = vy[ip] + srcs->fysrc[isrc]*st[it]*medy[ip];
      vz[ip] = vz[ip] + srcs->fzsrc[isrc]*st[it]*medz[ip];
      }

   if(srcs->expl)
      {
      vx[ip] = vx[ip] + srcs->expl*(*dt)*st[it]*medx[ip];
      vy[ip] = vy[ip] + srcs->expl*(*dt)*st[it]*medy[ip];
      vz[ip] = vz[ip] + srcs->expl*(*dt)*st[it]*medz[ip];

      vx[ip-1] = vx[ip-1] - srcs->expl*(*dt)*st[it]*medx[ip-1];
      vz[ip-nx] = vz[ip-nx] - srcs->expl*(*dt)*st[it]*medz[ip-nx];

      vy  = pvptr[1]  +   nx*nz;
      medy = medptr[1] + 6*nx*nz;
      vy[ip] = vy[ip] - srcs->expl*(*dt)*st[it]*medy[ip];
      }

   if(srcs->adjoint == 1 && it - srcs->it[isrc] >= 0)
      {
      itx = it - srcs->it[isrc];
      adj = &(srcs->adjsrc[isrc]);

      if(itx < adj->nt)
         {
         vx[ip]  = vx[ip] + adj->sx[itx]*medx[ip];
         vy[ip]  = vy[ip] + adj->sy[itx]*medy[ip];
         vz[ip]  = vz[ip] + adj->sz[itx]*medz[ip];
	 }
      }
   }
else if(iflag == PSRC)
   {
   txx = pvptr[2] + 3*nx*nz;
   tyy = pvptr[2] + 4*nx*nz;
   tzz = pvptr[2] + 5*nx*nz;

   txy1 = pvptr[2] + 6*nx*nz;
   txz  = pvptr[2] + 7*nx*nz;
   tyz1 = pvptr[2] + 8*nx*nz;

   txy0 = pvptr[1] + 6*nx*nz;
   tyz0 = pvptr[1] + 8*nx*nz;

   if(srcs->psrc)
      {
      txx[ip] = txx[ip] + srcs->psrc*(*dt)*st[it];
      tyy[ip] = tyy[ip] + srcs->psrc*(*dt)*st[it];
      tzz[ip] = tzz[ip] + srcs->psrc*(*dt)*st[it];
      }
   if(srcs->eqsrc)
      {
      txx[ip] = txx[ip] + srcs->eqsrc*(*dt)*st[it];
      tyy[ip] = tyy[ip] - srcs->eqsrc*(*dt)*st[it];
      }
   if((srcs->dblcpl || srcs->pointmt) && it - srcs->doublec.itsrc[isrc] > 0)
      {
      itx = it - srcs->doublec.itsrc[isrc];
      mt = &(srcs->momten[isrc]);

      txx[ip]  = txx[ip] - mt->mxx[0]*st[itx];
      tyy[ip]  = tyy[ip] - mt->myy[0]*st[itx];
      tzz[ip]  = tzz[ip] - mt->mzz[0]*st[itx];

/* without averaging

*/

      txy1[ip]    = txy1[ip]    - mt->mxy[0]*st[itx];
      txz[ip]     = txz[ip]     - mt->mxz[0]*st[itx];
      tyz1[ip]    = tyz1[ip]    - mt->myz[0]*st[itx];

/* with averaging

This may cause problems in parallel code if source is right at a boundary between two nodes!!!

      txy0[ip-1]  = txy0[ip-1]  - qrt*mt->mxy[0]*st[itx];
      txy0[ip]    = txy0[ip]    - qrt*mt->mxy[0]*st[itx];
      txy1[ip-1]  = txy1[ip-1 ] - qrt*mt->mxy[0]*st[itx];
      txy1[ip]    = txy1[ip]    - qrt*mt->mxy[0]*st[itx];

      txz[ip-nxp] = txz[ip-nxp] - qrt*mt->mxz[0]*st[itx]; 
      txz[ip-nx]  = txz[ip-nx]  - qrt*mt->mxz[0]*st[itx];
      txz[ip-1]   = txz[ip-1]   - qrt*mt->mxz[0]*st[itx];
      txz[ip]     = txz[ip]     - qrt*mt->mxz[0]*st[itx];

      tyz0[ip-nx] = tyz0[ip-nx] - qrt*mt->myz[0]*st[itx]; 
      tyz0[ip]    = tyz0[ip]    - qrt*mt->myz[0]*st[itx];
      tyz1[ip-nx] = tyz1[ip-nx] - qrt*mt->myz[0]*st[itx];
      tyz1[ip]    = tyz1[ip]    - qrt*mt->myz[0]*st[itx];

*/
      }
   if(srcs->ffault == 2 && it - srcs->it[isrc] >= 0)
      {
      itx = it - srcs->it[isrc];
      mt = &(srcs->momten[isrc]);

      if(itx < mt->nt)
         {
         txx[ip]  = txx[ip] - mt->mxx[itx];
         tyy[ip]  = tyy[ip] - mt->myy[itx];
         tzz[ip]  = tzz[ip] - mt->mzz[itx];

/* without averaging */

         txy1[ip]    = txy1[ip]    - mt->mxy[itx];
         txz[ip]     = txz[ip]     - mt->mxz[itx];
         tyz1[ip]    = tyz1[ip]    - mt->myz[itx];

/* with averaging

This may cause problems in parallel code if source is right at a boundary between two nodes!!!

         txy0[ip-1]  = txy0[ip-1]  - qrt*mt->mxy[itx];
         txy0[ip]    = txy0[ip]    - qrt*mt->mxy[itx];
         txy1[ip-1]  = txy1[ip-1 ] - qrt*mt->mxy[itx];
         txy1[ip]    = txy1[ip]    - qrt*mt->mxy[itx];

         txz[ip-nxp] = txz[ip-nxp] - qrt*mt->mxz[itx]; 
         txz[ip-nx]  = txz[ip-nx]  - qrt*mt->mxz[itx];
         txz[ip-1]   = txz[ip-1]   - qrt*mt->mxz[itx];
         txz[ip]     = txz[ip]     - qrt*mt->mxz[itx];

         tyz0[ip-nx] = tyz0[ip-nx] - qrt*mt->myz[itx]; 
         tyz0[ip]    = tyz0[ip]    - qrt*mt->myz[itx];
         tyz1[ip-nx] = tyz1[ip-nx] - qrt*mt->myz[itx];
         tyz1[ip]    = tyz1[ip]    - qrt*mt->myz[itx];
*/
	 }
      }
   }
}

void add_srcP3(float *stfp,int nt,float *dt,float *pvf,float *medf,int isrc,int nx,int nz,int it,struct pntsrcs *srcs,int iflag,struct nodeinfo *ni)
{
struct momenttensor *mt;
struct adjointsource *adj;
float *vx, *vy, *vz, *medx, *medy, *medz;
float *vx0, *vy0, *vz0, *vx1, *vy1, *vz1, *vx2, *vy2, *vz2;
float *txx, *tyy, *tzz, *txy0, *txy1, *txz, *tyz0, *tyz1;
float f;
int itx, ipminusy, ip, ipmed, nxp;
float *st;

if(!it)
   return;

if(srcs->ffault != 2 && srcs->adjoint != 1)
   st = stfp + nt*srcs->stfindx[isrc];

ip = srcs->iy[isrc]*N_WAVE_VARS*nx*nz + srcs->iz[isrc]*nx + srcs->ix[isrc];
ipmed = srcs->iy[isrc]*N_MED_VARS*nx*nz + srcs->iz[isrc]*nx + srcs->ix[isrc];
ipminusy = (srcs->iy[isrc]-1)*N_WAVE_VARS*nx*nz + srcs->iz[isrc]*nx + srcs->ix[isrc];

nxp = nx + 1;

if(iflag == FSRC)
   {
   vx  = pvf + ip;
   vy  = pvf + ip  +   nx*nz;
   vz  = pvf + ip  + 2*nx*nz;
   medx = medf + ipmed + 5*nx*nz;  /* med = 1/rho */
   medy = medf + ipmed + 6*nx*nz;  /* med = 1/rho */
   medz = medf + ipmed + 7*nx*nz;  /* med = 1/rho */

   if(srcs->bforce)
      {
      vx[0] = vx[0] + srcs->fxsrc[isrc]*st[it]*medx[0];
      vy[0] = vy[0] + srcs->fysrc[isrc]*st[it]*medy[0];
      vz[0] = vz[0] + srcs->fzsrc[isrc]*st[it]*medz[0];
      }

   if(srcs->adjoint == 1 && it - srcs->it[isrc] >= 0)
      {
      itx = it - srcs->it[isrc];
      adj = &(srcs->adjsrc[isrc]);

      if(itx < adj->nt)
         {
         vx[0]  = vx[0] + adj->sx[itx]*medx[0];
         vy[0]  = vy[0] + adj->sy[itx]*medy[0];
         vz[0]  = vz[0] + adj->sz[itx]*medz[0];
	 }
      }
   }
else if(iflag == PSRC)
   {
   if((srcs->ffault == 2 || srcs->dblcpl || srcs->pointmt) && it - srcs->it[isrc] >= 0)
      {
      txx = pvf + ip + 3*nx*nz;
      tyy = pvf + ip + 4*nx*nz;
      tzz = pvf + ip + 5*nx*nz;

      txy1 = pvf + ip + 6*nx*nz;
      txz  = pvf + ip + 7*nx*nz;
      tyz1 = pvf + ip + 8*nx*nz;

      txy0 = pvf + ipminusy + 6*nx*nz;
      tyz0 = pvf + ipminusy + 8*nx*nz;

      itx = it - srcs->it[isrc];
      mt = &(srcs->momten[isrc]);

      if(itx < mt->nt)
         {
         txx[0]    = txx[0]    - srcs->MTavg[isrc].wxx_x0y0z0*mt->mxx[itx];
         tyy[0]    = tyy[0]    - srcs->MTavg[isrc].wyy_x0y0z0*mt->myy[itx];
         tzz[0]    = tzz[0]    - srcs->MTavg[isrc].wzz_x0y0z0*mt->mzz[itx];

         txy1[0]   = txy1[0]   - srcs->MTavg[isrc].wxy_x0y0z0*mt->mxy[itx];
         txy1[-1]  = txy1[-1]  - srcs->MTavg[isrc].wxy_xmy0z0*mt->mxy[itx];
         txy0[0]   = txy0[0]   - srcs->MTavg[isrc].wxy_x0ymz0*mt->mxy[itx];
         txy0[-1]  = txy0[-1]  - srcs->MTavg[isrc].wxy_xmymz0*mt->mxy[itx];

         txz[0]    = txz[0]    - srcs->MTavg[isrc].wxz_x0y0z0*mt->mxz[itx];
         txz[-1]   = txz[-1]   - srcs->MTavg[isrc].wxz_xmy0z0*mt->mxz[itx];
         txz[-nx]  = txz[-nx]  - srcs->MTavg[isrc].wxz_x0y0zm*mt->mxz[itx];
         txz[-nxp] = txz[-nxp] - srcs->MTavg[isrc].wxz_xmy0zm*mt->mxz[itx]; 

         tyz1[0]   = tyz1[0]   - srcs->MTavg[isrc].wyz_x0y0z0*mt->myz[itx];
         tyz0[0]   = tyz0[0]   - srcs->MTavg[isrc].wyz_x0ymz0*mt->myz[itx];
         tyz1[-nx] = tyz1[-nx] - srcs->MTavg[isrc].wyz_x0y0zm*mt->myz[itx];
         tyz0[-nx] = tyz0[-nx] - srcs->MTavg[isrc].wyz_x0ymzm*mt->myz[itx]; 
	 }
      }
   }
}

add_plane(st,pvptr,nx,nz,it,psrc)
struct planesrc *psrc;
float *st, **pvptr;
int nx, it, nz;
{
float *vp, *v1, *v2;
int offp, off1, off2, ix, ip, i;

if(!it)
   return;

if(psrc->p > 0.0)    /* p waves on vz */
   {
   off1 = 0;
   off2 = nx*nz;
   offp = 2*nx*nz;
   }
else if(psrc->sh > 0.0)    /* sh waves on vy */
   {
   off1 = 0;
   offp = nx*nz;
   off2 = 2*nx*nz;
   }
else if(psrc->sv > 0.0)    /* sv waves on vx */
   {
   offp = 0;
   off1 = nx*nz;
   off2 = 2*nx*nz;
   }

for(i=0;i<4;i++)
   {
   vp = pvptr[i] + offp;
   v1 = pvptr[i] + off1;
   v2 = pvptr[i] + off2;

   for(ix=0;ix<nx;ix++)
      {
      ip = psrc->intercept*nx + ix;
      vp[ip] = st[it];
      v1[ip] = 0.0;
      v2[ip] = 0.0;
      }
   }
}

getsource(st,dt,nt,tz,nb,intsrc,type,name,bfilt,flo,fhi,indx,id,tsh,stfdir,sfile)
float *st, *dt, *tz, *flo, *fhi, *tsh;
int nt, nb, intsrc, bfilt, indx, id;
char *type, *name, *stfdir, *sfile;
{
FILE *fp, *fopfile();
int it, j, nt6, itsh, tshift_st();
char str[256];

float *stmp;
int phase = 0;

if(indx == 0)
   fprintf(stderr,"**** Source time function is '%s'\n",type);

                /* generate source */

if(strncmp("-1",type,2) == 0)
   readsource(sfile,st,dt,nt);
else
   makesource(type,st,dt,nt,tz,nb);

normsrc(st,dt,nt);

if(bfilt)
   {
   itsh = tshift_st(st,dt,nt,flo,tsh);

   stmp = (float *) check_malloc (4*nt*sizeof(float));
   tfilter(st,dt,nt,bfilt,fhi,flo,phase,stmp);
   taper_front(st,nt,itsh);
   free(stmp);

   if(indx == 0)
      {
      fprintf(stderr,"       -a time delay of %10.5e sec has been applied to preserve causality\n",*tsh);
      fprintf(stderr,"       -STF bandpass: %6.2f Hz < f < %6.2f Hz\n",(*fhi),(*flo));
      }
   }

if(intsrc)
   int_src(st,dt,nt);

fprintf(stderr,"           STF[%d] tzero = %f\n",indx,(*tz));

if(id == 0)      /* only output STF to file if root node */
   {
   if(strncmp("-1",type,2) == 0)
      sprintf(str,"%s/%s.stf",stfdir,name);
   else
      sprintf(str,"%s/%s_%s%.2f",stfdir,name,type,*tz);

   fp = fopfile(str,"w");

   if(strncmp("-1",type,2) == 0)
      fprintf(fp,"%s %s\n",name,sfile);
   else
      fprintf(fp,"%s %s\n",name,type);

   fprintf(fp,"%d %f 0 0 0.0 0.0 0.0 0.0\n",nt,(*dt));

   nt6 = nt/6;
   for(it=0;it<nt6;it++)
      {
      for(j=0;j<6;j++)
         fprintf(fp,"%13.5e",st[6*it + j]);

      fprintf(fp,"\n");
      }

   if(6*nt6 != nt)
      {
      for(j=6*nt6;j<nt;j++)
         fprintf(fp,"%13.5e",st[j]);

      fprintf(fp,"\n");
      }

   fclose(fp);
   }
}

taper_frontback(st,nt,itsh)
float *st;
int nt, itsh;
{
float fac, da;
int it;
int itstart, nn;

nn = 3;
itstart = itsh/nn;
while(itstart < 5 && nn > 1)
   {
   nn--;
   itstart = itsh/nn;
   }

da = 3.14159/(float)(itstart);

for(it=0;it<itstart;it++)
   {
   fac = 0.5*(1.0 - cos(it*da));
   st[it] = fac*st[it];
   st[nt-1-it] = fac*st[nt-1-it];
   }
}

taper_front(st,nt,itsh)
float *st;
int nt, itsh;
{
float fac, da;
int it;
int itstart, nn;

nn = 3;
itstart = itsh/nn;
while(itstart < 5 && nn > 1)
   {
   nn--;
   itstart = itsh/nn;
   }

da = 3.14159/(float)(itstart);

for(it=0;it<itstart;it++)
   {
   fac = 0.5*(1.0 - cos(it*da));
   st[it] = fac*st[it];
   }
}

normsrc(s,dt,nt)
float *s;
float *dt;
int nt;
{
float sum;
int it;

sum = 0.0;
for(it=0;it<nt;it++)
   sum = (*dt)*s[it] + sum;

for(it=0;it<nt;it++)
   s[it] = s[it]/sum;
}

int_src(s,dt,nt)
float *s;
float *dt;
int nt;
{
int it;

s[0] = (*dt)*s[0];
for(it=1;it<nt;it++)
   s[it] = (*dt)*s[it] + s[it-1];
}

int tshift_st(s,dt,nt,flo,tsh)
float *s;
float *flo, *dt, *tsh;
int nt;
{
int it, itsh;

itsh = 3.0/((*dt)*(*flo));
*tsh = itsh*(*dt);

for(it=nt-1;it>=itsh;it--)
   s[it] = s[it-itsh];

for(it=0;it<itsh;it++)
   s[it] = 0.0;

return(itsh);
}

tshift_stOLD(s,dt,nt,flo)
float *s;
float *flo, *dt;
int nt;
{
int it, itsh;

itsh = 3.0/((*dt)*(*flo));

fprintf(stderr,"     A time delay of %10.5e sec has been applied to preserve causality\n\n",itsh*(*dt));

for(it=nt-1;it>=itsh;it--)
   s[it] = s[it-itsh];

for(it=0;it<itsh;it++)
   s[it] = 0.0;
}

makesource(type,s,dt,nt,a,nb)
float *s;
float *dt, *a;
int nt, nb;
char *type;
{
double exp(), exp_arg, cos_arg;
float t, t0, t2, tshift, alpha, beta;
int it;
char str[128];
float src_amp = 1.0e+05;
float one = 1.0;

int ite, it2;
int it0, it1;
float da;
int itstart = 1;

if(strncmp("gaus",type,4) == 0)
   {
   tshift = AMP*(*a);
   alpha = -1.0/((*a)*(*a));
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;

      s[it] = src_amp*exp(exp_arg);
      if(nb != 1)
         s[it] = -t*s[it];
      }
   }

else if(strncmp("tri",type,3) == 0)
   {
   ite = (int)((*a)/(*dt) + 0.5);
   if(ite%2)
      ite++;
   it2 = ite/2;
   da = (float)(src_amp/it2);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<it2;it++)
      s[it+itstart] = s[(ite+itstart)-it] = it*da;

   s[it2+itstart] = it2*da;

   for(it=ite+itstart+1;it<nt;it++)
      s[it] = 0.0;
   }

else if(strncmp("atri",type,4) == 0)
   {
   if(type[4] == '\0')
      alpha = 0.5;
   else
      alpha = 0.01*atoi(type+4);

   ite = (int)((*a)/(*dt) + 0.5);
   if(ite%2)
      ite++;

   it2 = ite*alpha;
   alpha = src_amp*(float)(2.0/ite);

   da = alpha/(float)(it2);
   for(it=0;it<=it2;it++)
      s[it+itstart] = s[it+itstart] + it*da;

   da = alpha/(float)(ite-it2);
   for(it=it2+1;it<=ite;it++)
      s[it+itstart] = s[it+itstart] + (ite-it)*da;
   }

else if(strncmp("2triGOOD",type,8) == 0)
   {
   if(type[4] == '\0')  /* set at 10% of total rise time */
      alpha = 0.1*(*a);
   else if(strncmp("-p",type+4,2) == 0)  /* percentage is given after '-p' */
      alpha = 0.01*(*a)*atoi(type+6);
   else
      alpha = atoi(type+4);

   it0 = (int)((alpha)/(*dt) + 0.5);
   if(it0 < 2)
      it0 = 2;
   it1 = (int)((*a)/(*dt) + 0.5);
   if(it1 < 4)
      it1 = 2;

   it2 = it0 + it0/2;

   alpha = 1.0;
   da = alpha/(float)(it0);

   for(it=0;it<it0;it++)
      s[it+itstart] = s[it+itstart] + it*da;

   for(it=it0;it<it2;it++)
      s[it+itstart] = s[it+itstart] + (2*it0-it)*da;

   alpha = 0.5*alpha;
   da = alpha/(float)(it1-it2);

   for(it=it2;it<it1;it++)
      s[it+itstart] = s[it+itstart] + alpha + (it2-it)*da;
   }

else if(strncmp("2tri",type,4) == 0)
   {
   /* DEFAULT: set t0 at 10% of total rise time, h at 0.2*A */

   alpha = 0.1*(*a);
   beta = 0.2;

   if(strncmp("-p",type+4,2) == 0)       /* t0 percentage given after '-p' */
      {
      alpha = 0.01*(*a)*atoi(type+6);

      if(strncmp("-h",type+8,2) == 0)    /* ampl. percentage given after '-h' */
         beta = 0.01*atoi(type+10);
      }

   if(strncmp("-h",type+4,2) == 0)       /* ampl. percentage given after '-h' */
      {
      beta = 0.01*atoi(type+6);

      if(strncmp("-p",type+8,2) == 0)    /* t0 percentage given after '-p' */
         alpha = 0.01*(*a)*atoi(type+10);
      }

   it0 = (int)((alpha)/(*dt) + 0.5);
   if(it0 < 2)
      it0 = 2;
   it1 = (int)((*a)/(*dt) + 0.5);
   if(it1 < 4)
      it1 = 2;

   it2 = (2 - beta)*it0;

   alpha = 1.0;
   da = alpha/(float)(it0);

   for(it=0;it<it0;it++)
      s[it+itstart] = it*da;

   for(it=it0;it<it2;it++)
      s[it+itstart] = (2*it0-it)*da;

   alpha = beta*alpha;
   da = alpha/(float)(it1-it2);

   for(it=it2;it<it1;it++)
      s[it+itstart] = alpha + (it2-it)*da;
   }

else if(strncmp("scec",type,4) == 0)
   {
   alpha = 1.0/((*a)*(*a));

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<nt-itstart;it++)
      {
      exp_arg = it*(*dt);

      s[it+itstart] = exp_arg*alpha*(exp(-exp_arg/(*a)));
      }
   }

else if(strncmp("cos",type,3) == 0)
   {
   alpha = 2.0*PI/(*a);
   ite = (*a)/(*dt);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)*alpha;

      s[it+itstart] = (1.0/(*a))*(one - cos(cos_arg));
      }

   for(it=ite+itstart;it<nt;it++)
      s[it] = 0.0;
   }

else if(strncmp("rick",type,4) == 0)
   {
   tshift = AMP*(*a);
   alpha = -0.5/((*a)*(*a));
   beta = PI/(*a);
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;
      cos_arg = beta*t;

      s[it] = src_amp*exp(exp_arg)*cos(cos_arg);
      }
   }

else if(strncmp("expcos",type,6) == 0)
   {
   ite = (int)((*a)/(*dt) + 0.5);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)/(*a);

      s[it+itstart] = src_amp*cos_arg*(one + cos(cos_arg*PI))*exp(-cos_arg);
      }

   for(it=ite;it<nt;it++)
      s[it] = 0.0;
   }

else if(strncmp("brune",type,5) == 0)
   {
   t0 = 0.2108*(*a);  /* x0 is t95 ~ rise time */
   t0 = 0.1250*(*a);
   t0 = 0.1600*(*a);
   for(it=0;it<nt-itstart;it++)
      {
      exp_arg = it*(*dt)/t0;

      s[it+itstart] = (exp_arg/t0)*(exp(-exp_arg));
      }
   }

else
   {
   fprintf(stderr,"**** STF '%s' is invalid choice, exiting...\n",type);
   exit(-99);
   }
}

makesourceOLD(type,s,dt,nt,a,nb,name)
float *s;
float *dt, *a;
int nt, nb;
char *type, *name;
{
double exp(), exp_arg, cos_arg;
float t, t2, tshift, alpha, beta, gamma;
int it;
char str[128];
float src_amp = 1.0e+05;
float one = 1.0;

int ite, it2;
float da;
int itstart = 1;

if(strncmp("gaus",type,4) == 0)
   {
   tshift = AMP*(*a);
   alpha = -1.0/((*a)*(*a));
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;

      s[it] = src_amp*exp(exp_arg);
      if(nb != 1)
         s[it] = -t*s[it];
      }
   }

else if(strncmp("tri",type,3) == 0)
   {
   ite = (int)((*a)/(*dt) + 0.5);
   if(ite%2)
      ite++;
   it2 = ite/2;
   da = (float)(src_amp/it2);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<it2;it++)
      s[it+itstart] = s[(ite+itstart)-it] = it*da;

   s[it2+itstart] = it2*da;

   for(it=ite+itstart+1;it<nt;it++)
      s[it] = 0.0;
   }

else if(strncmp("atri",type,4) == 0)
   {
   if(type[4] == '\0')
      alpha = 0.5;
   else
      alpha = 0.01*atoi(type+4);

   ite = (int)((*a)/(*dt) + 0.5);
   if(ite%2)
      ite++;

   it2 = ite*alpha;
   alpha = src_amp*(float)(2.0/ite);

   da = alpha/(float)(it2);
   for(it=0;it<=it2;it++)
      s[it+itstart] = s[it+itstart] + it*da;

   da = alpha/(float)(ite-it2);
   for(it=it2+1;it<=ite;it++)
      s[it+itstart] = s[it+itstart] + (ite-it)*da;
   }

else if(strncmp("scec",type,4) == 0)
   {
   alpha = 1.0/((*a)*(*a));

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<nt-itstart;it++)
      {
      exp_arg = it*(*dt);

      s[it+itstart] = exp_arg*alpha*(exp(-exp_arg/(*a)));
      }
   }

else if(strncmp("cos",type,3) == 0)
   {
   alpha = 2.0*PI/(*a);
   ite = (*a)/(*dt);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)*alpha;
 
      s[it+itstart] = src_amp*(one - cos(cos_arg));
      }

   for(it=ite+itstart;it<nt;it++)
      s[it] = 0.0;
   }

else if(strncmp("rick",type,4) == 0)
   {
   tshift = AMP*(*a);
   alpha = -0.5/((*a)*(*a));
   beta = PI/(*a);
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      t2 = t*t;
      exp_arg = alpha*t2;
      cos_arg = beta*t;

      s[it] = src_amp*exp(exp_arg)*cos(cos_arg);
      }
   }

else if(strncmp("gabor",type,5) == 0)
   {
   gamma = 11.0;
   alpha = PI/2.0;
   tshift = 0.45*gamma*(*a);
   beta = 2.0*PI/(*a);
   for(it=0;it<nt;it++)
      {
      t = (it-1)*(*dt) - tshift;
      exp_arg = beta*t/gamma;
      cos_arg = beta*t + alpha;

      s[it] = src_amp*exp(-exp_arg*exp_arg)*cos(cos_arg);
      }
   }

else if(strncmp("expcos",type,6) == 0)
   {
   ite = (int)((*a)/(*dt) + 0.5);

   for(it=0;it<itstart;it++)
      s[it] = 0.0;

   for(it=0;it<ite;it++)
      {
      cos_arg = it*(*dt)/(*a);

      s[it+itstart] = src_amp*cos_arg*(one + cos(cos_arg*PI))*exp(-cos_arg);
      }

   for(it=ite;it<nt;it++)
      s[it] = 0.0;
   }

else
   {
   fprintf(stderr,"**** STF '%s' is invalid choice, exiting...\n",type);
   exit(-99);
   }
}

init_plane_src(pvf,nx,ny,nz,iy,h,t0,dt,ps,nb)
struct planesrc *ps;
float *pvf, *t0, *dt, *h;
int nx, ny, nz, iy, nb;
{
float *vx, *vy, *vz;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float ts, r0vp, r0vs, t2, m, b, z0, d, a;
float a1, rad, yc, xs, difcon;
float dt2, h2;
double cosA, sinA, cosB, sinB, arg;
int ix, iz, i;
float two = 2.0;
float sAsBov, sAcBov, cAov, dxa, dya, dza;
float ux, uy, uz;
float l2m, lam, mu, invvs, invvp;

float bdamp, xdamp, ydamp, zdamp;

invvp = 1.0/ps->pvel;
invvs = 1.0/ps->svel;
l2m = ps->srcl2m;
lam = ps->srclam;
mu = ps->srcmu;

vx =  pvf;
vy =  pvf +   nx*nz;
vz =  pvf + 2*nx*nz;
txx = pvf + 3*nx*nz;
tyy = pvf + 4*nx*nz;
tzz = pvf + 5*nx*nz;
txy = pvf + 6*nx*nz;
txz = pvf + 7*nx*nz;
tyz = pvf + 8*nx*nz;

ts = 5.0*(*t0);
t2 = 1.0/((*t0)*(*t0));

if(ps->incidence < 0.0)
   {
   ps->azimuth = ps->azimuth + 180.0;
   ps->incidence = -(ps->incidence);
   }
while(ps->azimuth >= 360.0)
   ps->azimuth = ps->azimuth - 360.0;
while(ps->azimuth < 0.0)
   ps->azimuth = ps->azimuth + 360.0;

arg = (ps->incidence)*PI/180.0;
cosA = cos(arg);
sinA = sin(arg);
m = sinA/cosA;

arg = -(ps->azimuth)*PI/180.0;
cosB = cos(arg);
sinB = sin(arg);

difcon = 1.0;
yc = iy*(*h)*cosB;

h2 = 0.5*(*h);
dt2 = -0.5*(*dt);

/* do P waves */

if(ps->intercept > 0)
   r0vp = (ps->intercept)*(*h);
else
   r0vp = ts*(ps->pvel) + ps->bnd_zero_pad*(*h);

if(ps->azimuth < 90.0)
   b = (nz + nx*m*sinB)*(*h) - r0vp;
else if(ps->azimuth < 180.0)
   b = (nz + (nx*sinB + ny*cosB)*m)*(*h) - r0vp;
else if(ps->azimuth < 270.0)
   b = (nz + ny*m*cosB)*(*h) - r0vp;
else if(ps->azimuth < 360.0)
   b = nz*(*h) - r0vp;

ux = ps->p*sinA*sinB;
uy = -ps->p*sinA*cosB;
uz = ps->p*cosA;

sAsBov = sinA*sinB*invvp;
sAcBov = sinA*cosB*invvp;
cAov = cosA*invvp;

dxa = sAsBov;
dya = sAcBov;
dza = cAov;

ydamp = 1.0;
if(iy < ps->bnd_zero_pad || iy > ny - 1 - ps->bnd_zero_pad)
   ydamp = 0.0;
else if(iy < ps->bnd_zero_pad + ps->iypzero)
   {
   arg = PI*(iy-ps->bnd_zero_pad)/(float)(ps->iypzero);
   ydamp = 0.5 - 0.5*cos(arg);
   }
else if(iy > ny - 1 - (ps->bnd_zero_pad + ps->iypzero))
   {
   arg = PI*(ny - 1 - ps->bnd_zero_pad - iy)/(float)(ps->iypzero);
   ydamp = 0.5 - 0.5*cos(arg);
   }

for(iz=0;iz<nz;iz++)
   {
   zdamp = 1.0;
   if(iz >= nz - ps->bnd_zero_pad)
      zdamp = 0.0;
   else if(iz > nz - 1 - (ps->bnd_zero_pad + ps->izpzero))
      {
      arg = PI*(nz - 1 - ps->bnd_zero_pad - iz)/(float)(ps->izpzero);
      zdamp = 0.5 - 0.5*cos(arg);
      }

   z0 = iz*(*h);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < ps->bnd_zero_pad || ix > nx - 1 - ps->bnd_zero_pad)
         xdamp = 0.0;
      else if(ix < ps->bnd_zero_pad + ps->ixpzero)
         {
         arg = PI*(ix-ps->bnd_zero_pad)/(float)(ps->ixpzero);
         xdamp = 0.5 - 0.5*cos(arg);
         }
      else if(ix > nx - 1 - (ps->bnd_zero_pad + ps->ixpzero))
         {
         arg = PI*(nx - 1 - ps->bnd_zero_pad - ix)/(float)(ps->ixpzero);
         xdamp = 0.5 - 0.5*cos(arg);
         }

      bdamp = xdamp*ydamp*zdamp;

      i = ix + iz*nx;
      xs = ix*(*h)*sinB;
      rad = xs + yc;

      d = (z0 - b)*cosA + rad*sinA;

      a = d*invvp;
      if(a >= -ts && a <= ts)
         {
         a1 = a + h2*sAsBov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vx[i] += bdamp*difcon*ux*exp(arg);

         a1 = a + h2*sAcBov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vy[i] += bdamp*difcon*uy*exp(arg);

         a1 = a + h2*cAov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vz[i] += bdamp*difcon*uz*exp(arg);

         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txx[i] += bdamp*difcon*(l2m*dxa*ux + lam*(dya*uy + dza*uz))*exp(arg);

         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tyy[i] += bdamp*difcon*(l2m*dya*uy + lam*(dxa*ux + dza*uz))*exp(arg);

         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tzz[i] += bdamp*difcon*(l2m*dza*uz + lam*(dxa*ux + dya*uy))*exp(arg);

         a1 = a + dt2 + h2*(sAsBov + sAcBov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txy[i] += bdamp*difcon*mu*(dya*ux + dxa*uy)*exp(arg);

         a1 = a + dt2 + h2*(sAsBov + cAov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txz[i] += bdamp*difcon*mu*(dza*ux + dxa*uz)*exp(arg);

         a1 = a + dt2 + h2*(sAcBov + cAov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tyz[i] += bdamp*difcon*mu*(dza*uy + dya*uz)*exp(arg);
         }
      }   
   }   

/* do S waves */

if(ps->intercept > 0)
   r0vs = (ps->intercept)*(*h);
else
   r0vs = ts*(ps->svel) + ps->bnd_zero_pad*(*h);

if(ps->azimuth < 90.0)
   b = (nz + nx*m*sinB)*(*h) - r0vs;
else if(ps->azimuth < 180.0)
   b = (nz + (nx*sinB + ny*cosB)*m)*(*h) - r0vs;
else if(ps->azimuth < 270.0)
   b = (nz + ny*m*cosB)*(*h) - r0vs;
else if(ps->azimuth < 360.0)
   b = nz*(*h) - r0vs;

ux = -ps->sh*cosB - ps->sv*cosA*sinB;
uy = ps->sh*sinB - ps->sv*cosA*cosB;
uz = ps->sv*sinA;

sAsBov = sinA*sinB*invvs;
sAcBov = sinA*cosB*invvs;
cAov = cosA*invvs;

dxa = sAsBov;
dya = sAcBov;
dza = cAov;
 
for(iz=0;iz<nz;iz++)
   {
   zdamp = 1.0;
   if(iz >= nz - ps->bnd_zero_pad)
      zdamp = 0.0;
   else if(iz > nz - 1 - (ps->bnd_zero_pad + ps->izpzero))
      {
      arg = PI*(nz - 1 - ps->bnd_zero_pad - iz)/(float)(ps->izpzero);
      zdamp = 0.5 - 0.5*cos(arg);
      }

   z0 = iz*(*h);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < ps->bnd_zero_pad || ix > nx - 1 - ps->bnd_zero_pad)
         xdamp = 0.0;
      else if(ix < ps->bnd_zero_pad + ps->ixpzero)
         {
         arg = PI*(ix-ps->bnd_zero_pad)/(float)(ps->ixpzero);
         xdamp = 0.5 - 0.5*cos(arg);
         }
      else if(ix > nx - 1 - (ps->bnd_zero_pad + ps->ixpzero))
         {
         arg = PI*(nx - 1 - ps->bnd_zero_pad - ix)/(float)(ps->ixpzero);
         xdamp = 0.5 - 0.5*cos(arg);
         }

      bdamp = xdamp*ydamp*zdamp;

      i = ix + iz*nx;
      xs = ix*(*h)*sinB;
      rad = xs + yc;
 
      d = (z0 - b)*cosA + rad*sinA;
 
      a = d*invvs;
      if(a >= -ts && a <= ts)
         {
         a1 = a + h2*sAsBov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vx[i] += bdamp*difcon*ux*exp(arg);
 
         a1 = a + h2*sAcBov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vy[i] += bdamp*difcon*uy*exp(arg);
 
         a1 = a + h2*cAov;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         vz[i] += bdamp*difcon*uz*exp(arg);
 
         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txx[i] += bdamp*difcon*(l2m*dxa*ux + lam*(dya*uy + dza*uz))*exp(arg);

         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tyy[i] += bdamp*difcon*(l2m*dya*uy + lam*(dxa*ux + dza*uz))*exp(arg);

         a1 = a + dt2;
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tzz[i] += bdamp*difcon*(l2m*dza*uz + lam*(dxa*ux + dya*uy))*exp(arg);

         a1 = a + dt2 + h2*(sAsBov + sAcBov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txy[i] += bdamp*difcon*mu*(dya*ux + dxa*uy)*exp(arg);

         a1 = a + dt2 + h2*(sAsBov + cAov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         txz[i] += bdamp*difcon*mu*(dza*ux + dxa*uz)*exp(arg);

         a1 = a + dt2 + h2*(sAcBov + cAov);
         if(nb != 1)
            difcon = -two*a1*t2;
         arg = -a1*a1*t2;
         tyz[i] += bdamp*difcon*mu*(dza*uy + dya*uz)*exp(arg);
         }
      }   
   }
}
/*SN: inserted ny1: to convert global system to local system */

void get_dc_par(struct pntsrcs *psrc,float *dt,float *tz,int ny1)
{
struct doublecouple *dc;
float sum, invdt, *tsrc, *rtime, rt, *rtime_uniq;
int isrc, j, nstf, *rtflag;
float modelrot = 0.0;
 
dc = &psrc->doublec;

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc        = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime       = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((psrc->nsource)*sizeof(int));

dc->strike  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->dip     = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->rake    = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->momwgts = (float *) check_malloc ((psrc->nsource)*sizeof(float));

dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   rtflag[isrc] = 0;
   rtime[isrc] = *tz;
   }
 
mstpar("xsrc","vd",psrc->ix);
mstpar("ysrc","vd",psrc->iy);
mstpar("zsrc","vd",psrc->iz);
mstpar("moment","f",&dc->moment);
mstpar("strike","vf",dc->strike);
mstpar("dip","vf",dc->dip);
mstpar("rake","vf",dc->rake);
getpar("rtime","vf",rtime);

getpar("modelrot","f",&modelrot);

getpar("relative_slip","d",&psrc->relative_slip);
getpar("absolute_slip","d",&psrc->absolute_slip);

/*global to local conversion */
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   dc->strike[isrc] = dc->strike[isrc] - modelrot;
   psrc->iy[isrc] = psrc->iy[isrc] - ny1;
   }

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc);
   mstpar("momwgts","vf",dc->momwgts);

   invdt = 1.0/(*dt);
   sum = 0.0;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      dc->itsrc[isrc] = tsrc[isrc]*invdt;
      sum = sum + dc->momwgts[isrc];
      }

   sum = 1.0/sum;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      dc->momwgts[isrc] = sum*dc->momwgts[isrc];
   }
else
   {
   dc->itsrc[0] = 0;
   dc->momwgts[0] = 1.0;
   }

free(tsrc);
free(rtime);
free(rtime_uniq);
free(rtflag);
}

void get_pointmt_par(struct pntsrcs *psrc,float *dt,float *tz,int ny1)
{
struct doublecouple *dc;
struct momenttensor *mt;
float sum, invdt, *tsrc, *rtime, rt, *rtime_uniq;
int isrc, j, nstf, *rtflag;

float modelrot = 0.0;

float *Mrr, *Mtt, *Mpp, *Mrt, *Mrp, *Mtp;
float *Mnn, *Mee, *Mdd, *Mne, *Mnd, *Med;
float arg, cosA, sinA, cos2A, sin2A;

float rperd = RPERD;
 
dc = &psrc->doublec;

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc        = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime       = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((psrc->nsource)*sizeof(int));

dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   rtflag[isrc] = 0;
   rtime[isrc] = *tz;
   }
 
mstpar("xsrc","vd",psrc->ix);
mstpar("ysrc","vd",psrc->iy);
mstpar("zsrc","vd",psrc->iz);
getpar("rtime","vf",rtime);

getpar("modelrot","f",&modelrot);

/*global to local conversion */
for(isrc=0;isrc<psrc->nsource;isrc++)
   psrc->iy[isrc] = psrc->iy[isrc] - ny1;

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc);

   invdt = 1.0/(*dt);
   for(isrc=0;isrc<psrc->nsource;isrc++)
      dc->itsrc[isrc] = tsrc[isrc]*invdt;
   }
else
   dc->itsrc[0] = 0;

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));

Mrr = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mpp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtp = (float *) check_malloc ((psrc->nsource)*sizeof(float));

Mnn = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mee = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mdd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mne = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mnd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Med = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   Mrr[isrc] = 1.0e-15;
   Mtt[isrc] = 1.0e-15;
   Mpp[isrc] = 1.0e-15;
   Mrt[isrc] = 1.0e-15;
   Mrp[isrc] = 1.0e-15;
   Mtp[isrc] = 1.0e-15;

   Mnn[isrc] = 1.0e-15;
   Mee[isrc] = 1.0e-15;
   Mdd[isrc] = 1.0e-15;
   Mne[isrc] = 1.0e-15;
   Mnd[isrc] = 1.0e-15;
   Med[isrc] = 1.0e-15;
   }

getpar("Mrr","vf",Mrr);
getpar("Mtt","vf",Mtt);
getpar("Mpp","vf",Mpp);
getpar("Mrt","vf",Mrt);
getpar("Mrp","vf",Mrp);
getpar("Mtp","vf",Mtp);

getpar("Mnn","vf",Mnn);
getpar("Mee","vf",Mee);
getpar("Mdd","vf",Mdd);
getpar("Mne","vf",Mne);
getpar("Mnd","vf",Mnd);
getpar("Med","vf",Med);

if(Mrr[0] > 1.0e-14 || -Mrr[0] > 1.0e-14)
   {
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      Mnn[isrc] = Mtt[isrc];
      Mee[isrc] = Mpp[isrc];
      Mdd[isrc] = Mrr[isrc];
      Mne[isrc] = -Mtp[isrc];
      Mnd[isrc] = Mrt[isrc];
      Med[isrc] = -Mrp[isrc];
      }
   }

arg = (90.0 + modelrot)*rperd;
cosA = cos(arg);
sinA = sin(arg);

cos2A = cosA*cosA - sinA*sinA;
sin2A = 2.0*sinA*cosA;

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   mt = &(psrc->momten[isrc]);
   mt->mxx = (float *) check_malloc (sizeof(float));
   mt->myy = (float *) check_malloc (sizeof(float));
   mt->mzz = (float *) check_malloc (sizeof(float));
   mt->mxy = (float *) check_malloc (sizeof(float));
   mt->mxz = (float *) check_malloc (sizeof(float));
   mt->myz = (float *) check_malloc (sizeof(float));

   mt->mxx[0] = Mnn[isrc]*cosA*cosA + Mee[isrc]*sinA*sinA + Mne[isrc]*sin2A;
   mt->myy[0] = Mnn[isrc]*sinA*sinA + Mee[isrc]*cosA*cosA - Mne[isrc]*sin2A;
   mt->mzz[0] = Mdd[isrc];
   mt->mxy[0] = 0.5*(Mee[isrc] - Mnn[isrc])*sin2A + Mne[isrc]*cos2A;
   mt->mxz[0] = Mnd[isrc]*cosA + Med[isrc]*sinA;
   mt->myz[0] = -Mnd[isrc]*sinA + Med[isrc]*cosA;
   }

free(tsrc);
free(rtime);
free(rtime_uniq);
free(rtflag);

free(Mrr);
free(Mtt);
free(Mpp);
free(Mrt);
free(Mrp);
free(Mtp);

free(Mnn);
free(Mee);
free(Mdd);
free(Mne);
free(Mnd);
free(Med);
}

void get_ff_par(struct pntsrcs *psrc,float *dt,float *h,int fs,float *tz)
{
FILE *fpr, *fopfile();
struct doublecouple *dc;
float sum, invdt, *tsrc, momwt, *rtime, rt, *rtime_uniq;
int i, isrc, nsrc_file, nstf, j, n, *rtflag;
char faultfile[512];
char string[512];

float zap = 0.0;

/* Set zap=-1.0 10/07/97, RWG ===>
   
   Originally, I used the variable 'zap' to check for zero weight sources,
   which I then did not use in the simulations, ie., they were not
   stored in the source arrays in order to save memory.  Although, this
   is the most efficient calculation method, it destroys some information
   which can be used to check total slip amounts, average slip, etc.
   Therefore, I decided to specify 'zap=-1.0' to ensure that all sources,
   including zero weight sources will be used in the simulation.  This
   change has absolutely no effect on the outcome of the calculations.

*/

zap = -1.0;

dc = &psrc->doublec;
 
getpar("relative_slip","d",&psrc->relative_slip);
getpar("absolute_slip","d",&psrc->absolute_slip);

if(psrc->absolute_slip == 1)
   {
   psrc->relative_slip = 0;
   mstpar("area","f",&psrc->area);
   }
else if(psrc->relative_slip == 1)
   {
   mstpar("moment","f",&dc->moment);
   mstpar("area","f",&psrc->area);
   }
else
   mstpar("moment","f",&dc->moment);

mstpar("faultfile","s",faultfile);

    /* First, determine number of non-zero source in faultfile */

fpr = fopfile(faultfile,"r");
fgets(string,512,fpr);
sscanf(string,"%d",&nsrc_file);

isrc = 0;
for(i=0;i<nsrc_file;i++)
   {
   fgets(string,512,fpr);
   sscanf(string,"%*f %*f %*f %*f %*f %*f %f %*f",&momwt);

   if(momwt > zap) /* don't store zero-weight sources */
      isrc++;
   }
fclose(fpr);
psrc->nsource = isrc;

    /* Allocate memory only for non-zero sources */
    /* add 1 for reading in next value */

psrc->xs = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
psrc->ys = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
psrc->zs = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
psrc->ix = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));

tsrc        = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
rtime       = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));

dc->strike  = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
dc->dip     = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
dc->rake    = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));
dc->momwgts = (float *) check_malloc ((1 + psrc->nsource)*sizeof(float));

dc->itsrc   = (int *) check_malloc ((1 + psrc->nsource)*sizeof(int));

    /* Now, read in sources, skipping those with zero weight */

fpr = fopfile(faultfile,"r");
fgets(string,512,fpr);
sscanf(string,"%d",&nsrc_file);

isrc = 0;
for(i=0;i<nsrc_file;i++)
   {
   fgets(string,512,fpr);
   n = sscanf(string,"%f %f %f %f %f %f %f %f %f",&psrc->xs[isrc],
					&psrc->ys[isrc],
					&psrc->zs[isrc],
					&dc->strike[isrc],
					&dc->dip[isrc],
					&dc->rake[isrc],
					&dc->momwgts[isrc],
					&tsrc[isrc],
                                        &rt);
   rtflag[isrc] = 0;
   if(n == 9)
      rtime[isrc] = rt;
   else
      rtime[isrc] = *tz;

   psrc->ix[isrc] = (int)(psrc->xs[isrc]/(*h));
   psrc->iy[isrc] = (int)(psrc->ys[isrc]/(*h));

   if(fs)
      psrc->iz[isrc] = (int)(psrc->zs[isrc]/(*h) + 0.5) + 1;
   else
      psrc->iz[isrc] = (int)(psrc->zs[isrc]/(*h));

   if(dc->momwgts[isrc] > zap) /* don't store zero-weight sources */
      isrc++;
   }
fclose(fpr);

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

if(psrc->nsource > 1)
   {
   invdt = 1.0/(*dt);
   sum = 0.0;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      dc->itsrc[isrc] = tsrc[isrc]*invdt;
      sum = sum + dc->momwgts[isrc];
      }

   if(psrc->absolute_slip == 0)
      {
      sum = 1.0/sum;
      for(isrc=0;isrc<psrc->nsource;isrc++)
         dc->momwgts[isrc] = sum*dc->momwgts[isrc];
      }
   }
else
   {
   dc->itsrc[0] = 0;
   dc->momwgts[0] = 1.0;
   }

free(tsrc);
free(rtime);
free(rtime_uniq);
free(rtflag);
}

void init_dc(struct pntsrcs *psrc,float *dt)
{
struct doublecouple *dc;
struct momenttensor *mt;
float cosD, sinD, cosL, sinL, cosA, sinA;
float cos2D, sin2D, cos2A, sin2A;
float sfac, arg, hcm, invh4;
int i;

float half = 0.5;
float two = 2.0;
float normf = 1.0e-05; /* convert dyne-cm to Nt-km */
 
dc = &psrc->doublec;

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));

hcm = normf/(dc->h);
hcm = (*dt)*hcm*hcm*hcm;

invh4 = hcm*(dc->moment)*(normf);

fprintf(stderr,"**** Moment-tensor source\n");
fprintf(stderr,"            total sources= %d\n",psrc->nsource);

for(i=0;i<psrc->nsource;i++)
   {
   mt = &(psrc->momten[i]);
   mt->mxx = (float *) check_malloc (sizeof(float));
   mt->myy = (float *) check_malloc (sizeof(float));
   mt->mzz = (float *) check_malloc (sizeof(float));
   mt->mxy = (float *) check_malloc (sizeof(float));
   mt->mxz = (float *) check_malloc (sizeof(float));
   mt->myz = (float *) check_malloc (sizeof(float));

   sfac = invh4*dc->momwgts[i];

   arg = dc->strike[i]*RPERD;
   cosA = cos(arg);
   sinA = sin(arg);

   cos2A = cosA*cosA - sinA*sinA;
   sin2A = two*sinA*cosA;

   arg = dc->dip[i]*RPERD;
   cosD = cos(arg);
   sinD = sin(arg);

   cos2D = cosD*cosD - sinD*sinD;
   sin2D = two*sinD*cosD;

   arg = dc->rake[i]*RPERD;
   cosL = cos(arg);
   sinL = sin(arg);
 
   mt->mxx[0] = sfac*(sinD*cosL*sin2A - sin2D*sinL*cosA*cosA);
   mt->myy[0] = -sfac*(sinD*cosL*sin2A + sin2D*sinL*sinA*sinA);
   mt->mzz[0] = sfac*sin2D*sinL;

   mt->mxy[0] = -sfac*(sinD*cosL*cos2A + half*sin2D*sinL*sin2A);
   mt->mxz[0] = -sfac*(cosD*cosL*sinA - cos2D*sinL*cosA);
   mt->myz[0] = sfac*(cosD*cosL*cosA + cos2D*sinL*sinA);

   if(dc->strike[i] < -999.0)  /* explosion */
      {
      mt->mxx[0] = sfac;
      mt->myy[0] = sfac;
      mt->mzz[0] = sfac;

      mt->mxy[0] = 0.0;
      mt->mxz[0] = 0.0;
      mt->myz[0] = 0.0;
      }
   }
}

void init_pointmt(struct pntsrcs *psrc,float *dt,float *h)
{
struct momenttensor *mt;
float tmom, hcm;
int i;

float norm20 = 1.0e-20;
float normf = 1.0e-05; /* convert dyne-cm to Nt-km */

hcm = normf/(*h);
hcm = (*dt)*hcm*hcm*hcm;

fprintf(stderr,"**** Point moment-tensor source\n");
fprintf(stderr,"            total sources= %d\n\n",psrc->nsource);

for(i=0;i<psrc->nsource;i++)
   {
   mt = &(psrc->momten[i]);

   tmom = sqrt( 0.5*(norm20*mt->mxx[0])*(norm20*mt->mxx[0])
              + 0.5*(norm20*mt->myy[0])*(norm20*mt->myy[0])
              + 0.5*(norm20*mt->mzz[0])*(norm20*mt->mzz[0])
                  + (norm20*mt->mxy[0])*(norm20*mt->mxy[0])
                  + (norm20*mt->mxz[0])*(norm20*mt->mxz[0])
                  + (norm20*mt->myz[0])*(norm20*mt->myz[0]) );

   tmom = tmom/norm20;
   fprintf(stderr,"              %d)      Mo= %.5e dyne-cm\n",i,tmom);
   fprintf(stderr,"                       Mw= %.2f\n",(2.0/3.0)*log10(tmom)-10.7);

   mt->mxx[0] = normf*hcm*mt->mxx[0];
   mt->myy[0] = normf*hcm*mt->myy[0];
   mt->mzz[0] = normf*hcm*mt->mzz[0];
   mt->mxy[0] = normf*hcm*mt->mxy[0];
   mt->mxz[0] = normf*hcm*mt->mxz[0];
   mt->myz[0] = normf*hcm*mt->myz[0];
   }
}

tfilter(y,dt,nt,order,fhi,flo,phase,space)
float *y, *dt, *fhi, *flo, *space;
int nt, order, phase;
{
struct complex *q, *p, *tmpptr;
double trig_arg;
int j, it, i;
float wplo, wphi, cosA, sinA;
float fnyq;
float are, aim;
float one = 1.0;
float two = 2.0;

p = (struct complex *) space;
q = (struct complex *) (space + 2*nt);
cxzero(p,nt);
cxzero(q,nt);

for(it=0;it<nt;it++)
   p[it].re = y[it];

fnyq = one/(two*(*dt));

if((*fhi))              
   {
   if((*fhi) < fnyq)
      wphi = tan(PI*(*fhi)*(*dt));
   else
      wphi = 1.0e+10;
   }   

if((*flo) < fnyq)
   wplo = tan(PI*(*flo)*(*dt));

           /*  forward pass  */ 

for(j=0;j<order;j++)            
   {
   trig_arg = 0.5*(two*j + one)*PI/order;
   cosA = cos(trig_arg);
   sinA = -sin(trig_arg);

   if((*flo) < fnyq)    /* low-pass filter */
      {
      are = wplo*sinA;
      aim = -wplo*cosA;
      lp_filter(q,p,nt,&are,&aim,1);

      tmpptr = p; p = q; q = tmpptr; 
      }

   if((*fhi) > 0.0)    /* high-pass filter */
      {
      are = wphi*sinA;
      aim = -wphi*cosA;
      hp_filter(q,p,nt,&are,&aim,1);

      tmpptr = p; p = q; q = tmpptr; 
      }
   }

              /*  reverse pass to obtain zero-phase response  */

if(phase == 0)                                                   
   {
   for(j=0;j<order;j++)
      {
      trig_arg = 0.5*(two*j + one)*PI/order;
      cosA = cos(trig_arg);
      sinA = -sin(trig_arg);

      if((*flo) < fnyq)    /* low-pass filter */
         {
         are = wplo*sinA;
         aim = -wplo*cosA;
         lp_filter(q,p,nt,&are,&aim,-1);

         tmpptr = p; p = q; q = tmpptr;  
         }

      if((*fhi) > 0.0)    /* high-pass filter */
         {
         are = wphi*sinA;
         aim = -wphi*cosA;
         hp_filter(q,p,nt,&are,&aim,-1);

         tmpptr = p; p = q; q = tmpptr;  
         }
      }
   }

for(it=0;it<nt;it++)
   y[it] = p[it].re;
}

cxzero(p,n)
struct complex *p;
int n;
{
float zap = 0.0;
int i;

for(i=0;i<n;i++)
   p[i].re = p[i].im = zap;
}

hp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

n1 = n - 1;    
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;   

q[k].re = (are*p[k].re - aim*p[k].im);
q[k].im = (are*p[k].im + aim*p[k].re);
 
while(n1--)
   {
   tmpre = (p[i].re - p[k].re) + bre*q[k].re - bim*q[k].im;
   tmpim = (p[i].im - p[k].im) + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;                         
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;                     
   k = i - sgn;
   }
}

lp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim, cre, cim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

cre = -(*alpha);
cim = -(*beta);

n1 = n - 1;     
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;   

tmpre = cre*p[k].re - cim*p[k].im;
tmpim = cre*p[k].im + cim*p[k].re;

q[k].re = are*tmpre - aim*tmpim;   
q[k].im = are*tmpim + aim*tmpre;

while(n1--)                      
   {
   tmpre = cre*(p[i].re + p[k].re) - cim*(p[i].im + p[k].im)
      + bre*q[k].re - bim*q[k].im;
   tmpim = cre*(p[i].im + p[k].im) + cim*(p[i].re + p[k].re)
      + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;                     
   k = i - sgn;
   }
}

int global2local_srcs(struct pntsrcs *srcs,struct nodeinfo *ni,int fs)
{
int k, isrc;

isrc = 0;
for(k=0;k<srcs->nsource;k++)
   {
   srcs->ix[k] = srcs->ix[k] - ni->nx1;
   srcs->iy[k] = srcs->iy[k] - ni->ny1;
   srcs->iz[k] = srcs->iz[k] - ni->nz1;

   if(srcs->ix[k] >= ni->ixminus && srcs->ix[k] <= ni->ixplus
         && srcs->iy[k] >= ni->iyminus && srcs->iy[k] <= ni->iyplus
            && srcs->iz[k] >= ni->izminus && srcs->iz[k] <= ni->izplus)
      {
      srcs->ix[isrc] = srcs->ix[k];
      srcs->iy[isrc] = srcs->iy[k];
      srcs->iz[isrc] = srcs->iz[k];

      if(fs && ni->minusId_z < 0 && srcs->iz[isrc] < 1)
         srcs->iz[isrc] = 1;

      isrc++;
      }
   }

return(isrc);
}

globalSource2Local(psrc,ny1)
struct pntsrcs *psrc;
int ny1;
{
  int isrc;
  for(isrc=0;isrc<psrc->nsource;isrc++)
  {
	psrc->iy[isrc]-=ny1;
  }
}

slip2mom(psrc,h,mds,medf,nx,ny1,ny2,ny,nz,outf,dt,xfrmt,tfrmt,nodeType)
struct pntsrcs *psrc;
struct modelstorage *mds;
float *medf, *dt, *h;
int nx, ny1, ny2, ny, nz, nodeType;
char *outf, *xfrmt, *tfrmt;
{
FILE *fpw, *fopfile();
struct doublecouple *dc;
float *mptr, *mu;
int iy, j, isrc, iloc, nloc;
int tnloc;
float mscl, cosA, sinA;
double sum, avgs, tdbl;
float *xs, *ys, *zs, *xf, *yf, *zf;
float *rgdloc, *rgdsrc, *slipss, *slipds, *slipv,  *ts;
float *rt;
float *stk, *dip, *rak;
float mfit, dd, xx, yy, zz;
float subf_area;

int iyleft, iyright, localiy;

char frmt[128];

float conv = 1.0e-20;
float dperr = 1.0/RPERD;

xs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ys = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
xf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
yf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ts = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
stk = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dip = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rgdloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rgdsrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipss = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipds = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   slipss[isrc] = 0.0;
   slipds[isrc] = 0.0;
   ts[isrc] = 1.0e+15;
   }

dc = &psrc->doublec;

init_model_seek(mds,MEDIA_FIELD);
mptr = medf;

mfit = 0.0;
nloc = 0;

iyleft = 0;
iyright = ny-1;

if(nodeType != PARL_HEAD)
   iyleft = 2;
if(nodeType != PARL_TAIL)
   iyright = ny-3;

iyleft = iyleft + ny1;
iyright = iyright + ny1;

/*SN:change iy to localiy, which is local index specific to the machine now */
/* because we take one extra plane , which is indexed as -1 in body/tail node */

/*
   nloc = local (ie, within this CPU) # of unique source locations
   tnloc = global # of unique source locations as resolved with other CPUs
*/

for(localiy=0;localiy<ny;localiy++)
   {
   iy = localiy + ny1;

   mptr = reed_model(mds,mptr,MEDIA_FIELD);
   mu = mptr + 12*nx*nz;  /* raw mu */

   if(iy >= iyleft && iy <= iyright)
      {
      for(isrc=0;isrc<psrc->nsource;isrc++) /* find sources in this y-plane */
         {
         if(iy == psrc->iy[isrc])
            {
	    j = psrc->ix[isrc] + psrc->iz[isrc]*nx;
	    cosA = cos(dc->rake[isrc]*RPERD);
	    sinA = sin(dc->rake[isrc]*RPERD);

            xx = psrc->xs[isrc]; /* use exact location */
            yy = psrc->ys[isrc];
            zz = psrc->zs[isrc];

	    rgdsrc[isrc] = mu[j];

            for(iloc=0;iloc<nloc;iloc++)
               {
               if(xx == xs[iloc] && yy == ys[iloc] && zz == zs[iloc])
                  break;
               }

	    if(iloc == nloc)  /* new source location */
	       {
	       xs[iloc] = xx;
	       ys[iloc] = yy;
	       zs[iloc] = zz;
	       rgdloc[iloc] = rgdsrc[isrc];
	       nloc++;

               xf[iloc] = (*h)*psrc->ix[isrc];
               yf[iloc] = (*h)*psrc->iy[isrc];
               zf[iloc] = (*h)*(psrc->iz[isrc] - 1);

               xx = xs[iloc] - xf[iloc];
               yy = ys[iloc] - yf[iloc];
               zz = zs[iloc] - zf[iloc];

               dd = sqrt(xx*xx + yy*yy + zz*zz);
               if(dd > mfit)
                  mfit = dd;
	       }

	    slipss[iloc] = slipss[iloc] + cosA*dc->momwgts[isrc];
	    slipds[iloc] = slipds[iloc] + sinA*dc->momwgts[isrc];

	    stk[iloc] = dc->strike[isrc];
	    dip[iloc] = dc->dip[isrc];

            if(((*dt)*dc->itsrc[isrc]) < ts[iloc]) /* use minimum tstart */
	       ts[iloc] = (*dt)*dc->itsrc[isrc];

	    rt[iloc] = psrc->rtime[psrc->stfindx[isrc]];

	 /* convert to moment (kindof) */
	    dc->momwgts[isrc] = rgdsrc[isrc]*dc->momwgts[isrc];
	    } /*if iy=isrc */
         } /*for isrc=.... to */
      }
   } /* for iy=... to */

if(nloc > 0)
   {
   rak = (float *) check_malloc (nloc*sizeof(float));
   slipv = (float *) check_malloc (nloc*sizeof(float));
   }

mpi_global_val(&nloc,&tnloc,"Nloc",1,MPI_INT,MPI_SUM);

subf_area = psrc->area/(float)(tnloc);

/*
   rak = final rake value for vector slip
   slipv = total relative vector slip for each source location
   sum = total relative moment for all sources
*/

sum = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   {
   if(slipss[iloc] != 0.0)
      {
      rak[iloc] = atan(slipds[iloc]/slipss[iloc]);
      if(slipss[iloc] < 0.0)
	 rak[iloc] = rak[iloc] + PI;
      }
   else if(slipds[iloc] < 0.0)
      rak[iloc] = -0.5*PI;
   else
      rak[iloc] = 0.5*PI;

   rak[iloc] = rak[iloc]*dperr;

   slipv[iloc] = sqrt(slipss[iloc]*slipss[iloc] + slipds[iloc]*slipds[iloc]);
   sum = sum + slipv[iloc]*rgdloc[iloc];
   }

/* here, the moment calculation should be more sophisticated
        should collect all the sum!
*/

tdbl = sum;
mpi_global_val(&tdbl,&sum,"Momsum",1,MPI_DOUBLE,MPI_SUM);

sum = 1.0/sum;

/* normalize each source moment by total relative moment */
for(isrc=0;isrc<psrc->nsource;isrc++)
   dc->momwgts[isrc] = sum*dc->momwgts[isrc];

/* compute scaling factor and scale relative slip to absolute slip */

if(psrc->relative_slip == 1)
   mscl = (conv*(dc->moment))*sum/(subf_area);

else if(psrc->absolute_slip == 1)
   {
   mscl = 1.0;
   dc->moment = (1.0/sum)*(1.0/conv)*subf_area;
   }

avgs = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   {
   slipv[iloc] = mscl*slipv[iloc];
   avgs = avgs + slipv[iloc];
   }

tdbl = avgs;
mpi_global_val(&tdbl,&avgs,"Avgs",1,MPI_DOUBLE,MPI_SUM);

avgs = avgs/(tnloc);
fprintf(stderr,"     total moment= %13.5e (dyne-cm)\n",dc->moment);
fprintf(stderr,"     average slip= %13.5le (cm)\n",avgs);

if(outf[0] != '\0')
   {
   fpw = fopfile(outf,"w");

   fprintf(fpw,"%d global source locations; %d local source locations\n",tnloc,nloc);
   fprintf(fpw,"%10.4e average slip (cm); %10.4e total moment (dyne-cm)\n",avgs,dc->moment);

   if(nloc)
      {
      fprintf(fpw,"%10.4e maximum misfit to FD grid\n",mfit);

      fprintf(fpw,"  x(exact)   y(exact)   z(exact)    x(grid)    y(grid)    z(grid)  strike    dip   rake  slip (cm)   time (s)  trise (s)\n");
      fprintf(fpw,"------------------------------------------------------------------------------------------------------------------------\n");

      sprintf(frmt,"%%%s %%%s %%%s %%%s %%%s %%%s  %%6.1f %%6.1f %%6.1f %%%s %%%s %%%s\n",xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,tfrmt,tfrmt,tfrmt);

      for(iloc=0;iloc<nloc;iloc++)
         {
         fprintf(fpw,frmt,xs[iloc],
                             ys[iloc],
                             zs[iloc],
                             xf[iloc],
                             yf[iloc],
                             zf[iloc],
                             stk[iloc],
                             dip[iloc],
                             rak[iloc],
                             slipv[iloc],
                             ts[iloc],
                             rt[iloc]);
         }
      }
   fclose(fpw);
   }

globalSource2Local(psrc,ny1);     /* transform to local indexing */
init_model_seek(mds,MEDIA_FIELD);

free(xs);
free(ys);
free(zs);
free(xf);
free(yf);
free(zf);
free(ts);
free(stk);
free(dip);
free(rgdloc);
free(rgdsrc);
free(slipss);
free(slipds);

if(nloc > 0)
   {
   free(rak);
   free(slipv);
   }
}

void get_srf_par(struct pntsrcs *psrc,struct runparams *rp,int bfilt,float *flo,float *fhi,float *tsh)
{
FILE *fpw, *fpr, *fopfile();
int iyleft, iyright;
float modellon, modellat, modelrot;
float *stf1, *stf2, *stf3, *space;
float lon, lat, dep, stk, dip, rake, area, tinit, dt_stf, slip1, slip2, slip3;
float de, dn, fnt, da, xx, yy;
int nt1, nt2, nt3, ntmax, ntpad, ntrsmp, ntout, resamp;
int isrc, j, jj, kk, it, nsrc;
char str[1024], pword[32], *rchar;

int gnt;
float tol = 1.0e-02;

struct momenttensor *mt;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
float latavg;
float cosR, sinR, kmlon, kmlat;
float invh;

double g0, b0;
int xy2ll = 0;
int ll2xy = 1;

int inside = 1;
float zap = 0.0;

float nyq_perc = 1.0;
float tap_perc = 0.0;

int order = 4;

int itsh;
int phase = 0;

iyleft = rp->ny1 + 2;
if(rp->ny1 == 0)
   iyleft = 0;

iyright = rp->ny2 - 3;
if(rp->ny2 == rp->globny)
   iyright = rp->ny2 - 1;

psrc->modelrot = rp->modelrot;

invh = 1.0/(rp->h);

cosR = cos(rp->modelrot*rperd);
sinR = sin(rp->modelrot*rperd);

itsh = 0;
if(bfilt)
   {
   itsh = 3.0/((rp->dt)*(*flo));
   *tsh = itsh*(rp->dt);
   }

fpr = fopfile(psrc->faultfile,"r");

/* 09/22/05
   For now, simply assume ASCII input and scan down to find "POINTS" line
*/

fgets(str,1024,fpr);
sscanf(str,"%s",pword);
while(strncmp(pword,"POINTS",6) != 0)
   {
   rchar = fgets(str,1024,fpr);
   if(rchar == NULL)
      {
      fprintf(stderr,"Unexpected EOF in %s, exiting...\n",psrc->faultfile);
      exit(-99);
      }
   sscanf(str,"%s",pword);
   }

sscanf(str,"%*s %d",&nsrc);

psrc->momten = (struct momenttensor *) check_malloc (nsrc*sizeof(struct momenttensor));

psrc->xs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ys = (float *) check_malloc (nsrc*sizeof(float));
psrc->zs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ix = (int *) check_malloc (nsrc*sizeof(int));
psrc->iy = (int *) check_malloc (nsrc*sizeof(int));
psrc->iz = (int *) check_malloc (nsrc*sizeof(int));
psrc->it = (int *) check_malloc (nsrc*sizeof(int));

space = NULL;
stf1 = NULL;
stf2 = NULL;
stf3 = NULL;

isrc = 0;
for(j=0;j<nsrc;j++)
   {
   if(j==nsrc || fgets(str,1024,fpr) == NULL)
      break;

   sscanf(str,"%f %f %f %f %f %f %f %f",&lon,
                                           &lat,
                                           &dep,
                                           &stk,
                                           &dip,
                                           &area,
                                           &tinit,
                                           &dt_stf);
   fgets(str,1024,fpr);
   sscanf(str,"%f %f %d %f %d %f %d",&rake,
                                        &slip1,
                                        &nt1,
                                        &slip2,
                                        &nt2,
                                        &slip3,
                                        &nt3);

   ntmax = nt1;
   if(nt2 > ntmax)
      ntmax = nt2;
   if(nt3 > ntmax)
      ntmax = nt3;

   if(ntmax)
      {
      stf1 = (float *)check_realloc((void *)stf1,ntmax*sizeof(float));
      stf2 = (float *)check_realloc((void *)stf2,ntmax*sizeof(float));
      stf3 = (float *)check_realloc((void *)stf3,ntmax*sizeof(float));
      }

   zero(stf1,ntmax);
   zero(stf2,ntmax);
   zero(stf3,ntmax);

   for(it=0;it<nt1;it++)
      fscanf(fpr,"%f",&stf1[it]);
   for(it=0;it<nt2;it++)
      fscanf(fpr,"%f",&stf2[it]);
   for(it=0;it<nt3;it++)
      fscanf(fpr,"%f",&stf3[it]);

/* get rouge newline character */
   if(nt1 || nt2 || nt3)
      fgets(str,1024,fpr);

/* convert source location to grid coordinates */

   if(rp->geoproj == 0)
      {
      de = rp->kmlon*(lon - rp->modellon);
      dn = rp->kmlat*(lat - rp->modellat);

      psrc->xs[isrc] = (de*rp->cosR - dn*rp->sinR) - rp->xshift;
      psrc->ys[isrc] = (-de*rp->sinR - dn*rp->cosR) - rp->yshift;
      }
   else if(rp->geoproj == 1)
      {
      gcproj(&xx,&yy,&lon,&lat,&rp->erad,&rp->g0,&rp->b0,rp->amat,rp->ainv,ll2xy);

      psrc->xs[isrc] = xx;
      psrc->ys[isrc] = yy;
      }

   psrc->ix[isrc] = (int)(psrc->xs[isrc]*invh + 0.5);
   psrc->iy[isrc] = (int)(psrc->ys[isrc]*invh + 0.5);

   psrc->zs[isrc] = dep;
   if(rp->freesurf)
      psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5) + 1;
   else
      psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5);

   inside = 1;
   if((psrc->ix[isrc] < 0 || psrc->ix[isrc] >= rp->nx) ||
         (psrc->iy[isrc] < 0 || psrc->iy[isrc] >= rp->globny) ||
            (psrc->iz[isrc] < 1 || psrc->iz[isrc] >= rp->nz))
      {
      fprintf(stderr,"**** source point %10.4f %10.4f is outside model box\n",lon,lat);
      fprintf(stderr,"     ix= %d iy= %d iz= %d\n",psrc->ix[isrc],psrc->iy[isrc],psrc->iz[isrc]);
      inside = 0;
      }

   psrc->it[isrc] = (int)(tinit/(rp->dt) + 0.5);

/* check if point is within model sub region for this node */

   if((inside == 1) && (psrc->iy[isrc] >= iyleft && psrc->iy[isrc] <= iyright))
      {
      mt = &(psrc->momten[isrc]);

      mt->lon = lon;
      mt->lat = lat;
      mt->dep = dep;
      mt->stk = stk;
      mt->dip = dip;
      mt->rake = rake;
      mt->area = area;
      mt->tinit = tinit;
      mt->slip1 = slip1;
      mt->slip2 = slip2;
      mt->slip3 = slip3;
      mt->dt_stf = dt_stf;
      mt->nt_stf = ntmax;

      if(ntmax == 0)
	 {
         mt->nt = 0;
	 mt->stf1 = NULL;
	 mt->stf2 = NULL;
	 mt->stf3 = NULL;

	 mt->mxx = NULL;
	 mt->myy = NULL;
	 mt->mzz = NULL;
	 mt->mxy = NULL;
	 mt->mxz = NULL;
	 mt->myz = NULL;
	 }
      else if(((rp->dt) > 0.999*dt_stf && (rp->dt) < 1.001*dt_stf))
         {
	 mt->nt = ntmax + 2*itsh;
	 mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	 mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	 for(it=0;it<itsh;it++)
	    {
	    mt->stf1[it] = 0.0;
	    mt->stf2[it] = 0.0;
	    mt->stf3[it] = 0.0;
	    mt->stf1[it+ntmax+itsh] = 0.0;
	    mt->stf2[it+ntmax+itsh] = 0.0;
	    mt->stf3[it+ntmax+itsh] = 0.0;
	    }
	 for(it=0;it<ntmax;it++)
	    {
	    mt->stf1[it+itsh] = stf1[it];
	    mt->stf2[it+itsh] = stf2[it];
	    mt->stf3[it+itsh] = stf3[it];
	    }
	 }
      else  /* need to resample STFs */
	 {
	 ntpad = 2*ntmax;
	 fnt = ntpad*dt_stf/(rp->dt);
	 gnt = (int)(fnt + 0.5);

         while(nt_tol(fnt,gnt) > tol)
            {
            ntpad++;
            fnt = ntpad*dt_stf/(rp->dt);
	    gnt = (int)(fnt + 0.5);
            }

         ntrsmp = (int)(fnt);
	 ntout = (int)(ntmax*dt_stf/(rp->dt));

         if((rp->dt) < dt_stf)
            {
            resamp = 1;

            if(ntout > ntrsmp)
               ntout = ntrsmp;

            space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
            stf1 = (float *)check_realloc((void *)stf1,2*ntrsmp*sizeof(float));
            stf2 = (float *)check_realloc((void *)stf2,2*ntrsmp*sizeof(float));
            stf3 = (float *)check_realloc((void *)stf3,2*ntrsmp*sizeof(float));
            }
         else
            {
            resamp = -1;

            if(ntout > ntpad)
               ntout = ntpad;

            space = (float *) check_realloc (space,2*ntpad*sizeof(float));
            stf1 = (float *)check_realloc(stf1,2*ntpad*sizeof(float));
            stf2 = (float *)check_realloc(stf2,2*ntpad*sizeof(float));
            stf3 = (float *)check_realloc(stf3,2*ntpad*sizeof(float));
            }

         resample(stf1,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
         resample(stf2,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
         resample(stf3,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);

	 mt->nt = ntout + 2*itsh;
	 mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	 mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	 for(it=0;it<itsh;it++)
	    {
	    mt->stf1[it] = 0.0;
	    mt->stf2[it] = 0.0;
	    mt->stf3[it] = 0.0;
	    mt->stf1[it+ntout+itsh] = 0.0;
	    mt->stf2[it+ntout+itsh] = 0.0;
	    mt->stf3[it+ntout+itsh] = 0.0;
	    }
	 for(it=0;it<ntout;it++)
	    {
	    mt->stf1[it+itsh] = stf1[it];
	    mt->stf2[it+itsh] = stf2[it];
	    mt->stf3[it+itsh] = stf3[it];
	    }
	 }

      if(bfilt && mt->nt != 0)
         {
         space = (float *) check_realloc (space,4*mt->nt*sizeof(float));

         tfilter(mt->stf1,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
         tfilter(mt->stf2,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
         tfilter(mt->stf3,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);

         taper_frontback(mt->stf1,mt->nt,itsh);
         taper_frontback(mt->stf2,mt->nt,itsh);
         taper_frontback(mt->stf3,mt->nt,itsh);
         }

      isrc++;
      }
   }
fclose(fpr);

free(space);
free(stf1);
free(stf2);
free(stf3);

psrc->nsource = isrc;

if(psrc->nsource)
   {
   psrc->momten = (struct momenttensor *)check_realloc(psrc->momten,(psrc->nsource)*sizeof(struct momenttensor));

   psrc->xs = (float *)check_realloc(psrc->xs,(psrc->nsource)*sizeof(float));
   psrc->ys = (float *)check_realloc(psrc->ys,(psrc->nsource)*sizeof(float));
   psrc->zs = (float *)check_realloc(psrc->zs,(psrc->nsource)*sizeof(float));
   psrc->ix = (int *)check_realloc(psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *)check_realloc(psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *)check_realloc(psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *)check_realloc(psrc->it,(psrc->nsource)*sizeof(int));
   }
else
   {
   free(psrc->momten);
   free(psrc->xs);
   free(psrc->ys);
   free(psrc->zs);
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   }
}

void init_mt(struct pntsrcs *psrc,float *h,struct modelstorage *mds,float *medf,int nx,int ny1,int ny2,int ny,int nz,char *outf,float *dt,char *xfrmt,char *tfrmt,int nodeType)
{
FILE *fpw, *fopfile();
struct doublecouple *dc;
float *mptr, *lam, *mu, *rho;
int iy, j, isrc, iloc, nloc;
int tnloc, it;
float  cosS, sinS, cosD, sinD, cosL, sinL, arg;
double sum, avgs, tdbl;
float *xs, *ys, *zs, *xf, *yf, *zf;
float *slipss, *slipds, *slipnm, *slipv,  *ts;
float *l2mloc, *lamloc, *muloc;
float *l2msrc, *lamsrc, *musrc;
float *stf1, *stf2, *stf3;
float *mxx, *myy, *mzz, *mxy, *mxz, *myz;
float *Mloc_xx, *Mloc_yy, *Mloc_zz, *Mloc_xy, *Mloc_xz, *Mloc_yz;
float l2mcnv, lamcnv, mucnv, vx, vy, vz, ux, uy, uz;
float *rt, *area;
float *stk, *dip, *rak, moment;
float mfit, dd, xx, yy, zz;
struct momenttensor *mt;
float conv_fac;
char frmt[128];
int iyleft, iyright, localiy;

float dperr = 1.0/RPERD;
float rt2inv = 0.707106781;

/*
   Recall, global units are:

   wavefield velocity is-
    vx,vy,vz       => cm/s

   seismic velocity and density are-
    Vp,Vs     => km/s
    density   => gm/cm^3

   Lame parameters and stress are-
    Tnn,Tij,lambda,mu => km*km*gm/(s*s*cm*cm*cm)
                      = 1e+09 Nt/(m^2)
                      = 1e+09 Pa
                      = 1 GPa

   distance and time measures-
    h         => km
    dt        => s

   sources are CGS-
    moment    => dyne-cm
    force     => dyne
    slip      => cm
    slip-rate => cm/s

   The SRF slip parameters are given as:

    slip      => cm
    slip-rate => cm/s
    area      => cm*cm

   'Slip' and 'slip-rate' are consistent with source CGS units, but 'area'
   which is a distance measure, needs to be converted to km^2.

   Thus, we apply a conversion factor of 1.0e-10.
*/

/*
   Also add factors here of-
    dt for time integration
    1/h^3 for volume normalization
*/

conv_fac = 1.0e-10*(*dt)/((*h)*(*h)*(*h));

xs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ys = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
xf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
yf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ts = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
stk = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dip = (float *) check_malloc ((psrc->nsource)*sizeof(float));
l2mloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
lamloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
muloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
l2msrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
lamsrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
musrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipss = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipds = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipnm = (float *) check_malloc ((psrc->nsource)*sizeof(float));
area = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xx = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_yy = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_zz = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xy = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xz = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_yz = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   slipss[isrc] = 0.0;
   slipds[isrc] = 0.0;
   slipnm[isrc] = 0.0;
   ts[isrc] = 1.0e+15;

   Mloc_xx[isrc] = 0.0;
   Mloc_yy[isrc] = 0.0;
   Mloc_zz[isrc] = 0.0;
   Mloc_xy[isrc] = 0.0;
   Mloc_xz[isrc] = 0.0;
   Mloc_yz[isrc] = 0.0;
   }

init_model_seek(mds,MEDIA_FIELD);
mptr = medf;

mfit = 0.0;
nloc = 0;

iyleft = 0;
iyright = ny-1;

if(nodeType != PARL_HEAD)
   iyleft = 2;
if(nodeType != PARL_TAIL)
   iyright = ny-3;

iyleft = iyleft + ny1;
iyright = iyright + ny1;

/*
   nloc = local (ie, within this CPU) # of unique source locations
   tnloc = global # of unique source locations as resolved with other CPUs
*/

for(localiy=0;localiy<ny;localiy++)
   {
   iy = localiy + ny1;

   mptr = reed_model(mds,mptr,MEDIA_FIELD);
   rho = mptr + 7*nx*nz;  /* raw density */
   lam = mptr + 11*nx*nz;  /* raw lambda */
   mu  = mptr + 12*nx*nz;  /* raw mu */

   if(iy >= iyleft && iy <= iyright)
      {
      for(isrc=0;isrc<psrc->nsource;isrc++) /* find sources in this y-plane */
         {
         if(iy == psrc->iy[isrc])
            {
            mt = &(psrc->momten[isrc]);

	    j = psrc->ix[isrc] + psrc->iz[isrc]*nx;
	    l2msrc[isrc] = lam[j] + 2.0*mu[j];
	    lamsrc[isrc] = lam[j];
	    musrc[isrc] = mu[j];

            mxx = mt->mxx;
            myy = mt->myy;
            mzz = mt->mzz;
            mxy = mt->mxy;
            mxz = mt->mxz;
            myz = mt->myz;

	    stf1 = mt->stf1;
	    stf2 = mt->stf2;
	    stf3 = mt->stf3;

            arg = ((mt->stk) - 90.0 - (psrc->modelrot))*RPERD;
            cosS = cos(arg);
            sinS = sin(arg);

            cosD = cos(mt->dip*RPERD);
            sinD = sin(mt->dip*RPERD);

            cosL = cos(mt->rake*RPERD);
            sinL = sin(mt->rake*RPERD);

	    /* unit normals */
	    vx = -sinD*sinS;
	    vy =  sinD*cosS;
	    vz = -cosD;

	    /* convert from CGS to Nt-km per unit volume */
	    l2mcnv = l2msrc[isrc]*mt->area*conv_fac;
	    lamcnv = lamsrc[isrc]*mt->area*conv_fac;
	    mucnv = musrc[isrc]*mt->area*conv_fac;

	    for(it=0;it<mt->nt;it++)
	       {
	       ux = -(stf3[it]*sinD - cosD*(stf1[it]*sinL + stf2[it]*cosL))*sinS
	               + (stf1[it]*cosL - stf2[it]*sinL)*cosS;

	       uy =  (stf3[it]*sinD - cosD*(stf1[it]*sinL + stf2[it]*cosL))*cosS
	               + (stf1[it]*cosL - stf2[it]*sinL)*sinS;

	       uz = -stf3[it]*cosD - (stf1[it]*sinL + stf2[it]*cosL)*sinD;

	       mxx[it] = l2mcnv*vx*ux + lamcnv*(vy*uy + vz*uz);
	       myy[it] = l2mcnv*vy*uy + lamcnv*(vx*ux + vz*uz);
	       mzz[it] = l2mcnv*vz*uz + lamcnv*(vx*ux + vy*uy);

	       mxy[it] = mucnv*(vx*uy + vy*ux);
	       mxz[it] = mucnv*(vx*uz + vz*ux);
	       myz[it] = mucnv*(vy*uz + vz*uy);
	       }

            free(mt->stf1);
            free(mt->stf2);
            free(mt->stf3);

            xx = psrc->xs[isrc]; /* use exact location */
            yy = psrc->ys[isrc];
            zz = psrc->zs[isrc];

            for(iloc=0;iloc<nloc;iloc++)
               {
               if(xx == xs[iloc] && yy == ys[iloc] && zz == zs[iloc])
                  break;
               }

	    if(iloc == nloc)  /* new source location */
	       {
	       xs[iloc] = xx;
	       ys[iloc] = yy;
	       zs[iloc] = zz;
	       l2mloc[iloc] = l2msrc[isrc];
	       lamloc[iloc] = lamsrc[isrc];
	       muloc[iloc] = musrc[isrc];
	       nloc++;

               xf[iloc] = (*h)*psrc->ix[isrc];
               yf[iloc] = (*h)*psrc->iy[isrc];
               zf[iloc] = (*h)*(psrc->iz[isrc] - 1);

               xx = xs[iloc] - xf[iloc];
               yy = ys[iloc] - yf[iloc];
               zz = zs[iloc] - zf[iloc];

               dd = sqrt(xx*xx + yy*yy + zz*zz);
               if(dd > mfit)
                  mfit = dd;
	       }

	    slipss[iloc] = slipss[iloc] + cosL*mt->slip1 - sinL*mt->slip2;
	    slipds[iloc] = slipds[iloc] + sinL*mt->slip1 + cosL*mt->slip2;
	    slipnm[iloc] = slipnm[iloc] + mt->slip3;

	    ux = -(mt->slip3*sinD - cosD*(mt->slip1*sinL + mt->slip2*cosL))*sinS
	               + (mt->slip1*cosL - mt->slip2*sinL)*cosS;

	    uy =  (mt->slip3*sinD - cosD*(mt->slip1*sinL + mt->slip2*cosL))*cosS
	               + (mt->slip1*cosL - mt->slip2*sinL)*sinS;

	    uz = -mt->slip3*cosD - (mt->slip1*sinL + mt->slip2*cosL)*sinD;

	    Mloc_xx[iloc] = Mloc_xx[iloc]
	                  + l2mloc[iloc]*vx*ux*mt->area
			  + lamloc[iloc]*(vy*uy + vz*uz)*mt->area;
	    Mloc_yy[iloc] = Mloc_yy[iloc]
	                  + l2mloc[iloc]*vy*uy*mt->area
			  + lamloc[iloc]*(vx*ux + vz*uz)*mt->area;
	    Mloc_zz[iloc] = Mloc_zz[iloc]
	                  + l2mloc[iloc]*vz*uz*mt->area
			  + lamloc[iloc]*(vx*ux + vy*uy)*mt->area;

	    Mloc_xy[iloc] = Mloc_xy[iloc] + muloc[iloc]*(vx*uy + vy*ux)*mt->area;
	    Mloc_xz[iloc] = Mloc_xz[iloc] + muloc[iloc]*(vx*uz + vz*ux)*mt->area;
	    Mloc_yz[iloc] = Mloc_yz[iloc] + muloc[iloc]*(vy*uz + vz*uy)*mt->area;

	    area[iloc] = mt->area;

	    stk[iloc] = mt->stk;
	    dip[iloc] = mt->dip;

            if(((*dt)*psrc->it[isrc]) < ts[iloc]) /* use minimum tstart */
	       ts[iloc] = (*dt)*psrc->it[isrc];

	    rt[iloc] = (*dt)*mt->nt;
	    }
         }
      }
   }

if(nloc > 0)
   {
   rak = (float *) check_malloc (nloc*sizeof(float));
   slipv = (float *) check_malloc (nloc*sizeof(float));
   }

mpi_global_val(&nloc,&tnloc,"Nloc",1,MPI_INT,MPI_SUM);

/*
   rak = final rake value for vector slip
   slipv = total relative vector slip for each source location
   sum = total relative moment for all sources
*/

sum = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   {
   if(slipss[iloc] != 0.0)
      {
      rak[iloc] = atan(slipds[iloc]/slipss[iloc]);
      if(slipss[iloc] < 0.0)
	 rak[iloc] = rak[iloc] + PI;
      }
   else if(slipds[iloc] < 0.0)
      rak[iloc] = -0.5*PI;
   else
      rak[iloc] = 0.5*PI;

   rak[iloc] = rak[iloc]*dperr;

   slipv[iloc] = sqrt(slipss[iloc]*slipss[iloc] + slipds[iloc]*slipds[iloc]);

/* purely double couple
   sum = sum + slipv[iloc]*muloc[iloc]*area[iloc]*1.0e+10;
*/

/*
   full moment tensor: Mo = sqrt(0.5*[MM]),
   [MM]=double dot product of moment tensor, Dahlen and Tromp
*/
   sum = sum + rt2inv*sqrt(Mloc_xx[iloc]*Mloc_xx[iloc]
                         + Mloc_yy[iloc]*Mloc_yy[iloc]
	                 + Mloc_zz[iloc]*Mloc_zz[iloc]
	                 + 2.0*Mloc_xy[iloc]*Mloc_xy[iloc]
	                 + 2.0*Mloc_xz[iloc]*Mloc_xz[iloc]
	                 + 2.0*Mloc_yz[iloc]*Mloc_yz[iloc])*1.0e+10;
   }

tdbl = sum;
mpi_global_val(&tdbl,&sum,"Momsum",1,MPI_DOUBLE,MPI_SUM);

moment = sum;

avgs = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   avgs = avgs + slipv[iloc];

tdbl = avgs/(tnloc);
mpi_global_val(&tdbl,&avgs,"Avgs",1,MPI_DOUBLE,MPI_SUM);

fprintf(stderr,"     total moment= %13.5e (dyne-cm)\n",moment);
fprintf(stderr,"     average slip= %13.5le (cm)\n",avgs);
fflush(stderr);

if(outf[0] != '\0')
   {
   fpw = fopfile(outf,"w");

   fprintf(fpw,"%d global source locations; %d local source locations\n",tnloc,nloc);
   fprintf(fpw,"%10.4e average slip (cm); %10.4e total moment (dyne-cm)\n",avgs,moment);

   if(nloc)
      {
      fprintf(fpw,"%10.4e maximum misfit to FD grid\n",mfit);

      fprintf(fpw,"  x(exact)   y(exact)   z(exact)    x(grid)    y(grid)    z(grid)  strike    dip   rake  slip (cm)   time (s)  trise (s) area(km*km) mu() ix iy iz\n");
      fprintf(fpw,"------------------------------------------------------------------------------------------------------------------------\n");

      sprintf(frmt,"%%%s %%%s %%%s %%%s %%%s %%%s  %%6.1f %%6.1f %%6.1f %%%s %%%s %%%s %%13.5e %%13.5e %%5d %%5d %%5d\n",xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,tfrmt,tfrmt,tfrmt);

      for(iloc=0;iloc<nloc;iloc++)
         {
         fprintf(fpw,frmt,xs[iloc],
                             ys[iloc],
                             zs[iloc],
                             xf[iloc],
                             yf[iloc],
                             zf[iloc],
                             stk[iloc],
                             dip[iloc],
                             rak[iloc],
                             slipv[iloc],
                             ts[iloc],
                             rt[iloc],
                             area[iloc],
                             muloc[iloc],
                             (int)(xf[iloc]/(*h) + 0.5),
                             (int)(yf[iloc]/(*h) + 0.5),
                             (int)(zf[iloc]/(*h) + 1.5));
         }
      }
   fclose(fpw);
   }

globalSource2Local(psrc,ny1);     /* transform to local indexing */
init_model_seek(mds,MEDIA_FIELD);

free(xs);
free(ys);
free(zs);
free(xf);
free(yf);
free(zf);
free(ts);
free(stk);
free(dip);
free(l2mloc);
free(lamloc);
free(muloc);
free(l2msrc);
free(lamsrc);
free(musrc);
free(slipss);
free(slipds);
free(slipnm);

free(rt);
free(area);
free(Mloc_xx);
free(Mloc_yy);
free(Mloc_zz);
free(Mloc_xy);
free(Mloc_xz);
free(Mloc_yz);

if(nloc > 0)
   {
   free(rak);
   free(slipv);
   }
}

void write_mt(struct pntsrcs *psrc,float *dt,float *stf,int ny,int nt,int eflag,int myid)
{
struct momenttensor *mt;
float *mxx, *myy, *mzz, *mxy, *mxz, *myz;
int it, isrc, fd;
FILE *fpw, *fopfile();
char outfile[1024];

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(psrc->iy[isrc] > 0 && psrc->iy[isrc] < ny)       
      {
      mt = &(psrc->momten[isrc]);

      mxx = mt->mxx;
      myy = mt->myy;
      mzz = mt->mzz;
      mxy = mt->mxy;
      mxz = mt->mxz;
      myz = mt->myz;

      if(psrc->ffault == 2)
         {
	 sprintf(outfile,"Stf/mxx.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"mxx mxx %13.5e\n",mxx[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",mxx[it]);

	 fclose(fpw);
	 close(fd);

	 sprintf(outfile,"Stf/myy.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"myy myy %13.5e\n",myy[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",myy[it]);

	 fclose(fpw);
	 close(fd);

	 sprintf(outfile,"Stf/mzz.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"mzz mzz %13.5e\n",mzz[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",mzz[it]);

	 fclose(fpw);
	 close(fd);

	 sprintf(outfile,"Stf/mxy.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"mxy mxy %13.5e\n",mxy[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",mxy[it]);

	 fclose(fpw);
	 close(fd);

	 sprintf(outfile,"Stf/mxz.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"mxz mxz %13.5e\n",mxz[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",mxz[it]);

	 fclose(fpw);
	 close(fd);

	 sprintf(outfile,"Stf/myz.%.5d.%.4d",myid,isrc);
	 fd = croptrfile(outfile);
	 fpw = fdopen(fd,"w");

         fprintf(fpw,"myz myz %13.5e\n",myz[1]*(*dt));
         fprintf(fpw,"%d %.5e\n",mt->nt,(*dt));
         for(it=0;it<mt->nt;it++)
            fprintf(fpw,"%13.5e\n",myz[it]);

	 fclose(fpw);
	 close(fd);
         }
      else
         {
         fprintf(stderr,"mxx mxx %13.5e\n",mxx[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",mxx[0]*stf[it]);

         fprintf(stderr,"myy myy %13.5e\n",myy[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",myy[0]*stf[it]);

         fprintf(stderr,"mzz mzz %13.5e\n",mzz[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",mzz[0]*stf[it]);

         fprintf(stderr,"mxy mxy %13.5e\n",mxy[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",mxy[0]*stf[it]);

         fprintf(stderr,"mxz mxz %13.5e\n",mxz[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",mxz[0]*stf[it]);

         fprintf(stderr,"myz myz %13.5e\n",myz[0]);
         fprintf(stderr,"%d %.5e\n",nt,(*dt));
         for(it=0;it<nt;it++)
            fprintf(stderr,"%13.5e\n",myz[0]*stf[it]);
         }
      }
   }

mpi_global_val(&psrc->iy[0],&it,"iy",1,MPI_INT,MPI_MAX);

if(eflag != 0)
   {
   fflush(stderr);
   exit(-1);
   }
}

int readsource(char *sfile,float *st,float *dt,int nt)
{
FILE *fpr, *fopfile();
float dt0, *st0, *space, fnt;
int it, nt0;
char str[512];

int gnt, ntpad, ntrsmp, ntout, resamp;

float nyq_perc = 1.0;
float tap_perc = 0.0;

int order = 4;

float tol = 1.0e-02;
int ierr = 0;

space = NULL;

fpr = fopfile(sfile,"r");

fgets(str,512,fpr);
fgets(str,512,fpr);
sscanf(str,"%d %f",&nt0,&dt0);

st0 = (float *) check_malloc (nt0*sizeof(float));

for(it=0;it<nt0;it++)
   fscanf(fpr,"%f",&st0[it]);

fclose(fpr);

zero(st,nt);
ntout = nt0;

if((*dt) <= 0.999*dt0 || (*dt) >= 1.001*dt0)		/* need to resample STFs */
   {
   ntpad = 2*nt0;
   fnt = ntpad*dt0/(*dt);
   gnt = (int)(fnt + 0.5);

   while(nt_tol(fnt,gnt) > tol)
      {
      ntpad++;
      fnt = ntpad*dt0/(*dt);
      gnt = (int)(fnt + 0.5);
      }

   ntrsmp = (int)(fnt);
   ntout = (int)(nt0*dt0/(*dt));

   if((*dt) < dt0)
      {
      resamp = 1;

      if(ntout > ntrsmp)
         ntout = ntrsmp;

      space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
      st0 = (float *)check_realloc((void *)st0,2*ntrsmp*sizeof(float));
      }
   else
      {
      resamp = -1;

      if(ntout > ntpad)
         ntout = ntpad;

      space = (float *) check_realloc (space,2*ntpad*sizeof(float));
      st0 = (float *)check_realloc(st0,2*ntpad*sizeof(float));
      }

   resample(st0,nt0,&dt0,resamp,ntpad,ntrsmp,dt,space,order,&nyq_perc,&tap_perc);
   }

if(ntout > nt)
   ntout = nt;

for(it=0;it<ntout;it++)
   st[it] = st0[it];

free(space);
free(st0);

return(ierr);
}

void get_dc_parP3(struct pntsrcs *psrc,float *tz,float *dt)
{
struct doublecouple *dc;
float sum, invdt, *tsrc, *rtime, *rtime_uniq;
int isrc, j, nstf, *rtflag;

dc = &psrc->doublec;

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc        = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime       = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((psrc->nsource)*sizeof(int));

dc->strike  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->dip     = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->rake    = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->momwgts = (float *) check_malloc ((psrc->nsource)*sizeof(float));

dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   rtflag[isrc] = 0;
   rtime[isrc] = *tz;
   }
 
mstpar("xsrc","vd",psrc->ix);
mstpar("ysrc","vd",psrc->iy);
mstpar("zsrc","vd",psrc->iz);
mstpar("moment","f",&dc->moment);
mstpar("strike","vf",dc->strike);
mstpar("dip","vf",dc->dip);
mstpar("rake","vf",dc->rake);
getpar("rtime","vf",rtime);

getpar("modelrot","f",&psrc->modelrot);

getpar("relative_slip","d",&psrc->relative_slip);
getpar("absolute_slip","d",&psrc->absolute_slip);

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc);
   mstpar("momwgts","vf",dc->momwgts);

   invdt = 1.0/(*dt);
   sum = 0.0;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      dc->itsrc[isrc] = tsrc[isrc]*invdt;
      sum = sum + dc->momwgts[isrc];
      }

   sum = 1.0/sum;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      dc->momwgts[isrc] = sum*dc->momwgts[isrc];
   }
else
   {
   dc->itsrc[0] = 0;
   dc->momwgts[0] = 1.0;
   }

free(tsrc);
free(rtime);
free(rtime_uniq);
free(rtflag);
}

void get_pointmt_parP3(struct pntsrcs *psrc,float *tz,float *dt)
{
struct doublecouple *dc;
struct momenttensor *mt;
float sum, invdt, *tsrc, *rtime, *rtime_uniq;
int isrc, j, nstf, *rtflag;

float *Mrr, *Mtt, *Mpp, *Mrt, *Mrp, *Mtp;
float *Mnn, *Mee, *Mdd, *Mne, *Mnd, *Med;
float arg, cosA, sinA, cos2A, sin2A;

float rperd = RPERD;

dc = &psrc->doublec;

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc        = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime       = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((psrc->nsource)*sizeof(int));

dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   rtflag[isrc] = 0;
   rtime[isrc] = *tz;
   }
 
mstpar("xsrc","vd",psrc->ix);
mstpar("ysrc","vd",psrc->iy);
mstpar("zsrc","vd",psrc->iz);
getpar("rtime","vf",rtime);

getpar("modelrot","f",&psrc->modelrot);

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc);

   invdt = 1.0/(*dt);
   for(isrc=0;isrc<psrc->nsource;isrc++)
      dc->itsrc[isrc] = tsrc[isrc]*invdt;
   }
else
   dc->itsrc[0] = 0;

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));

Mrr = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mpp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtp = (float *) check_malloc ((psrc->nsource)*sizeof(float));

Mnn = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mee = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mdd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mne = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mnd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Med = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   Mrr[isrc] = 1.0e-15;
   Mtt[isrc] = 1.0e-15;
   Mpp[isrc] = 1.0e-15;
   Mrt[isrc] = 1.0e-15;
   Mrp[isrc] = 1.0e-15;
   Mtp[isrc] = 1.0e-15;

   Mnn[isrc] = 1.0e-15;
   Mee[isrc] = 1.0e-15;
   Mdd[isrc] = 1.0e-15;
   Mne[isrc] = 1.0e-15;
   Mnd[isrc] = 1.0e-15;
   Med[isrc] = 1.0e-15;
   }

getpar("Mrr","vf",Mrr);
getpar("Mtt","vf",Mtt);
getpar("Mpp","vf",Mpp);
getpar("Mrt","vf",Mrt);
getpar("Mrp","vf",Mrp);
getpar("Mtp","vf",Mtp);

getpar("Mnn","vf",Mnn);
getpar("Mee","vf",Mee);
getpar("Mdd","vf",Mdd);
getpar("Mne","vf",Mne);
getpar("Mnd","vf",Mnd);
getpar("Med","vf",Med);

if(Mrr[0] > 1.0e-14 || -Mrr[0] > 1.0e-14)
   {
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      Mnn[isrc] = Mtt[isrc];
      Mee[isrc] = Mpp[isrc];
      Mdd[isrc] = Mrr[isrc];
      Mne[isrc] = -Mtp[isrc];
      Mnd[isrc] = Mrt[isrc];
      Med[isrc] = -Mrp[isrc];
      }
   }

arg = (90.0 + psrc->modelrot)*rperd;
cosA = cos(arg);
sinA = sin(arg);

cos2A = cosA*cosA - sinA*sinA;
sin2A = 2.0*sinA*cosA;

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   mt = &(psrc->momten[isrc]);
   mt->mxx = (float *) check_malloc (sizeof(float));
   mt->myy = (float *) check_malloc (sizeof(float));
   mt->mzz = (float *) check_malloc (sizeof(float));
   mt->mxy = (float *) check_malloc (sizeof(float));
   mt->mxz = (float *) check_malloc (sizeof(float));
   mt->myz = (float *) check_malloc (sizeof(float));

   mt->mxx[0] = Mnn[isrc]*cosA*cosA + Mee[isrc]*sinA*sinA + Mne[isrc]*sin2A;
   mt->myy[0] = Mnn[isrc]*sinA*sinA + Mee[isrc]*cosA*cosA - Mne[isrc]*sin2A;
   mt->mzz[0] = Mdd[isrc];
   mt->mxy[0] = 0.5*(Mee[isrc] - Mnn[isrc])*sin2A + Mne[isrc]*cos2A;
   mt->mxz[0] = Mnd[isrc]*cosA + Med[isrc]*sinA;
   mt->myz[0] = -Mnd[isrc]*sinA + Med[isrc]*cosA;
   }

free(tsrc);
free(rtime);
free(rtime_uniq);
free(rtflag);

free(Mrr);
free(Mtt);
free(Mpp);
free(Mrt);
free(Mrp);
free(Mtp);

free(Mnn);
free(Mee);
free(Mdd);
free(Mne);
free(Mnd);
free(Med);
}

void get_bforce_parP3(struct pntsrcs *psrc,float *tz,struct runparamsP3 *rpars)
{
float *rtime, *rtime_uniq;
int isrc, j, nstf, *rtflag;
struct nodeinfo *ni;

ni = &(rpars->ni);

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));

rtime       = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_uniq  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtflag      = (int *) check_malloc ((psrc->nsource)*sizeof(int));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   rtflag[isrc] = 0;
   rtime[isrc] = *tz;
   }

mstpar("xsrc","vd",psrc->ix);
mstpar("ysrc","vd",psrc->iy);
mstpar("zsrc","vd",psrc->iz);
getpar("rtime","vf",rtime);

nstf = 0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   if(rtflag[isrc] == 0)
      {
      psrc->stfindx[isrc] = nstf;
      rtime_uniq[nstf] = rtime[isrc];

      for(j=isrc+1;j<psrc->nsource;j++)
         {
         if(rtime[j] == rtime[isrc])
            {
            rtflag[j] = 1;
            psrc->stfindx[j] = nstf;
            }
         }

      nstf++;
      }
   }

psrc->nstf = nstf;
psrc->rtime = (float *) check_malloc ((psrc->nstf)*sizeof(float));

for(j=0;j<psrc->nstf;j++)
   psrc->rtime[j] = rtime_uniq[j];

psrc->fxsrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
psrc->fysrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
psrc->fzsrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));

/*
   First check to see if moment source components are
   specified for reciprocal GF calculations.  If they are
   then load them into the corresponding body force
   components.  If not, then get the body forces directly.

   Note that the moment sources are NOT normalized by the
   factor (1/h).  This is implicit.  Later, when the output
   strains are calculated, the factor of (h) that should
   be applied is not.  These two factors cancel one another,
   thus there is no need to put them in.  See also wrout()
   in iofunc.c
*/

getpar("xmom","f",&rpars->xmom);
getpar("ymom","f",&rpars->ymom);
getpar("zmom","f",&rpars->zmom);

if(rpars->xmom != 0.0 || rpars->ymom != 0.0 || rpars->zmom != 0.0)
   {
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      psrc->fxsrc[isrc] = rpars->xmom;
      psrc->fysrc[isrc] = rpars->ymom;
      psrc->fzsrc[isrc] = rpars->zmom;
      }
   }
else
   {
   mstpar("fxsrc","vf",psrc->fxsrc);
   mstpar("fysrc","vf",psrc->fysrc);
   mstpar("fzsrc","vf",psrc->fzsrc);
   }

free(rtime);
free(rtime_uniq);
free(rtflag);
}

void init_dcP3(struct pntsrcs *psrc,struct runparamsP3 *rpars)
{
struct doublecouple *dc;
struct momenttensor *mt;
float cosD, sinD, cosL, sinL, cosA, sinA;
float cos2D, sin2D, cos2A, sin2A;
float sfac, arg, hcm, invh4, tmom;
float dt, h;
int i, js;

float half = 0.5;
float two = 2.0;
float normf = 1.0e-05; /* convert dyne-cm to Nt-km */
struct nodeinfo *ni;

ni = &(rpars->ni);
dc = &psrc->doublec;
h = rpars->h;
dt = rpars->dt;

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));

hcm = normf/h;
hcm = dt*hcm*hcm*hcm;

invh4 = hcm*(dc->moment)*(normf);

fprintf(stderr,"**** Point double-couple source\n");
fprintf(stderr,"            total sources (all nodes)= %d\n\n",psrc->nsource);
fprintf(stderr,"            sources active in this node:\n");

js = 0;
for(i=0;i<psrc->nsource;i++)
   {
/* extend to +1 for averaging (nov 2018) */
   if(psrc->ix[i] >= ni->ixminus && psrc->ix[i] <= (ni->ixplus + 1)
         && psrc->iy[i] >= ni->iyminus && psrc->iy[i] <= (ni->iyplus + 1)
               && psrc->iz[i] >= ni->izminus && psrc->iz[i] <= (ni->izplus + 1))
      {

         /*global to local conversion */
      dc->strike[i] = dc->strike[i] - psrc->modelrot;
      psrc->ix[i] = psrc->ix[i] - ni->nx1;
      psrc->iy[i] = psrc->iy[i] - ni->ny1;
      psrc->iz[i] = psrc->iz[i] - ni->nz1;

      if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[i] < 1)
         psrc->iz[i] = 1;

      js++;

      mt = &(psrc->momten[i]);
      mt->mxx = (float *) check_malloc (sizeof(float));
      mt->myy = (float *) check_malloc (sizeof(float));
      mt->mzz = (float *) check_malloc (sizeof(float));
      mt->mxy = (float *) check_malloc (sizeof(float));
      mt->mxz = (float *) check_malloc (sizeof(float));
      mt->myz = (float *) check_malloc (sizeof(float));

      sfac = invh4*dc->momwgts[i];

      arg = dc->strike[i]*RPERD;
      cosA = cos(arg);
      sinA = sin(arg);

      cos2A = cosA*cosA - sinA*sinA;
      sin2A = two*sinA*cosA;

      arg = dc->dip[i]*RPERD;
      cosD = cos(arg);
      sinD = sin(arg);

      cos2D = cosD*cosD - sinD*sinD;
      sin2D = two*sinD*cosD;

      arg = dc->rake[i]*RPERD;
      cosL = cos(arg);
      sinL = sin(arg);
 
      mt->mxx[0] = sfac*(sinD*cosL*sin2A - sin2D*sinL*cosA*cosA);
      mt->myy[0] = -sfac*(sinD*cosL*sin2A + sin2D*sinL*sinA*sinA);
      mt->mzz[0] = sfac*sin2D*sinL;

      mt->mxy[0] = -sfac*(sinD*cosL*cos2A + half*sin2D*sinL*sin2A);
      mt->mxz[0] = -sfac*(cosD*cosL*sinA - cos2D*sinL*cosA);
      mt->myz[0] = sfac*(cosD*cosL*cosA + cos2D*sinL*sinA);

      if(dc->strike[i]+psrc->modelrot < -999.0)  /* explosion */
         {
         mt->mxx[0] = sfac;
         mt->myy[0] = sfac;
         mt->mzz[0] = sfac;

         mt->mxy[0] = 0.0;
         mt->mxz[0] = 0.0;
         mt->myz[0] = 0.0;
         }

      tmom = dc->moment*dc->momwgts[i];
      fprintf(stderr,"            %d) Mo= %.5e dyne-cm (Mw= %.2f)\n",js,tmom,(2.0/3.0)*log10(tmom)-10.7);
      }
   else
      {
/*
   source not active in this node, but store location as negative
   of global index for later reference (e.g., SGT headers)
*/
      psrc->ix[i] = -psrc->ix[i];
      psrc->iy[i] = -psrc->iy[i];
      psrc->iz[i] = -psrc->iz[i];
      }
   }

if(js == 0)
   {
   fprintf(stderr,"            NONE\n");
   psrc->nsource = 0;
   }

fflush(stderr);
}

void init_pointmtP3(struct pntsrcs *psrc,struct runparamsP3 *rpars)
{
struct momenttensor *mt;
float tmom, hcm;
float dt, h;
int i, js;

float norm20 = 1.0e-20;
float normf = 1.0e-05; /* convert dyne-cm to dyne-km */
struct nodeinfo *ni;

ni = &(rpars->ni);
h = rpars->h;
dt = rpars->dt;

hcm = normf/h;
hcm = dt*hcm*hcm*hcm;

fprintf(stderr,"**** Point moment-tensor source\n");
fprintf(stderr,"            total sources (all nodes)= %d\n\n",psrc->nsource);
fprintf(stderr,"            sources active in this node:\n");

js = 0;
for(i=0;i<psrc->nsource;i++)
   {
/* extend to +1 for averaging (nov 2018) */
   if(psrc->ix[i] >= ni->ixminus && psrc->ix[i] <= (ni->ixplus + 1)
         && psrc->iy[i] >= ni->iyminus && psrc->iy[i] <= (ni->iyplus + 1)
               && psrc->iz[i] >= ni->izminus && psrc->iz[i] <= (ni->izplus + 1))
      {

         /*global to local conversion */
      psrc->ix[i] = psrc->ix[i] - ni->nx1;
      psrc->iy[i] = psrc->iy[i] - ni->ny1;
      psrc->iz[i] = psrc->iz[i] - ni->nz1;

      if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[i] < 1)
         psrc->iz[i] = 1;

      js++;

      mt = &(psrc->momten[i]);

      tmom = sqrt( 0.5*(norm20*mt->mxx[0])*(norm20*mt->mxx[0])
                 + 0.5*(norm20*mt->myy[0])*(norm20*mt->myy[0])
                 + 0.5*(norm20*mt->mzz[0])*(norm20*mt->mzz[0])
                     + (norm20*mt->mxy[0])*(norm20*mt->mxy[0])
                     + (norm20*mt->mxz[0])*(norm20*mt->mxz[0])
                     + (norm20*mt->myz[0])*(norm20*mt->myz[0]) );

      tmom = tmom/norm20;
      fprintf(stderr,"            %d) Mo= %.5e dyne-cm (Mw= %.2f)\n",js,tmom,(2.0/3.0)*log10(tmom)-10.7);

      mt->mxx[0] = normf*hcm*mt->mxx[0];
      mt->myy[0] = normf*hcm*mt->myy[0];
      mt->mzz[0] = normf*hcm*mt->mzz[0];
      mt->mxy[0] = normf*hcm*mt->mxy[0];
      mt->mxz[0] = normf*hcm*mt->mxz[0];
      mt->myz[0] = normf*hcm*mt->myz[0];
      }
   else
      {
/*
   source not active in this node, but store location as negative
   of global index for later reference (e.g., SGT headers)
*/
      psrc->ix[i] = -psrc->ix[i];
      psrc->iy[i] = -psrc->iy[i];
      psrc->iz[i] = -psrc->iz[i];
      }
   }

if(js == 0)
   {
   fprintf(stderr,"            NONE\n");
   psrc->nsource = 0;
   }

fflush(stderr);
}

void init_bforceP3(struct pntsrcs *psrc,struct runparamsP3 *rpars)
{
struct momenttensor *mt;
float tmom, hcm;
float dt, h;
int i, js;

float normf = 1.0e-05; /* convert dyne-cm to Nt-km */
struct nodeinfo *ni;

ni = &(rpars->ni);
h = rpars->h;
dt = rpars->dt;

hcm = normf/h;
hcm = dt*hcm*hcm*hcm;

fprintf(stderr,"**** Point body-force source\n");
fprintf(stderr,"            total sources (all nodes)= %d\n\n",psrc->nsource);
fprintf(stderr,"            sources active in this node:\n");

js = 0;
for(i=0;i<psrc->nsource;i++)
   {
   psrc->fxsrc[i] = psrc->fxsrc[i]/((float)(1.0*psrc->nsource));
   psrc->fysrc[i] = psrc->fysrc[i]/((float)(1.0*psrc->nsource));
   psrc->fzsrc[i] = psrc->fzsrc[i]/((float)(1.0*psrc->nsource));

   if(psrc->ix[i] >= ni->ixminus && psrc->ix[i] <= ni->ixplus
         && psrc->iy[i] >= ni->iyminus && psrc->iy[i] <= ni->iyplus
	       && psrc->iz[i] >= ni->izminus && psrc->iz[i] <= ni->izplus)
      {

         /*global to local conversion */
      psrc->ix[i] = psrc->ix[i] - ni->nx1;
      psrc->iy[i] = psrc->iy[i] - ni->ny1;
      psrc->iz[i] = psrc->iz[i] - ni->nz1;

      if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[i] < 1)
         psrc->iz[i] = 1;

      js++;

      fprintf(stderr,"            %d) (%d,%d,%d) Fx= %.5e dyne, Fy= %.5e dyne, Fz= %.5e dyne\n",js,psrc->ix[i]+ni->nx1,psrc->iy[i]+ni->ny1,psrc->iz[i]+ni->nz1,psrc->fxsrc[i],psrc->fysrc[i],psrc->fzsrc[i]);

      psrc->fxsrc[i] = psrc->fxsrc[i]*hcm;
      psrc->fysrc[i] = psrc->fysrc[i]*hcm;
      psrc->fzsrc[i] = psrc->fzsrc[i]*hcm;

/* 7/29/98 - RWG:

   For a body force located exactly at the free-surface, the amplitude
   of the radiated field is low by about a factor of two.  This was
   explored and checked empirically with reciprocal experiments and
   comparisons with FK runs.  I'm not sure why this happens, but I suspect
   that it has to do with some mismatch of the free-surface boundary
   condition that occurs when a source is located there.  For the zero-stress
   free-surface formulation centered on the normal stress nodes, this
   condition applies to the x and y body forces when zsrc=1.

   The soultion is to multiply these force components by two.

   See also sgtheader information.
*/

      if(rpars->freesurf && ni->minusId_z < 0 && psrc->iz[i] == 1)
         {
         psrc->fxsrc[i] = 2.0*psrc->fxsrc[i];
         psrc->fysrc[i] = 2.0*psrc->fysrc[i];
         }
      }
   else
      {
/*
   source not active in this node, but store location as negative
   of global index for later reference (e.g., SGT headers)
*/
      psrc->ix[i] = -psrc->ix[i];
      psrc->iy[i] = -psrc->iy[i];
      psrc->iz[i] = -psrc->iz[i];
      }
   }

if(js == 0)
   {
   fprintf(stderr,"            NONE\n");
   psrc->nsource = 0;
   }

fflush(stderr);
}

/* NOT SEGMENT*/
void get_srf_parP3_NOTSEGMENT(struct pntsrcs *psrc,struct runparamsP3 *rp,int bfilt,float *flo,float *fhi,float *tsh)
{
FILE *fpw, *fpr, *fopfile();
int iyleft, iyright;
float modellon, modellat, modelrot;
float *stf1, *stf2, *stf3, *space;
float lon, lat, dep, stk, dip, rake, area, tinit, dt_stf, slip1, slip2, slip3;
float de, dn, fnt, da, xx, yy;
int nt1, nt2, nt3, ntmax, ntpad, ntrsmp, ntout, resamp;
int isrc, j, jj, kk, it, nsrc;
char str[1024], pword[32], *rchar;

int gnt;
float tol = 1.0e-02;

struct momenttensor *mt;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
float latavg;
float cosR, sinR, kmlon, kmlat;
float invh;

double g0, b0;
int xy2ll = 0;
int ll2xy = 1;

int inside = 1;
float zap = 0.0;

float nyq_perc = 1.0;
float tap_perc = 0.0;

int order = 4;

int itsh;
int phase = 0;

struct nodeinfo *ni;

ni = &(rp->ni);

psrc->modelrot = rp->modelrot;

invh = 1.0/(rp->h);

cosR = cos(rp->modelrot*rperd);
sinR = sin(rp->modelrot*rperd);

itsh = 0;
if(bfilt)
   {
   itsh = 3.0/((rp->dt)*(*flo));
   *tsh = itsh*(rp->dt);
   }

fpr = fopfile(psrc->faultfile,"r");

/* 09/22/05
   For now, simply assume ASCII input and scan down to find "POINTS" line
*/

fgets(str,1024,fpr);
sscanf(str,"%s",pword);
while(strncmp(pword,"POINTS",6) != 0)
   {
   rchar = fgets(str,1024,fpr);
   if(rchar == NULL)
      {
      fprintf(stderr,"Unexpected EOF in %s, exiting...\n",psrc->faultfile);
      exit(-99);
      }
   sscanf(str,"%s",pword);
   }

sscanf(str,"%*s %d",&nsrc);

psrc->momten = (struct momenttensor *) check_malloc (nsrc*sizeof(struct momenttensor));

psrc->xs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ys = (float *) check_malloc (nsrc*sizeof(float));
psrc->zs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ix = (int *) check_malloc (nsrc*sizeof(int));
psrc->iy = (int *) check_malloc (nsrc*sizeof(int));
psrc->iz = (int *) check_malloc (nsrc*sizeof(int));
psrc->it = (int *) check_malloc (nsrc*sizeof(int));

space = NULL;
stf1 = NULL;
stf2 = NULL;
stf3 = NULL;

isrc = 0;
for(j=0;j<nsrc;j++)
   {
   if(j==nsrc || fgets(str,1024,fpr) == NULL)
      break;

   sscanf(str,"%f %f %f %f %f %f %f %f",&lon,
                                           &lat,
                                           &dep,
                                           &stk,
                                           &dip,
                                           &area,
                                           &tinit,
                                           &dt_stf);
   fgets(str,1024,fpr);
   sscanf(str,"%f %f %d %f %d %f %d",&rake,
                                        &slip1,
                                        &nt1,
                                        &slip2,
                                        &nt2,
                                        &slip3,
                                        &nt3);

   ntmax = nt1;
   if(nt2 > ntmax)
      ntmax = nt2;
   if(nt3 > ntmax)
      ntmax = nt3;

   if(ntmax)
      {
      stf1 = (float *)check_realloc((void *)stf1,ntmax*sizeof(float));
      stf2 = (float *)check_realloc((void *)stf2,ntmax*sizeof(float));
      stf3 = (float *)check_realloc((void *)stf3,ntmax*sizeof(float));
      }

   zero(stf1,ntmax);
   zero(stf2,ntmax);
   zero(stf3,ntmax);

   for(it=0;it<nt1;it++)
      fscanf(fpr,"%f",&stf1[it]);
   for(it=0;it<nt2;it++)
      fscanf(fpr,"%f",&stf2[it]);
   for(it=0;it<nt3;it++)
      fscanf(fpr,"%f",&stf3[it]);

/* get rouge newline character */
   if(nt1 || nt2 || nt3)
      fgets(str,1024,fpr);

/* convert source location to grid coordinates */

   if(rp->geoproj == 0)
      {
      de = rp->kmlon*(lon - rp->modellon);
      dn = rp->kmlat*(lat - rp->modellat);

      psrc->xs[isrc] = (de*rp->cosR - dn*rp->sinR) - rp->xshift;
      psrc->ys[isrc] = (-de*rp->sinR - dn*rp->cosR) - rp->yshift;
      }
   else if(rp->geoproj == 1)
      {
      gcproj(&xx,&yy,&lon,&lat,&rp->erad,&rp->g0,&rp->b0,rp->amat,rp->ainv,ll2xy);

      psrc->xs[isrc] = xx;
      psrc->ys[isrc] = yy;
      }

   psrc->ix[isrc] = (int)(psrc->xs[isrc]*invh + 0.5);
   psrc->iy[isrc] = (int)(psrc->ys[isrc]*invh + 0.5);

   psrc->zs[isrc] = dep;
   if(rp->freesurf)
      psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5) + 1;
   else
      psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5);

   inside = 1;
   if((psrc->ix[isrc] < 0 || psrc->ix[isrc] >= ni->globnx) ||
         (psrc->iy[isrc] < 0 || psrc->iy[isrc] >= ni->globny) ||
            (psrc->iz[isrc] < 1 || psrc->iz[isrc] >= ni->globnz))
      {
      fprintf(stderr,"**** source point %10.4f %10.4f is outside global model box\n",lon,lat);
      fprintf(stderr,"     ix= %d iy= %d iz= %d\n",psrc->ix[isrc],psrc->iy[isrc],psrc->iz[isrc]);
      inside = 0;
      }

   psrc->it[isrc] = (int)(tinit/(rp->dt) + 0.5);

/* check if point is within model sub region for this node */

   if((inside == 1)
      && psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= ni->ixplus
         && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= ni->iyplus
            && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= ni->izplus)
      {
      mt = &(psrc->momten[isrc]);

      mt->lon = lon;
      mt->lat = lat;
      mt->dep = dep;
      mt->stk = stk;
      mt->dip = dip;
      mt->rake = rake;
      mt->area = area;
      mt->tinit = tinit;
      mt->slip1 = slip1;
      mt->slip2 = slip2;
      mt->slip3 = slip3;
      mt->dt_stf = dt_stf;
      mt->nt_stf = ntmax;

      if(ntmax == 0)
	 {
         mt->nt = 0;
	 mt->stf1 = NULL;
	 mt->stf2 = NULL;
	 mt->stf3 = NULL;

	 mt->mxx = NULL;
	 mt->myy = NULL;
	 mt->mzz = NULL;
	 mt->mxy = NULL;
	 mt->mxz = NULL;
	 mt->myz = NULL;
	 }
      else if(((rp->dt) > 0.999*dt_stf && (rp->dt) < 1.001*dt_stf))
         {
	 mt->nt = ntmax + 2*itsh;
	 mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	 mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	 for(it=0;it<itsh;it++)
	    {
	    mt->stf1[it] = 0.0;
	    mt->stf2[it] = 0.0;
	    mt->stf3[it] = 0.0;
	    mt->stf1[it+ntmax+itsh] = 0.0;
	    mt->stf2[it+ntmax+itsh] = 0.0;
	    mt->stf3[it+ntmax+itsh] = 0.0;
	    }
	 for(it=0;it<ntmax;it++)
	    {
	    mt->stf1[it+itsh] = stf1[it];
	    mt->stf2[it+itsh] = stf2[it];
	    mt->stf3[it+itsh] = stf3[it];
	    }
	 }
      else  /* need to resample STFs */
	 {
	 ntpad = 2*ntmax;
	 fnt = ntpad*dt_stf/(rp->dt);
	 gnt = (int)(fnt + 0.5);

         while(nt_tol(fnt,gnt) > tol)
            {
            ntpad++;
            fnt = ntpad*dt_stf/(rp->dt);
	    gnt = (int)(fnt + 0.5);
            }

         ntrsmp = (int)(fnt);
	 ntout = (int)(ntmax*dt_stf/(rp->dt));

         if((rp->dt) < dt_stf)
            {
            resamp = 1;

            if(ntout > ntrsmp)
               ntout = ntrsmp;

            space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
            stf1 = (float *)check_realloc((void *)stf1,2*ntrsmp*sizeof(float));
            stf2 = (float *)check_realloc((void *)stf2,2*ntrsmp*sizeof(float));
            stf3 = (float *)check_realloc((void *)stf3,2*ntrsmp*sizeof(float));
            }
         else
            {
            resamp = -1;

            if(ntout > ntpad)
               ntout = ntpad;

            space = (float *) check_realloc (space,2*ntpad*sizeof(float));
            stf1 = (float *)check_realloc(stf1,2*ntpad*sizeof(float));
            stf2 = (float *)check_realloc(stf2,2*ntpad*sizeof(float));
            stf3 = (float *)check_realloc(stf3,2*ntpad*sizeof(float));
            }

         resample(stf1,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
         resample(stf2,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
         resample(stf3,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);

	 mt->nt = ntout + 2*itsh;
	 mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	 mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	 mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	 for(it=0;it<itsh;it++)
	    {
	    mt->stf1[it] = 0.0;
	    mt->stf2[it] = 0.0;
	    mt->stf3[it] = 0.0;
	    mt->stf1[it+ntout+itsh] = 0.0;
	    mt->stf2[it+ntout+itsh] = 0.0;
	    mt->stf3[it+ntout+itsh] = 0.0;
	    }
	 for(it=0;it<ntout;it++)
	    {
	    mt->stf1[it+itsh] = stf1[it];
	    mt->stf2[it+itsh] = stf2[it];
	    mt->stf3[it+itsh] = stf3[it];
	    }
	 }

      if(bfilt && mt->nt != 0)
         {
         space = (float *) check_realloc (space,4*mt->nt*sizeof(float));

         tfilter(mt->stf1,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
         tfilter(mt->stf2,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
         tfilter(mt->stf3,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);

         taper_frontback(mt->stf1,mt->nt,itsh);
         taper_frontback(mt->stf2,mt->nt,itsh);
         taper_frontback(mt->stf3,mt->nt,itsh);
         }

      isrc++;
      }
   }
fclose(fpr);

free(space);
free(stf1);
free(stf2);
free(stf3);

psrc->nsource = isrc;

if(psrc->nsource)
   {
   psrc->momten = (struct momenttensor *)check_realloc(psrc->momten,(psrc->nsource)*sizeof(struct momenttensor));

   psrc->xs = (float *)check_realloc(psrc->xs,(psrc->nsource)*sizeof(float));
   psrc->ys = (float *)check_realloc(psrc->ys,(psrc->nsource)*sizeof(float));
   psrc->zs = (float *)check_realloc(psrc->zs,(psrc->nsource)*sizeof(float));
   psrc->ix = (int *)check_realloc(psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *)check_realloc(psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *)check_realloc(psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *)check_realloc(psrc->it,(psrc->nsource)*sizeof(int));
   }
else
   {
   free(psrc->momten);
   free(psrc->xs);
   free(psrc->ys);
   free(psrc->zs);
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   }
}

void init_srfP3(struct pntsrcs *psrc,float *medf,char *outf,char *xfrmt,char *tfrmt,struct runparamsP3 *rpars)
{
FILE *fpw, *fopfile();
struct doublecouple *dc;
float *mptr, *lam, *mu, *rho;
int iy, j, isrc, iloc, nloc;
int tot_nloc, it;
float  cosS, sinS, cosD, sinD, cosL, sinL, arg;
double sum, avgs, tdbl;
float *xs, *ys, *zs, *xf, *yf, *zf;
float *slipss, *slipds, *slipnm, *slipv,  *ts;
float *l2mloc, *lamloc, *muloc;
float *l2msrc, *lamsrc, *musrc;
float *stf1, *stf2, *stf3;
float *mxx, *myy, *mzz, *mxy, *mxz, *myz;
float *Mloc_xx, *Mloc_yy, *Mloc_zz, *Mloc_xy, *Mloc_xz, *Mloc_yz;
double *Mxx_src, *Myy_src, *Mzz_src, *Mxy_src, *Mxz_src, *Myz_src;
float l2mcnv, lamcnv, mucnv, vx, vy, vz, ux, uy, uz;
float *rt, *area;
float *stk, *dip, *rak, moment;
float mfit, dd, xx, yy, zz;
struct momenttensor *mt;
double conv_fac, rt2inv, conv_dynecm;
char frmt[128];
int nx, ny, nz;

float dperr = 1.0/RPERD;

int loc_iy;
float h, dt;
struct nodeinfo *ni;

ni = &(rpars->ni);
h = rpars->h;
dt = rpars->dt;

if(psrc->MTaverage == 1)
   fprintf(stderr,"         - moment-tensor averaged to center source at normal stress node\n");
else
   fprintf(stderr,"         - moment-tensor not averaged; shear components offset by 1/2 grid\n");

fflush(stderr);

/*
   Recall, global units are:

   wavefield velocity is-
    vx,vy,vz       => cm/s

   seismic velocity and density are-
    Vp,Vs     => km/s
    density   => gm/cm^3

   Lame parameters and stress are-
    Tnn,Tij,lambda,mu => km*km*gm/(s*s*cm*cm*cm)
                      = 1e+09 Nt/(m^2)
                      = 1e+09 Pa
                      = 1 GPa
2015-01-26: Above is WRONG, should be (see main_mpi.c)

    lambda,mu => km*km*gm/(s*s*cm*cm*cm)                      
    Tnn,Tij   => [km*km*gm/(s*s*cm*cm*cm)]*(cm/km)  

    The (cm/km) term in the above should be strain, but the units don't
    cancel (displacement=cm, grid_spacing=km). I guess this is 10^-5 strain(?)

    In any case, following this through:

     Tnn,Tij   = km*cm*gm/(s*s*cm*cm*cm)
               = 1e+04 Nt/(m^2)                                
               = 1e+04 Pa                                      
               = 1e-02 MPa

   distance and time measures-
    h         => km
    dt        => s

   sources are CGS-
    moment    => dyne-cm
    force     => dyne
    slip      => cm
    slip-rate => cm/s

   The SRF slip parameters are given as:

    slip      => cm
    slip-rate => cm/s
    area      => cm*cm

   'Slip' and 'slip-rate' are consistent with source CGS units, but 'area'
   which is a distance measure, needs to be converted to km^2.

   Thus, we apply a conversion factor of 1.0e-10.
*/

/*
   Also add factors here of-
    dt for time integration
    1/h^3 for volume normalization
*/

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

rt2inv = 1.0/sqrt(2.0);
conv_fac = 1.0e-10*(double)(dt)/(1.0*(double)(h)*(double)(h)*(double)(h));
conv_dynecm = 1.0e+20*(1.0*(double)(h)*(double)(h)*(double)(h));

psrc->MTavg = (struct MT_average_weights *) check_malloc ((psrc->nsource)*sizeof(struct MT_average_weights));

xs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ys = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zs = (float *) check_malloc ((psrc->nsource)*sizeof(float));
xf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
yf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
zf = (float *) check_malloc ((psrc->nsource)*sizeof(float));
ts = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
stk = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dip = (float *) check_malloc ((psrc->nsource)*sizeof(float));
l2mloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
lamloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
muloc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
l2msrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
lamsrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
musrc = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipss = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipds = (float *) check_malloc ((psrc->nsource)*sizeof(float));
slipnm = (float *) check_malloc ((psrc->nsource)*sizeof(float));
area = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xx = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_yy = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_zz = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xy = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_xz = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mloc_yz = (float *) check_malloc ((psrc->nsource)*sizeof(float));

Mxx_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));
Myy_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));
Mzz_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));
Mxy_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));
Mxz_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));
Myz_src = (double *) check_malloc ((psrc->nsource)*sizeof(double));

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   slipss[isrc] = 0.0;
   slipds[isrc] = 0.0;
   slipnm[isrc] = 0.0;
   ts[isrc] = 1.0e+15;

   Mloc_xx[isrc] = 0.0;
   Mloc_yy[isrc] = 0.0;
   Mloc_zz[isrc] = 0.0;
   Mloc_xy[isrc] = 0.0;
   Mloc_xz[isrc] = 0.0;
   Mloc_yz[isrc] = 0.0;

   Mxx_src[isrc] = 0.0;
   Myy_src[isrc] = 0.0;
   Mzz_src[isrc] = 0.0;
   Mxy_src[isrc] = 0.0;
   Mxz_src[isrc] = 0.0;
   Myz_src[isrc] = 0.0;
   }

mfit = 0.0;
nloc = 0;

/*
   nloc = local (ie, within this CPU) # of unique source locations
   tot_nloc = global # of unique source locations as resolved with other CPUs
*/

for(iy=ni->iyminus;iy<=(ni->iyplus+1);iy++)
   {
   loc_iy = iy - ni->ny1;

   mptr = medf + loc_iy*N_MED_VARS*nx*nz;
   rho = mptr + 7*nx*nz;  /* raw density */
   lam = mptr + 11*nx*nz;  /* raw lambda */
   mu  = mptr + 12*nx*nz;  /* raw mu */

   for(isrc=0;isrc<psrc->nsource;isrc++) /* find sources in this y-plane */
      {
      if(iy == psrc->iy[isrc])
         {
         mt = &(psrc->momten[isrc]);

	 /* j must be in local coords. */
	 j = (psrc->ix[isrc]-ni->nx1) + (psrc->iz[isrc]-ni->nz1)*nx;

	 l2msrc[isrc] = lam[j] + 2.0*mu[j];
	 lamsrc[isrc] = lam[j];
	 musrc[isrc] = mu[j];

         mxx = mt->mxx;
         myy = mt->myy;
         mzz = mt->mzz;
         mxy = mt->mxy;
         mxz = mt->mxz;
         myz = mt->myz;

	 stf1 = mt->stf1;
	 stf2 = mt->stf2;
	 stf3 = mt->stf3;

         arg = ((mt->stk) - 90.0 - (psrc->modelrot))*RPERD;
         cosS = cos(arg);
         sinS = sin(arg);

         cosD = cos(mt->dip*RPERD);
         sinD = sin(mt->dip*RPERD);

         cosL = cos(mt->rake*RPERD);
         sinL = sin(mt->rake*RPERD);

	    /* unit normals */
	 vx = -sinD*sinS;
	 vy =  sinD*cosS;
	 vz = -cosD;

	    /* convert from CGS to Nt-km per unit volume */
	 l2mcnv = l2msrc[isrc]*mt->area*conv_fac;
	 lamcnv = lamsrc[isrc]*mt->area*conv_fac;
	 mucnv = musrc[isrc]*mt->area*conv_fac;

	 for(it=0;it<mt->nt;it++)
	    {
	    ux = -(stf3[it]*sinD - cosD*(stf1[it]*sinL + stf2[it]*cosL))*sinS
	               + (stf1[it]*cosL - stf2[it]*sinL)*cosS;

	    uy =  (stf3[it]*sinD - cosD*(stf1[it]*sinL + stf2[it]*cosL))*cosS
	               + (stf1[it]*cosL - stf2[it]*sinL)*sinS;

	    uz = -stf3[it]*cosD - (stf1[it]*sinL + stf2[it]*cosL)*sinD;

	    mxx[it] = l2mcnv*vx*ux + lamcnv*(vy*uy + vz*uz);
	    myy[it] = l2mcnv*vy*uy + lamcnv*(vx*ux + vz*uz);
	    mzz[it] = l2mcnv*vz*uz + lamcnv*(vx*ux + vy*uy);

	    mxy[it] = mucnv*(vx*uy + vy*ux);
	    mxz[it] = mucnv*(vx*uz + vz*ux);
	    myz[it] = mucnv*(vy*uz + vz*uy);
	    }

         free(mt->stf1);
         free(mt->stf2);
         free(mt->stf3);

/* RWG-20181128
 * Compute weights for MT components, cryptic but efficient for averaging of shear
 * components so as to center source at normal stress node
*/

         get_MTavg_coefs(psrc,isrc,ni);

/* only count statistics if within responsible volume for this CPU (nov 2018) */

         if(psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= ni->ixplus
               && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= ni->iyplus
                  && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= ni->izplus)
	    {
	    for(it=0;it<mt->nt;it++)
	       {
               Mxx_src[isrc] = Mxx_src[isrc] + (double)(mxx[it]);
               Myy_src[isrc] = Myy_src[isrc] + (double)(myy[it]);
               Mzz_src[isrc] = Mzz_src[isrc] + (double)(mzz[it]);
               Mxy_src[isrc] = Mxy_src[isrc] + (double)(mxy[it]);
               Mxz_src[isrc] = Mxz_src[isrc] + (double)(mxz[it]);
               Myz_src[isrc] = Myz_src[isrc] + (double)(myz[it]);
	       }

            xx = psrc->xs[isrc]; /* use exact location */
            yy = psrc->ys[isrc];
            zz = psrc->zs[isrc];

            for(iloc=0;iloc<nloc;iloc++)
               {
               if(xx == xs[iloc] && yy == ys[iloc] && zz == zs[iloc])
                  break;
               }

	    if(iloc == nloc)  /* new source location */
	       {
	       xs[iloc] = xx;
	       ys[iloc] = yy;
	       zs[iloc] = zz;
	       l2mloc[iloc] = l2msrc[isrc];
	       lamloc[iloc] = lamsrc[isrc];
	       muloc[iloc] = musrc[isrc];
	       nloc++;

               xf[iloc] = h*psrc->ix[isrc];
               yf[iloc] = h*psrc->iy[isrc];
               zf[iloc] = h*(psrc->iz[isrc] - 1);

               xx = xs[iloc] - xf[iloc];
               yy = ys[iloc] - yf[iloc];
               zz = zs[iloc] - zf[iloc];

               dd = sqrt(xx*xx + yy*yy + zz*zz);
               if(dd > mfit)
                  mfit = dd;
	       }

	    slipss[iloc] = slipss[iloc] + cosL*mt->slip1 - sinL*mt->slip2;
	    slipds[iloc] = slipds[iloc] + sinL*mt->slip1 + cosL*mt->slip2;
	    slipnm[iloc] = slipnm[iloc] + mt->slip3;

	    ux = -(mt->slip3*sinD - cosD*(mt->slip1*sinL + mt->slip2*cosL))*sinS
	               + (mt->slip1*cosL - mt->slip2*sinL)*cosS;

	    uy =  (mt->slip3*sinD - cosD*(mt->slip1*sinL + mt->slip2*cosL))*cosS
	               + (mt->slip1*cosL - mt->slip2*sinL)*sinS;

	    uz = -mt->slip3*cosD - (mt->slip1*sinL + mt->slip2*cosL)*sinD;

	    Mloc_xx[iloc] = Mloc_xx[iloc]
	                  + l2mloc[iloc]*vx*ux*mt->area
			  + lamloc[iloc]*(vy*uy + vz*uz)*mt->area;
	    Mloc_yy[iloc] = Mloc_yy[iloc]
	                  + l2mloc[iloc]*vy*uy*mt->area
			  + lamloc[iloc]*(vx*ux + vz*uz)*mt->area;
	    Mloc_zz[iloc] = Mloc_zz[iloc]
	                  + l2mloc[iloc]*vz*uz*mt->area
			  + lamloc[iloc]*(vx*ux + vy*uy)*mt->area;

	    Mloc_xy[iloc] = Mloc_xy[iloc] + muloc[iloc]*(vx*uy + vy*ux)*mt->area;
	    Mloc_xz[iloc] = Mloc_xz[iloc] + muloc[iloc]*(vx*uz + vz*ux)*mt->area;
	    Mloc_yz[iloc] = Mloc_yz[iloc] + muloc[iloc]*(vy*uz + vz*uy)*mt->area;

	    area[iloc] = mt->area;

	    stk[iloc] = mt->stk;
	    dip[iloc] = mt->dip;

            if((dt*psrc->it[isrc]) < ts[iloc]) /* use minimum tstart */
	       ts[iloc] = dt*psrc->it[isrc];

	    rt[iloc] = dt*mt->nt;
            }
	 }
      }
   }

if(nloc > 0)
   {
   rak = (float *) check_malloc (nloc*sizeof(float));
   slipv = (float *) check_malloc (nloc*sizeof(float));
   }

mpi_global_val(&nloc,&tot_nloc,"Nloc",1,MPI_INT,MPI_SUM);

/*
   rak = final rake value for vector slip
   slipv = total relative vector slip for each source location
   sum = total relative moment for all sources

   full moment tensor: Mo = sqrt(0.5*[MM]),
   [MM]=double dot product of moment tensor, Dahlen and Tromp
*/

sum = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   {
   if(slipss[iloc] != 0.0)
      {
      rak[iloc] = atan(slipds[iloc]/slipss[iloc]);
      if(slipss[iloc] < 0.0)
	 rak[iloc] = rak[iloc] + PI;
      }
   else if(slipds[iloc] < 0.0)
      rak[iloc] = -0.5*PI;
   else
      rak[iloc] = 0.5*PI;

   rak[iloc] = rak[iloc]*dperr;

   slipv[iloc] = sqrt(slipss[iloc]*slipss[iloc] + slipds[iloc]*slipds[iloc]);

/* purely double couple
   sum = sum + slipv[iloc]*muloc[iloc]*area[iloc]*1.0e+10;
*/

/*
   full moment tensor: Mo = sqrt(0.5*[MM]),
   [MM]=double dot product of moment tensor, Dahlen and Tromp
*/
   sum = sum + rt2inv*sqrt(Mloc_xx[iloc]*Mloc_xx[iloc]
                         + Mloc_yy[iloc]*Mloc_yy[iloc]
	                 + Mloc_zz[iloc]*Mloc_zz[iloc]
	                 + 2.0*Mloc_xy[iloc]*Mloc_xy[iloc]
	                 + 2.0*Mloc_xz[iloc]*Mloc_xz[iloc]
	                 + 2.0*Mloc_yz[iloc]*Mloc_yz[iloc])*1.0e+10;
   }

tdbl = sum;
mpi_global_val(&tdbl,&sum,"Mom(sum)",1,MPI_DOUBLE,MPI_SUM);

sum = 0.0;
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   sum = sum + rt2inv*sqrt(Mxx_src[isrc]*Mxx_src[isrc]
                         + Myy_src[isrc]*Myy_src[isrc]
	                 + Mzz_src[isrc]*Mzz_src[isrc]
	                 + 2.0*Mxy_src[isrc]*Mxy_src[isrc]
	                 + 2.0*Mxz_src[isrc]*Mxz_src[isrc]
	                 + 2.0*Myz_src[isrc]*Myz_src[isrc])*conv_dynecm;
   }

tdbl = sum;
mpi_global_val(&tdbl,&sum,"Mom(src)",1,MPI_DOUBLE,MPI_SUM);

moment = sum;

avgs = 0.0;
for(iloc=0;iloc<nloc;iloc++)
   avgs = avgs + slipv[iloc];

tdbl = avgs/(tot_nloc);
mpi_global_val(&tdbl,&avgs,"Avgs",1,MPI_DOUBLE,MPI_SUM);

fprintf(stderr,"     total moment= %13.5e (dyne-cm)\n",moment);
fprintf(stderr,"     average slip= %13.5le (cm)\n",avgs);
fflush(stderr);

if(outf[0] != '\0')
   {
   fpw = fopfile(outf,"w");

   fprintf(fpw,"%d global source locations; %d local source locations\n",tot_nloc,nloc);
   fprintf(fpw,"%10.4e average slip (cm); %10.4e total moment (dyne-cm)\n",avgs,moment);

   if(nloc)
      {
      fprintf(fpw,"%10.4e maximum misfit to FD grid\n",mfit);

      fprintf(fpw,"  x(exact)   y(exact)   z(exact)    x(grid)    y(grid)    z(grid)  strike    dip   rake  slip (cm)   time (s)  trise (s) area(km*km) mu() ix iy iz\n");
      fprintf(fpw,"------------------------------------------------------------------------------------------------------------------------\n");

      sprintf(frmt,"%%%s %%%s %%%s %%%s %%%s %%%s  %%6.1f %%6.1f %%6.1f %%%s %%%s %%%s %%13.5e %%13.5e %%5d %%5d %%5d\n",xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,xfrmt,tfrmt,tfrmt,tfrmt);

      for(iloc=0;iloc<nloc;iloc++)
         {
         fprintf(fpw,frmt,xs[iloc],
                             ys[iloc],
                             zs[iloc],
                             xf[iloc],
                             yf[iloc],
                             zf[iloc],
                             stk[iloc],
                             dip[iloc],
                             rak[iloc],
                             slipv[iloc],
                             ts[iloc],
                             rt[iloc],
                             area[iloc],
                             muloc[iloc],
                             (int)(xf[iloc]/h + 0.5),
                             (int)(yf[iloc]/h + 0.5),
                             (int)(zf[iloc]/h + 1.5));
         }
      }
   fclose(fpw);
   }

/* transform to local indexing */
for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   psrc->ix[isrc] = psrc->ix[isrc] - ni->nx1;
   psrc->iy[isrc] = psrc->iy[isrc] - ni->ny1;
   psrc->iz[isrc] = psrc->iz[isrc] - ni->nz1;

   if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[isrc] < 1)
      psrc->iz[isrc] = 1;
   }

free(xs);
free(ys);
free(zs);
free(xf);
free(yf);
free(zf);
free(ts);
free(stk);
free(dip);
free(l2mloc);
free(lamloc);
free(muloc);
free(l2msrc);
free(lamsrc);
free(musrc);
free(slipss);
free(slipds);
free(slipnm);

free(rt);
free(area);
free(Mloc_xx);
free(Mloc_yy);
free(Mloc_zz);
free(Mloc_xy);
free(Mloc_xz);
free(Mloc_yz);

free(Mxx_src);
free(Myy_src);
free(Mzz_src);
free(Mxy_src);
free(Mxz_src);
free(Myz_src);

if(nloc > 0)
   {
   free(rak);
   free(slipv);
   }
}

void get_adjsrc_parP3(struct pntsrcs *psrc,struct runparamsP3 *rp)
{
struct seisheader shead;
int fdr, *mypoints;
float *pvbuf;
int isrc, j, it, nsrc;
int i, n, nt_max, nt_in;
int n_comps, vx_pos, vy_pos, vz_pos;
off_t i_off, cur_off, n_reed;
char t_str[16];

struct nodeinfo *ni;
struct adjointsource *adj;

int inside = 1;
int n_comps_max = 9;

ni = &(rp->ni);
psrc->modelrot = rp->modelrot;

n_comps = 0;
vx_pos = -1;
vy_pos = -1;
vz_pos = -1;

fprintf(stderr,"\n      adj_tfastest= %d\n",psrc->adjoint_tfastest);

i = 0;
while(psrc->adjoint_comps[i] != '\0' && n_comps < n_comps_max)
   {
   n = 0;
   while(psrc->adjoint_comps[i+n] != ',' && psrc->adjoint_comps[i+n] != '\0')
      n++;

   strncpy(t_str,&psrc->adjoint_comps[i],n);
   t_str[n] = '\0';

   fprintf(stderr,"      adj_comp %d: %s\n",n_comps,t_str);

   if(strncmp(t_str,"vx",2) == 0 || strncmp(t_str,"Vx",2) == 0 || strncmp(t_str,"VX",2) == 0)
      vx_pos = n_comps;
   else if(strncmp(t_str,"vy",2) == 0 || strncmp(t_str,"Vy",2) == 0 || strncmp(t_str,"VY",2) == 0)
      vy_pos = n_comps;
   else if(strncmp(t_str,"vz",2) == 0 || strncmp(t_str,"Vz",2) == 0 || strncmp(t_str,"VZ",2) == 0)
      vz_pos = n_comps;

   n_comps++;

   if(psrc->adjoint_comps[i+n] == ',')
      n++;

   i = i + n;
   }

/*
fprintf(stderr,"vx_pos= %d\n",vx_pos);
fprintf(stderr,"vy_pos= %d\n",vy_pos);
fprintf(stderr,"vz_pos= %d\n",vz_pos);
*/

fflush(stderr);

fdr = opfile_ro(psrc->adjoint_file);

n_reed = reed(fdr,&nsrc,sizeof(int));

psrc->adjsrc = (struct adjointsource *) check_malloc (nsrc*sizeof(struct adjointsource));

psrc->xs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ys = (float *) check_malloc (nsrc*sizeof(float));
psrc->zs = (float *) check_malloc (nsrc*sizeof(float));
psrc->ix = (int *) check_malloc (nsrc*sizeof(int));
psrc->iy = (int *) check_malloc (nsrc*sizeof(int));
psrc->iz = (int *) check_malloc (nsrc*sizeof(int));
psrc->it = (int *) check_malloc (nsrc*sizeof(int));

mypoints = (int *) check_malloc (nsrc*sizeof(int));

nt_max = -1;
isrc = 0;
for(j=0;j<nsrc;j++)
   {
   n_reed = reed(fdr,&shead,sizeof(struct seisheader));

   psrc->ix[isrc] = shead.ix;
   psrc->iy[isrc] = shead.iy;
   psrc->iz[isrc] = shead.iz;

   if(psrc->iz[isrc] < 1)
      psrc->iz[isrc] = 1;

   inside = 1;
   if((psrc->ix[isrc] < 0 || psrc->ix[isrc] >= ni->globnx) ||
         (psrc->iy[isrc] < 0 || psrc->iy[isrc] >= ni->globny) ||
            (psrc->iz[isrc] < 1 || psrc->iz[isrc] >= ni->globnz))
      {
      fprintf(stderr,"**** source point %10.4f %10.4f is outside global model box\n",shead.slon,shead.slat);
      fprintf(stderr,"     ix= %d iy= %d iz= %d\n",psrc->ix[isrc],psrc->iy[isrc],psrc->iz[isrc]);
      inside = 0;
      }

   /*
   if((rp->dt) < 0.999*shead.dt || (rp->dt) > 1.001*shead.dt)
      {
XXXXX      NEED to Resample input    XXXXX
      }
*/

   /*
   if((rp->modelrot) < 0.999*shead.modelrot || (rp->modelrot) > 1.001*shead.modelrot)
      {
XXXXX      X & Y component orientations may be different    XXXXX
      }
*/

   psrc->it[isrc] = 0;

   psrc->xs[isrc] = (shead.h)*(psrc->ix[isrc]);
   psrc->ys[isrc] = (shead.h)*(psrc->iy[isrc]);

   if(rp->freesurf)
      psrc->zs[isrc] = (shead.h)*(psrc->iz[isrc]-1);
   else
      psrc->zs[isrc] = (shead.h)*(psrc->iz[isrc]);

/* check if point is within model sub region for this node */

   if((inside == 1)
      && psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= ni->ixplus
         && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= ni->iyplus
            && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= ni->izplus)
      {
      mypoints[isrc] = j;  /* j gives offset within data field */
      nt_in = shead.nt;

      if(nt_in > nt_max)
         nt_max = nt_in;

      isrc++;
      }
   }

/* done with headers, now compress memory to only those within this node */

psrc->nsource = isrc;

if(psrc->nsource)
   {
   psrc->adjsrc = (struct adjointsource *)check_realloc(psrc->adjsrc,(psrc->nsource)*sizeof(struct adjointsource));

   psrc->xs = (float *)check_realloc(psrc->xs,(psrc->nsource)*sizeof(float));
   psrc->ys = (float *)check_realloc(psrc->ys,(psrc->nsource)*sizeof(float));
   psrc->zs = (float *)check_realloc(psrc->zs,(psrc->nsource)*sizeof(float));
   psrc->ix = (int *)check_realloc(psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *)check_realloc(psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *)check_realloc(psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *)check_realloc(psrc->it,(psrc->nsource)*sizeof(int));

   mypoints = (int *)check_realloc(mypoints,(psrc->nsource)*sizeof(int));

   if(psrc->adjoint_tfastest == 0)   /*   (comps,srcs,time) (fastest -> slowest) */
      pvbuf = (float *) check_malloc (n_comps*nsrc*sizeof(float));
   else                              /*   (time,comps,srcs) (fastest -> slowest) */
      pvbuf = (float *) check_malloc (n_comps*nt_in*sizeof(float));

   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      adj = &(psrc->adjsrc[isrc]);

      adj->nt = rp->nt;

      adj->sx = (float *)check_malloc((adj->nt)*sizeof(float));
      adj->sy = (float *)check_malloc((adj->nt)*sizeof(float));
      adj->sz = (float *)check_malloc((adj->nt)*sizeof(float));

      for(it=0;it<adj->nt;it++)
         {
         adj->sx[it] = 0.0;
         adj->sy[it] = 0.0;
         adj->sz[it] = 0.0;
         }
      }

   if(psrc->adjoint_tfastest == 0)   /*   original way */
      {
      for(it=0;it<adj->nt;it++)
         {
         reed(fdr,pvbuf,n_comps*nsrc*sizeof(float));

         for(isrc=0;isrc<psrc->nsource;isrc++)
            {
	    i_off = mypoints[isrc]*n_comps;
            adj = &(psrc->adjsrc[isrc]);

	    if(vx_pos >= 0)
	       adj->sx[it] = pvbuf[i_off+vx_pos];
	    if(vy_pos >= 0)
	       adj->sy[it] = pvbuf[i_off+vy_pos];
	    if(vz_pos >= 0)
	       adj->sz[it] = pvbuf[i_off+vz_pos];
	    }
         }
      }
   else /* tfastest way */
      {
      cur_off = sizeof(int) + nsrc*sizeof(struct seisheader);

      for(isrc=0;isrc<psrc->nsource;isrc++)
         {
	 i_off = sizeof(int) + nsrc*sizeof(struct seisheader)
	           + mypoints[isrc]*n_comps*nt_in*sizeof(float) - cur_off;
	 lseek(fdr,i_off,SEEK_CUR);

         n_reed = reed(fdr,pvbuf,n_comps*nt_in*sizeof(float));
	 cur_off = cur_off + i_off + n_reed;

         adj = &(psrc->adjsrc[isrc]);
         for(it=0;it<adj->nt;it++)
            {
	    if(vx_pos >= 0)
	       adj->sx[it] = pvbuf[it+vx_pos*nt_in];
	    if(vy_pos >= 0)
	       adj->sy[it] = pvbuf[it+vy_pos*nt_in];
	    if(vz_pos >= 0)
	       adj->sz[it] = pvbuf[it+vz_pos*nt_in];
	    }
         }
      }

   free(pvbuf);
   }
else
   {
   free(psrc->adjsrc);
   free(psrc->xs);
   free(psrc->ys);
   free(psrc->zs);
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   }

free(mypoints);
close(fdr);
}

void init_adjsrcP3(struct pntsrcs *psrc,struct runparamsP3 *rpars)
{
int isrc, it;
float h, dt, hcm;
struct adjointsource *adj;
struct nodeinfo *ni;
float xfac, yfac, zfac;
float normf = 1.0e+05; /* convert km to cm */

ni = &(rpars->ni);
h = rpars->h;
dt = rpars->dt;

/* shouldn't need to do volume normalization
hcm = 1.0/(normf*h);
*/

hcm = dt;

/*
   Add factors of dt for time integration

  And adjust to local indexing
*/

fprintf(stderr,"\n");
fprintf(stderr,"     adjoint sources in this node= %d\n",psrc->nsource);
fflush(stderr);

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   psrc->ix[isrc] = psrc->ix[isrc] - ni->nx1;
   psrc->iy[isrc] = psrc->iy[isrc] - ni->ny1;
   psrc->iz[isrc] = psrc->iz[isrc] - ni->nz1;

   xfac = hcm;
   yfac = hcm;
   zfac = hcm;

   if(rpars->freesurf == 1 && ni->minusId_z < 0)
      {
      if(psrc->iz[isrc] < 1)
         psrc->iz[isrc] = 1;

      /*  free-surface adjustment, Not sure if this is needed?
      */
      if(psrc->iz[isrc] == 1)
         {
         xfac = 2.0*xfac;
         yfac = 2.0*yfac;
	 }
      }

   fprintf(stderr,"     source: %d\n",isrc);
   fprintf(stderr,"        ix_loc= %d\n",psrc->ix[isrc]);
   fprintf(stderr,"        iy_loc= %d\n",psrc->iy[isrc]);
   fprintf(stderr,"        iz_loc= %d\n",psrc->iz[isrc]);
   fprintf(stderr,"        it_loc= %d\n",psrc->it[isrc]);
   fflush(stderr);

   adj = &(psrc->adjsrc[isrc]);

   adj->ix = psrc->ix[isrc];
   adj->iy = psrc->iy[isrc];
   adj->iz = psrc->iz[isrc];

   for(it=0;it<adj->nt;it++)
      {
      adj->sx[it] = xfac*adj->sx[it];
      adj->sy[it] = yfac*adj->sy[it];
      adj->sz[it] = zfac*adj->sz[it];
      }
   }
}

/* SEGMENT*/
void get_srf_parP3(struct pntsrcs *psrc,struct runparamsP3 *rp,int bfilt,float *flo,float *fhi,float *tsh)
{
FILE *fpw, *fpr, *fopfile();
int iyleft, iyright;
float modellon, modellat, modelrot;
float *stf1, *stf2, *stf3, *space;
float lon, lat, dep, stk, dip, rake, area, tinit, dt_stf, slip1, slip2, slip3;
float de, dn, fnt, da, xx, yy;
int nt1, nt2, nt3, ntmax, ntpad, ntrsmp, ntout, resamp;
int isrc, j, jj, kk, it, nsrc;
char str[MAXLINE], pword[32], *rchar;

int gnt;
float tol = 1.0e-02;

struct momenttensor *mt;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc;
float latavg;
float cosR, sinR, kmlon, kmlat;
float invh;

double g0, b0;
int xy2ll = 0;
int ll2xy = 1;

int inside = 1;
float zap = 0.0;

float nyq_perc = 1.0;
float tap_perc = 0.0;

int order = 4;

int itsh;
int phase = 0;

struct nodeinfo *ni;
size_t blen, nt_tot;

ni = &(rp->ni);

psrc->modelrot = rp->modelrot;

invh = 1.0/(rp->h);

cosR = cos(rp->modelrot*rperd);
sinR = sin(rp->modelrot*rperd);

itsh = 0;
if(bfilt)
   {
   itsh = 3.0/((rp->dt)*(*flo));
   *tsh = itsh*(rp->dt);
   }

psrc->momten = NULL;

psrc->xs = NULL;
psrc->ys = NULL;
psrc->zs = NULL;
psrc->ix = NULL;
psrc->iy = NULL;
psrc->iz = NULL;
psrc->it = NULL;

space = NULL;
stf1 = NULL;
stf2 = NULL;
stf3 = NULL;

nt_tot = 0;

fpr = fopfile(psrc->faultfile,"r");

/* 09/22/05
   For now, simply assume ASCII input and scan down to find "POINTS" line
*/

fgets(str,MAXLINE,fpr);
sscanf(str,"%s",pword);
while(strncmp(pword,"POINTS",6) != 0)
   {
   rchar = fgets(str,MAXLINE,fpr);
   if(rchar == NULL)
      {
      fprintf(stderr,"Unexpected EOF in %s, exiting...\n",psrc->faultfile);
      exit(-99);
      }
   sscanf(str,"%s",pword);
   }

isrc = 0;
while(strncmp(pword,"POINTS",6) == 0) /* will handle multiple "POINTS" entries */
   {
   sscanf(str,"%*s %d",&nsrc);

   psrc->momten = (struct momenttensor *)check_realloc(psrc->momten,(isrc+nsrc)*sizeof(struct momenttensor));

   psrc->xs = (float *)check_realloc(psrc->xs,(isrc+nsrc)*sizeof(float));
   psrc->ys = (float *)check_realloc(psrc->ys,(isrc+nsrc)*sizeof(float));
   psrc->zs = (float *)check_realloc(psrc->zs,(isrc+nsrc)*sizeof(float));
   psrc->ix = (int *)check_realloc(psrc->ix,(isrc+nsrc)*sizeof(int));
   psrc->iy = (int *)check_realloc(psrc->iy,(isrc+nsrc)*sizeof(int));
   psrc->iz = (int *)check_realloc(psrc->iz,(isrc+nsrc)*sizeof(int));
   psrc->it = (int *)check_realloc(psrc->it,(isrc+nsrc)*sizeof(int));

   for(j=0;j<nsrc;j++)
      {
      if(fgets(str,MAXLINE,fpr) == NULL)
         break;

      sscanf(str,"%f %f %f %f %f %f %f %f",&lon,
                                           &lat,
                                           &dep,
                                           &stk,
                                           &dip,
                                           &area,
                                           &tinit,
                                           &dt_stf);
      fgets(str,MAXLINE,fpr);
      sscanf(str,"%f %f %d %f %d %f %d",&rake,
                                        &slip1,
                                        &nt1,
                                        &slip2,
                                        &nt2,
                                        &slip3,
                                        &nt3);

      ntmax = nt1;
      if(nt2 > ntmax)
         ntmax = nt2;
      if(nt3 > ntmax)
         ntmax = nt3;

      if(ntmax)
         {
         stf1 = (float *)check_realloc((void *)stf1,ntmax*sizeof(float));
         stf2 = (float *)check_realloc((void *)stf2,ntmax*sizeof(float));
         stf3 = (float *)check_realloc((void *)stf3,ntmax*sizeof(float));
         }

      zero(stf1,ntmax);
      zero(stf2,ntmax);
      zero(stf3,ntmax);

      for(it=0;it<nt1;it++)
         fscanf(fpr,"%f",&stf1[it]);
      for(it=0;it<nt2;it++)
         fscanf(fpr,"%f",&stf2[it]);
      for(it=0;it<nt3;it++)
         fscanf(fpr,"%f",&stf3[it]);

/* get rouge newline character */
      if(nt1 || nt2 || nt3)
         fgets(str,MAXLINE,fpr);

/* convert source location to grid coordinates */

      if(rp->geoproj == 0)
         {
         de = rp->kmlon*(lon - rp->modellon);
         dn = rp->kmlat*(lat - rp->modellat);

         psrc->xs[isrc] = (de*rp->cosR - dn*rp->sinR) - rp->xshift;
         psrc->ys[isrc] = (-de*rp->sinR - dn*rp->cosR) - rp->yshift;
         }
      else if(rp->geoproj == 1)
         {
         gcproj(&xx,&yy,&lon,&lat,&rp->erad,&rp->g0,&rp->b0,rp->amat,rp->ainv,ll2xy);

         psrc->xs[isrc] = xx;
         psrc->ys[isrc] = yy;
         }

      psrc->ix[isrc] = (int)(psrc->xs[isrc]*invh + 0.5);
      psrc->iy[isrc] = (int)(psrc->ys[isrc]*invh + 0.5);

      psrc->zs[isrc] = dep;
      if(rp->freesurf)
         psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5) + 1;
      else
         psrc->iz[isrc] = (int)(psrc->zs[isrc]*invh + 0.5);

      inside = 1;
      if((psrc->ix[isrc] < 0 || psrc->ix[isrc] >= ni->globnx) ||
            (psrc->iy[isrc] < 0 || psrc->iy[isrc] >= ni->globny) ||
               (psrc->iz[isrc] < 1 || psrc->iz[isrc] >= ni->globnz))
         {
         fprintf(stderr,"**** source point %10.4f %10.4f is outside global model box\n",lon,lat);
         fprintf(stderr,"     ix= %d iy= %d iz= %d\n",psrc->ix[isrc],psrc->iy[isrc],psrc->iz[isrc]);
         inside = 0;
         }

      psrc->it[isrc] = (int)(tinit/(rp->dt) + 0.5);

/* check if point is within model sub region for this node */

/* extend to +1 for averaging (nov 2018) */
      if((inside == 1)
         && psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= (ni->ixplus + 1)
            && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= (ni->iyplus + 1)
               && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= (ni->izplus + 1))
         {
         mt = &(psrc->momten[isrc]);

         mt->lon = lon;
         mt->lat = lat;
         mt->dep = dep;
         mt->stk = stk;
         mt->dip = dip;
         mt->rake = rake;
         mt->area = area;
         mt->tinit = tinit;
         mt->slip1 = slip1;
         mt->slip2 = slip2;
         mt->slip3 = slip3;
         mt->dt_stf = dt_stf;
         mt->nt_stf = ntmax;

         if(ntmax == 0)
	    {
            mt->nt = 0;
	    mt->stf1 = NULL;
	    mt->stf2 = NULL;
	    mt->stf3 = NULL;

	    mt->mxx = NULL;
	    mt->myy = NULL;
	    mt->mzz = NULL;
	    mt->mxy = NULL;
	    mt->mxz = NULL;
	    mt->myz = NULL;
	    }
         else if(((rp->dt) > 0.999*dt_stf && (rp->dt) < 1.001*dt_stf))
            {
	    mt->nt = ntmax + 2*itsh;
	    mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	    mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	    for(it=0;it<itsh;it++)
	       {
	       mt->stf1[it] = 0.0;
	       mt->stf2[it] = 0.0;
	       mt->stf3[it] = 0.0;
	       mt->stf1[it+ntmax+itsh] = 0.0;
	       mt->stf2[it+ntmax+itsh] = 0.0;
	       mt->stf3[it+ntmax+itsh] = 0.0;
	       }
	    for(it=0;it<ntmax;it++)
	       {
	       mt->stf1[it+itsh] = stf1[it];
	       mt->stf2[it+itsh] = stf2[it];
	       mt->stf3[it+itsh] = stf3[it];
	       }
	    }
         else  /* need to resample STFs */
	    {
	    ntpad = 2*ntmax;
	    fnt = ntpad*dt_stf/(rp->dt);
	    gnt = (int)(fnt + 0.5);

            while(nt_tol(fnt,gnt) > tol)
               {
               ntpad++;
               fnt = ntpad*dt_stf/(rp->dt);
	       gnt = (int)(fnt + 0.5);
               }

            ntrsmp = (int)(fnt);
	    ntout = (int)(ntmax*dt_stf/(rp->dt));

            if((rp->dt) < dt_stf)
               {
               resamp = 1;

               if(ntout > ntrsmp)
                  ntout = ntrsmp;

               space = (float *) check_realloc ((void *)space,2*ntrsmp*sizeof(float));
               stf1 = (float *)check_realloc((void *)stf1,2*ntrsmp*sizeof(float));
               stf2 = (float *)check_realloc((void *)stf2,2*ntrsmp*sizeof(float));
               stf3 = (float *)check_realloc((void *)stf3,2*ntrsmp*sizeof(float));
               }
            else
               {
               resamp = -1;

               if(ntout > ntpad)
                  ntout = ntpad;

               space = (float *) check_realloc (space,2*ntpad*sizeof(float));
               stf1 = (float *)check_realloc(stf1,2*ntpad*sizeof(float));
               stf2 = (float *)check_realloc(stf2,2*ntpad*sizeof(float));
               stf3 = (float *)check_realloc(stf3,2*ntpad*sizeof(float));
               }

            resample(stf1,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
            resample(stf2,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);
            resample(stf3,ntmax,&dt_stf,resamp,ntpad,ntrsmp,&(rp->dt),space,order,&nyq_perc,&tap_perc);

	    mt->nt = ntout + 2*itsh;
	    mt->stf1 = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->stf2 = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->stf3 = (float *)check_malloc((mt->nt)*sizeof(float));

	    mt->mxx = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->myy = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mzz = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mxy = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->mxz = (float *)check_malloc((mt->nt)*sizeof(float));
	    mt->myz = (float *)check_malloc((mt->nt)*sizeof(float));

	    for(it=0;it<itsh;it++)
	       {
	       mt->stf1[it] = 0.0;
	       mt->stf2[it] = 0.0;
	       mt->stf3[it] = 0.0;
	       mt->stf1[it+ntout+itsh] = 0.0;
	       mt->stf2[it+ntout+itsh] = 0.0;
	       mt->stf3[it+ntout+itsh] = 0.0;
	       }
	    for(it=0;it<ntout;it++)
	       {
	       mt->stf1[it+itsh] = stf1[it];
	       mt->stf2[it+itsh] = stf2[it];
	       mt->stf3[it+itsh] = stf3[it];
	       }
	    }

	 nt_tot = nt_tot + mt->nt;

         if(bfilt && mt->nt != 0)
            {
            space = (float *) check_realloc (space,4*mt->nt*sizeof(float));

            tfilter(mt->stf1,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
            tfilter(mt->stf2,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);
            tfilter(mt->stf3,&(rp->dt),mt->nt,bfilt,fhi,flo,phase,space);

            taper_frontback(mt->stf1,mt->nt,itsh);
            taper_frontback(mt->stf2,mt->nt,itsh);
            taper_frontback(mt->stf3,mt->nt,itsh);
            }

         isrc++;
         }
      }

   if(fgets(str,MAXLINE,fpr) == NULL)
      break;
   else
      sscanf(str,"%s",pword);
   }

fclose(fpr);

free(space);
free(stf1);
free(stf2);
free(stf3);

psrc->nsource = isrc;

if(psrc->nsource)
   {
   psrc->momten = (struct momenttensor *)check_realloc(psrc->momten,(psrc->nsource)*sizeof(struct momenttensor));

   psrc->xs = (float *)check_realloc(psrc->xs,(psrc->nsource)*sizeof(float));
   psrc->ys = (float *)check_realloc(psrc->ys,(psrc->nsource)*sizeof(float));
   psrc->zs = (float *)check_realloc(psrc->zs,(psrc->nsource)*sizeof(float));
   psrc->ix = (int *)check_realloc(psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *)check_realloc(psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *)check_realloc(psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *)check_realloc(psrc->it,(psrc->nsource)*sizeof(int));

   blen = psrc->nsource*(7*sizeof(float) + sizeof(struct momenttensor)) + 9*nt_tot*sizeof(float);

   fprintf(stderr,"\n");
   fprintf(stderr,"**** Memory stats for finite-fault source (approximate):\n");
   fprintf(stderr,"          total for source= %8.1f Mb\n\n",(float)(blen/1.0e+06));
   fflush(stderr);
   }
else
   {
   free(psrc->momten);
   free(psrc->xs);
   free(psrc->ys);
   free(psrc->zs);
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   }
}

void get_dc_parP3_ns(struct pntsrcs *psrc,struct runparamsP3 *rp,float *tz,float *dt)
{
struct doublecouple *dc;
float sum, invdt;
int isrc, i;

struct nodeinfo *ni;
int *ixs_in, *iys_in, *izs_in, *its_in;
float *tsrc_in, *rtime_in, *stk_in, *dip_in, *rak_in, *mom_in;

ni = &(rp->ni);
dc = &psrc->doublec;

ixs_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));
iys_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));
izs_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
its_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));

stk_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dip_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rak_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
mom_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   rtime_in[isrc] = *tz;
 
mstpar("xsrc","vd",ixs_in);
mstpar("ysrc","vd",iys_in);
mstpar("zsrc","vd",izs_in);
mstpar("strike","vf",stk_in);
mstpar("dip","vf",dip_in);
mstpar("rake","vf",rak_in);
getpar("rtime","vf",rtime_in);

mstpar("moment","f",&dc->moment);

getpar("modelrot","f",&psrc->modelrot);
getpar("relative_slip","d",&psrc->relative_slip);
getpar("absolute_slip","d",&psrc->absolute_slip);

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc_in);
   mstpar("momwgts","vf",mom_in);

   invdt = 1.0/(*dt);
   sum = 0.0;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      its_in[isrc] = tsrc_in[isrc]*invdt;
      sum = sum + mom_in[isrc];
      }

   sum = 1.0/sum;
   for(isrc=0;isrc<psrc->nsource;isrc++)
      mom_in[isrc] = sum*mom_in[isrc];
   }
else
   {
   its_in[0] = 0;
   mom_in[0] = 1.0;
   }

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->it = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->rtime = (float *) check_malloc ((psrc->nsource)*sizeof(float));

dc->strike  = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->dip     = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->rake    = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->momwgts = (float *) check_malloc ((psrc->nsource)*sizeof(float));
dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

isrc = 0;
for(i=0;i<psrc->nsource;i++)
   {
/* extend to +1 for averaging (nov 2018) */
   if(ixs_in[i] >= ni->ixminus && ixs_in[i] <= (ni->ixplus + 1)
         && iys_in[i] >= ni->iyminus && iys_in[i] <= (ni->iyplus + 1)
	    && izs_in[i] >= ni->izminus && izs_in[i] <= (ni->izplus + 1))
      {
      psrc->ix[isrc] = ixs_in[i];
      psrc->iy[isrc] = iys_in[i];
      psrc->iz[isrc] = izs_in[i];
      psrc->it[isrc] = its_in[i];
      psrc->stfindx[isrc] = isrc;
      psrc->rtime[isrc] = rtime_in[i];

      dc->strike[isrc] = stk_in[i];
      dc->dip[isrc] = dip_in[i];
      dc->rake[isrc] = rak_in[i];
      dc->momwgts[isrc] = mom_in[i];
      dc->itsrc[isrc] = its_in[i];

      isrc++;
      }
   }

psrc->nsource = isrc;
psrc->nstf = isrc;

if(psrc->nsource)
   {
   psrc->ix = (int *) check_realloc (psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *) check_realloc (psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *) check_realloc (psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *) check_realloc (psrc->it,(psrc->nsource)*sizeof(int));
   psrc->stfindx = (int *) check_realloc (psrc->stfindx,(psrc->nsource)*sizeof(int));
   psrc->rtime = (float *) check_realloc (psrc->rtime,(psrc->nsource)*sizeof(float));

   dc->strike  = (float *) check_realloc (dc->strike,(psrc->nsource)*sizeof(float));
   dc->dip     = (float *) check_realloc (dc->dip,(psrc->nsource)*sizeof(float));
   dc->rake    = (float *) check_realloc (dc->rake,(psrc->nsource)*sizeof(float));
   dc->momwgts = (float *) check_realloc (dc->momwgts,(psrc->nsource)*sizeof(float));
   dc->itsrc = (int *) check_realloc (dc->itsrc,(psrc->nsource)*sizeof(int));
   }
else
   {
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   free(psrc->stfindx);
   free(psrc->rtime);

   free(dc->strike);
   free(dc->dip);
   free(dc->rake);
   free(dc->momwgts);
   free(dc->itsrc);
   }

free(ixs_in);
free(iys_in);
free(izs_in);

free(tsrc_in);
free(rtime_in);
free(its_in);

free(stk_in);
free(dip_in);
free(rak_in);
free(mom_in);
}

void init_dcP3_stf(struct pntsrcs *psrc,struct runparamsP3 *rpars,float *stfunc,int nt)
{
struct doublecouple *dc;
struct momenttensor *mt;
float cosD, sinD, cosL, sinL, cosA, sinA;
float cos2D, sin2D, cos2A, sin2A;
float sfac, arg, hcm, invh4, tmom;
float dt, h, *stfp;
int i, it;

float half = 0.5;
float two = 2.0;
float normf = 1.0e-05; /* convert dyne-cm to Nt-km */
struct nodeinfo *ni;

ni = &(rpars->ni);
dc = &psrc->doublec;
h = rpars->h;
dt = rpars->dt;

if(psrc->MTaverage == 1)
   fprintf(stderr,"         - moment-tensor averaged to center source at normal stress node\n");
else
   fprintf(stderr,"         - moment-tensor not averaged; shear components offset by 1/2 grid\n");

fflush(stderr);

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));
psrc->MTavg = (struct MT_average_weights *) check_malloc ((psrc->nsource)*sizeof(struct MT_average_weights));

hcm = normf/h;
hcm = dt*hcm*hcm*hcm;

invh4 = hcm*(dc->moment)*(normf);

fprintf(stderr,"**** Point double-couple source\n");
fprintf(stderr,"            sources active in this node:\n");

for(i=0;i<psrc->nsource;i++)
   {
   dc->strike[i] = dc->strike[i] - psrc->modelrot; /* global to local conversion */
   mt = &(psrc->momten[i]);

   mt->nt = nt;
   mt->mxx = (float *) check_malloc (nt*sizeof(float));
   mt->myy = (float *) check_malloc (nt*sizeof(float));
   mt->mzz = (float *) check_malloc (nt*sizeof(float));
   mt->mxy = (float *) check_malloc (nt*sizeof(float));
   mt->mxz = (float *) check_malloc (nt*sizeof(float));
   mt->myz = (float *) check_malloc (nt*sizeof(float));

   sfac = invh4*dc->momwgts[i];

   arg = dc->strike[i]*RPERD;
   cosA = cos(arg);
   sinA = sin(arg);

   cos2A = cosA*cosA - sinA*sinA;
   sin2A = two*sinA*cosA;

   arg = dc->dip[i]*RPERD;
   cosD = cos(arg);
   sinD = sin(arg);

   cos2D = cosD*cosD - sinD*sinD;
   sin2D = two*sinD*cosD;

   arg = dc->rake[i]*RPERD;
   cosL = cos(arg);
   sinL = sin(arg);

   stfp = stfunc + nt*psrc->stfindx[i];

   if(dc->strike[i]+psrc->modelrot < -999.0)  /* explosion */
      {
      for(it=0;it<nt;it++)
         {
         mt->mxx[it] = sfac*stfp[it];
         mt->myy[it] = sfac*stfp[it];
         mt->mzz[it] = sfac*stfp[it];

         mt->mxy[it] = 0.0;
         mt->mxz[it] = 0.0;
         mt->myz[it] = 0.0;
         }
      }
   else
      {
      for(it=0;it<nt;it++)
         {
         mt->mxx[it] = sfac*(sinD*cosL*sin2A - sin2D*sinL*cosA*cosA)*stfp[it];
         mt->myy[it] = -sfac*(sinD*cosL*sin2A + sin2D*sinL*sinA*sinA)*stfp[it];
         mt->mzz[it] = sfac*sin2D*sinL*stfp[it];

         mt->mxy[it] = -sfac*(sinD*cosL*cos2A + half*sin2D*sinL*sin2A)*stfp[it];
         mt->mxz[it] = -sfac*(cosD*cosL*sinA - cos2D*sinL*cosA)*stfp[it];
         mt->myz[it] = sfac*(cosD*cosL*cosA + cos2D*sinL*sinA)*stfp[it];
         }
      }

/* only print statistics if completely within this node */
   if(psrc->ix[i] >= ni->ixminus && psrc->ix[i] <= ni->ixplus
         && psrc->iy[i] >= ni->iyminus && psrc->iy[i] <= ni->iyplus
	    && psrc->iz[i] >= ni->izminus && psrc->iz[i] <= ni->izplus)
      {
      tmom = dc->moment*dc->momwgts[i];
      fprintf(stderr,"            %d) Mo= %.5e dyne-cm (Mw= %.2f)\n",i,tmom,(2.0/3.0)*log10(tmom)-10.7);
      }

/* RWG-20181128
 * Compute weights for MT components, cryptic but efficient for averaging of shear 
 * components so as to center source at normal stress node 
*/

   get_MTavg_coefs(psrc,i,ni);
   }

/* transform to local indexing */
for(i=0;i<psrc->nsource;i++)
   {
   psrc->ix[i] = psrc->ix[i] - ni->nx1;
   psrc->iy[i] = psrc->iy[i] - ni->ny1;
   psrc->iz[i] = psrc->iz[i] - ni->nz1;

   if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[i] < 1)
      psrc->iz[i] = 1;
   }

if(psrc->nsource == 0)
   fprintf(stderr,"            NONE\n");

fflush(stderr);
}

void get_pointmt_parP3_ns(struct pntsrcs *psrc,struct runparamsP3 *rp,float *tz,float *dt)
{
struct doublecouple *dc;
struct momenttensor *mt;
float sum, invdt;
int isrc, i, nt;

float *Mrr, *Mtt, *Mpp, *Mrt, *Mrp, *Mtp;
float *Mnn, *Mee, *Mdd, *Mne, *Mnd, *Med;
float arg, cosA, sinA, cos2A, sin2A;

float rperd = RPERD;

struct nodeinfo *ni;
int *ixs_in, *iys_in, *izs_in, *its_in;
float *tsrc_in, *rtime_in;

ni = &(rp->ni);
dc = &psrc->doublec;
nt = rp->nt;

ixs_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));
iys_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));
izs_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));

tsrc_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
rtime_in = (float *) check_malloc ((psrc->nsource)*sizeof(float));
its_in = (int *) check_malloc ((psrc->nsource)*sizeof(int));

Mrr = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mpp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrt = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mrp = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mtp = (float *) check_malloc ((psrc->nsource)*sizeof(float));

Mnn = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mee = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mdd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mne = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Mnd = (float *) check_malloc ((psrc->nsource)*sizeof(float));
Med = (float *) check_malloc ((psrc->nsource)*sizeof(float));

for(isrc=0;isrc<psrc->nsource;isrc++)
   rtime_in[isrc] = *tz;

for(isrc=0;isrc<psrc->nsource;isrc++)
   {
   Mrr[isrc] = 1.0e-15;
   Mtt[isrc] = 1.0e-15;
   Mpp[isrc] = 1.0e-15;
   Mrt[isrc] = 1.0e-15;
   Mrp[isrc] = 1.0e-15;
   Mtp[isrc] = 1.0e-15;

   Mnn[isrc] = 1.0e-15;
   Mee[isrc] = 1.0e-15;
   Mdd[isrc] = 1.0e-15;
   Mne[isrc] = 1.0e-15;
   Mnd[isrc] = 1.0e-15;
   Med[isrc] = 1.0e-15;
   }

mstpar("xsrc","vd",ixs_in);
mstpar("ysrc","vd",iys_in);
mstpar("zsrc","vd",izs_in);
getpar("rtime","vf",rtime_in);
getpar("modelrot","f",&psrc->modelrot);

getpar("Mrr","vf",Mrr);
getpar("Mtt","vf",Mtt);
getpar("Mpp","vf",Mpp);
getpar("Mrt","vf",Mrt);
getpar("Mrp","vf",Mrp);
getpar("Mtp","vf",Mtp);

getpar("Mnn","vf",Mnn);
getpar("Mee","vf",Mee);
getpar("Mdd","vf",Mdd);
getpar("Mne","vf",Mne);
getpar("Mnd","vf",Mnd);
getpar("Med","vf",Med);

if(psrc->nsource > 1)
   {
   mstpar("tsrc","vf",tsrc_in);

   invdt = 1.0/(*dt);
   for(isrc=0;isrc<psrc->nsource;isrc++)
      its_in[isrc] = tsrc_in[isrc]*invdt;
   }
else
   its_in[0] = 0;

if(Mrr[0] > 1.0e-14 || -Mrr[0] > 1.0e-14)
   {
   for(isrc=0;isrc<psrc->nsource;isrc++)
      {
      Mnn[isrc] = Mtt[isrc];
      Mee[isrc] = Mpp[isrc];
      Mdd[isrc] = Mrr[isrc];
      Mne[isrc] = -Mtp[isrc];
      Mnd[isrc] = Mrt[isrc];
      Med[isrc] = -Mrp[isrc];
      }
   }

psrc->ix = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iy = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->iz = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->it = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->stfindx = (int *) check_malloc ((psrc->nsource)*sizeof(int));
psrc->rtime = (float *) check_malloc ((psrc->nsource)*sizeof(float));

psrc->momten = (struct momenttensor *) check_malloc ((psrc->nsource)*sizeof(struct momenttensor));

dc->itsrc = (int *) check_malloc ((psrc->nsource)*sizeof(int));

arg = (90.0 + psrc->modelrot)*rperd;
cosA = cos(arg);
sinA = sin(arg);

cos2A = cosA*cosA - sinA*sinA;
sin2A = 2.0*sinA*cosA;

isrc = 0;
for(i=0;i<psrc->nsource;i++)
   {
/* extend to +1 for averaging (nov 2018) */
   if(ixs_in[i] >= ni->ixminus && ixs_in[i] <= (ni->ixplus + 1)
         && iys_in[i] >= ni->iyminus && iys_in[i] <= (ni->iyplus + 1)
            && izs_in[i] >= ni->izminus && izs_in[i] <= (ni->izplus + 1))
      {
      psrc->ix[isrc] = ixs_in[i];
      psrc->iy[isrc] = iys_in[i];
      psrc->iz[isrc] = izs_in[i];
      psrc->it[isrc] = its_in[i];
      psrc->stfindx[isrc] = isrc;
      psrc->rtime[isrc] = rtime_in[i];

      dc->itsrc[isrc] = its_in[i];

      mt = &(psrc->momten[isrc]);

      mt->nt = nt;
      mt->mxx = (float *) check_malloc (nt*sizeof(float));
      mt->myy = (float *) check_malloc (nt*sizeof(float));
      mt->mzz = (float *) check_malloc (nt*sizeof(float));
      mt->mxy = (float *) check_malloc (nt*sizeof(float));
      mt->mxz = (float *) check_malloc (nt*sizeof(float));
      mt->myz = (float *) check_malloc (nt*sizeof(float));

/* 
 * store moment values in first element for now, later call to init_pointmtP3_stf()
 * will add in source time function
 */
      mt->mxx[0] = Mnn[i]*cosA*cosA + Mee[i]*sinA*sinA + Mne[i]*sin2A;
      mt->myy[0] = Mnn[i]*sinA*sinA + Mee[i]*cosA*cosA - Mne[i]*sin2A;
      mt->mzz[0] = Mdd[i];
      mt->mxy[0] = 0.5*(Mee[i] - Mnn[i])*sin2A + Mne[i]*cos2A;
      mt->mxz[0] = Mnd[i]*cosA + Med[i]*sinA;
      mt->myz[0] = -Mnd[i]*sinA + Med[i]*cosA;

      isrc++;
      }
   }

psrc->nsource = isrc;
psrc->nstf = isrc;
 
if(psrc->nsource)
   {
   psrc->ix = (int *) check_realloc (psrc->ix,(psrc->nsource)*sizeof(int));
   psrc->iy = (int *) check_realloc (psrc->iy,(psrc->nsource)*sizeof(int));
   psrc->iz = (int *) check_realloc (psrc->iz,(psrc->nsource)*sizeof(int));
   psrc->it = (int *) check_realloc (psrc->it,(psrc->nsource)*sizeof(int));
   psrc->stfindx = (int *) check_realloc (psrc->stfindx,(psrc->nsource)*sizeof(int));
   psrc->rtime = (float *) check_realloc (psrc->rtime,(psrc->nsource)*sizeof(float));

   psrc->momten = (struct momenttensor *) check_realloc (psrc->momten,(psrc->nsource)*sizeof(struct momenttensor));

   dc->itsrc = (int *) check_realloc (dc->itsrc,(psrc->nsource)*sizeof(int));
   }
else
   {
   free(psrc->ix);
   free(psrc->iy);
   free(psrc->iz);
   free(psrc->it);
   free(psrc->stfindx);
   free(psrc->rtime);

   free(psrc->momten);

   free(dc->itsrc);
   }

free(ixs_in);
free(iys_in);
free(izs_in);

free(tsrc_in);
free(rtime_in);
free(its_in);

free(Mrr);
free(Mtt);
free(Mpp);
free(Mrt);
free(Mrp);
free(Mtp);

free(Mnn);
free(Mee);
free(Mdd);
free(Mne);
free(Mnd);
free(Med);
}

void init_pointmtP3_stf(struct pntsrcs *psrc,struct runparamsP3 *rpars,float *stfunc,int nt)
{
struct momenttensor *mt;
float tmom, hcm;
float dt, h, *stfp;
float mxx_fac, myy_fac, mzz_fac, mxy_fac, mxz_fac, myz_fac;
int i, it;

float norm20 = 1.0e-20;
float normf = 1.0e-05; /* convert dyne-cm to dyne-km */
struct nodeinfo *ni;

ni = &(rpars->ni);
h = rpars->h;
dt = rpars->dt;

if(psrc->MTaverage == 1)
   fprintf(stderr,"         - moment-tensor averaged to center source at normal stress node\n");
else
   fprintf(stderr,"         - moment-tensor not averaged; shear components offset by 1/2 grid\n");

fflush(stderr);

psrc->MTavg = (struct MT_average_weights *) check_malloc ((psrc->nsource)*sizeof(struct MT_average_weights));

hcm = normf/h;
hcm = dt*hcm*hcm*hcm;

fprintf(stderr,"**** Point moment-tensor source\n");
fprintf(stderr,"            sources active in this node:\n");

for(i=0;i<psrc->nsource;i++)
   {
   mt = &(psrc->momten[i]);


/* only print statistics if completely within this node */
   if(psrc->ix[i] >= ni->ixminus && psrc->ix[i] <= ni->ixplus
         && psrc->iy[i] >= ni->iyminus && psrc->iy[i] <= ni->iyplus
	    && psrc->iz[i] >= ni->izminus && psrc->iz[i] <= ni->izplus)
      {
      tmom = sqrt( 0.5*(norm20*mt->mxx[0])*(norm20*mt->mxx[0])
                 + 0.5*(norm20*mt->myy[0])*(norm20*mt->myy[0])
                 + 0.5*(norm20*mt->mzz[0])*(norm20*mt->mzz[0])
                     + (norm20*mt->mxy[0])*(norm20*mt->mxy[0])
                     + (norm20*mt->mxz[0])*(norm20*mt->mxz[0])
                     + (norm20*mt->myz[0])*(norm20*mt->myz[0]) );

      tmom = tmom/norm20;
      fprintf(stderr,"            %d) Mo= %.5e dyne-cm (Mw= %.2f)\n",i,tmom,(2.0/3.0)*log10(tmom)-10.7);
      }

   mxx_fac = normf*hcm*mt->mxx[0];
   myy_fac = normf*hcm*mt->myy[0];
   mzz_fac = normf*hcm*mt->mzz[0];
   mxy_fac = normf*hcm*mt->mxy[0];
   mxz_fac = normf*hcm*mt->mxz[0];
   myz_fac = normf*hcm*mt->myz[0];

   stfp = stfunc + nt*psrc->stfindx[i];
   for(it=0;it<nt;it++)
      {
      mt->mxx[it] = mxx_fac*stfp[it];
      mt->myy[it] = myy_fac*stfp[it];
      mt->mzz[it] = mzz_fac*stfp[it];

      mt->mxy[it] = mxy_fac*stfp[it];
      mt->mxz[it] = mxz_fac*stfp[it];
      mt->myz[it] = myz_fac*stfp[it];

      }

/* RWG-20181128
 * Compute weights for MT components, cryptic but efficient for averaging of shear
 * components so as to center source at normal stress node
*/

   get_MTavg_coefs(psrc,i,ni);
   }

/* transform to local indexing */
for(i=0;i<psrc->nsource;i++)
   {
   psrc->ix[i] = psrc->ix[i] - ni->nx1;
   psrc->iy[i] = psrc->iy[i] - ni->ny1;
   psrc->iz[i] = psrc->iz[i] - ni->nz1;

   if(rpars->freesurf == 1 && ni->minusId_z < 0 && psrc->iz[i] < 1)
      psrc->iz[i] = 1;
   }

if(psrc->nsource == 0)
   fprintf(stderr,"            NONE\n");

fflush(stderr);
}

void get_MTavg_coefs(struct pntsrcs *psrc,int isrc,struct nodeinfo *ni)
{

/* first zero everything just to make sure (may not be necessary?) */

psrc->MTavg[isrc].wxx_x0y0z0 = 0.0;
psrc->MTavg[isrc].wyy_x0y0z0 = 0.0;
psrc->MTavg[isrc].wzz_x0y0z0 = 0.0;

psrc->MTavg[isrc].wxy_x0y0z0 = 0.0;
psrc->MTavg[isrc].wxy_xmy0z0 = 0.0;
psrc->MTavg[isrc].wxy_x0ymz0 = 0.0;
psrc->MTavg[isrc].wxy_xmymz0 = 0.0;

psrc->MTavg[isrc].wxz_x0y0z0 = 0.0;
psrc->MTavg[isrc].wxz_xmy0z0 = 0.0;
psrc->MTavg[isrc].wxz_x0y0zm = 0.0;
psrc->MTavg[isrc].wxz_xmy0zm = 0.0;

psrc->MTavg[isrc].wyz_x0y0z0 = 0.0;
psrc->MTavg[isrc].wyz_x0ymz0 = 0.0;
psrc->MTavg[isrc].wyz_x0y0zm = 0.0;
psrc->MTavg[isrc].wyz_x0ymzm = 0.0;

/* formulation below does no averaging */

if(psrc->MTaverage == 0)
   {
   if(psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= ni->ixplus
               && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= ni->iyplus
                  && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= ni->izplus)
      {
      psrc->MTavg[isrc].wxx_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wyy_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wzz_x0y0z0 = 1.0;

      psrc->MTavg[isrc].wxy_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wxy_xmy0z0 = 0.0;
      psrc->MTavg[isrc].wxy_x0ymz0 = 0.0;
      psrc->MTavg[isrc].wxy_xmymz0 = 0.0;

      psrc->MTavg[isrc].wxz_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wxz_xmy0z0 = 0.0;
      psrc->MTavg[isrc].wxz_x0y0zm = 0.0;
      psrc->MTavg[isrc].wxz_xmy0zm = 0.0;

      psrc->MTavg[isrc].wyz_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wyz_x0ymz0 = 0.0;
      psrc->MTavg[isrc].wyz_x0y0zm = 0.0;
      psrc->MTavg[isrc].wyz_x0ymzm = 0.0;
      }
   }

/* formulation below averages across adjacent nodes to center source at normal stress */

if(psrc->MTaverage == 1)
   {
   if(psrc->ix[isrc] >= ni->ixminus && psrc->ix[isrc] <= (ni->ixplus + 1)
               && psrc->iy[isrc] >= ni->iyminus && psrc->iy[isrc] <= (ni->iyplus + 1)
                  && psrc->iz[isrc] >= ni->izminus && psrc->iz[isrc] <= (ni->izplus + 1))
      {
      psrc->MTavg[isrc].wxx_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wyy_x0y0z0 = 1.0;
      psrc->MTavg[isrc].wzz_x0y0z0 = 1.0;

      psrc->MTavg[isrc].wxy_x0y0z0 = 0.25;
      psrc->MTavg[isrc].wxy_xmy0z0 = 0.25;
      psrc->MTavg[isrc].wxy_x0ymz0 = 0.25;
      psrc->MTavg[isrc].wxy_xmymz0 = 0.25;

      psrc->MTavg[isrc].wxz_x0y0z0 = 0.25;
      psrc->MTavg[isrc].wxz_xmy0z0 = 0.25;
      psrc->MTavg[isrc].wxz_x0y0zm = 0.25;
      psrc->MTavg[isrc].wxz_xmy0zm = 0.25;

      psrc->MTavg[isrc].wyz_x0y0z0 = 0.25;
      psrc->MTavg[isrc].wyz_x0ymz0 = 0.25;
      psrc->MTavg[isrc].wyz_x0y0zm = 0.25;
      psrc->MTavg[isrc].wyz_x0ymzm = 0.25;
      }
   }
}
