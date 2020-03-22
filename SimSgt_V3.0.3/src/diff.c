/*
   diff.c contains the following functions:

      setcoefs()
      diff()
      ord2()
      ord4()
      diffx()
      ord2x()
      ord4x()
      set_vminvmax_reedmod()
*/

#include "include.h"

setcoefs(order,coefs,dt,h,vsmin)
struct fdcoefs *coefs;
float *dt, *h, *vsmin;
int order;
{
float sum, f, con, vmn, vmx;
float lmin, fmax;

coefs->order = order;
coefs->ordx = ORDERX;      /* spatial differencing order below iz=izord2 */

sum = 0.0;
vmn = coefs->vmin;
vmx = coefs->vmax;
con = (*dt)/(*h);

coefs->dtoh = con;
coefs->c0 = 9.0/8.0;
coefs->c1 = -1.0/24.0;

if(coefs->order == 2)
   sum = 1.0;
else if(coefs->order == 4)
   sum = coefs->c0 - coefs->c1;

f = (*h)/(vmx*sum*sqrt(3.0));
if((*dt) > f)
   {
   fprintf(stderr,"**** Unstable run:\n");
   fprintf(stderr,"      vmax = %f\n",vmx);
   fprintf(stderr,"        dt = %f, must be less than %f\n",(*dt),f);
   fprintf(stderr,"        exiting...\n");
   exit(-1);
   }
else
   {
   fprintf(stderr,"**** Stable run:\n");
   fprintf(stderr,"      vmax = %f\n",vmx);
   fprintf(stderr,"        dt = %f, must be less than %f\n\n",(*dt),f);
   }

coefs->c0 = coefs->dtoh*coefs->c0;
coefs->c1 = coefs->dtoh*coefs->c1;

if(coefs->order == 2)
   lmin = 10*(*h);
else if(coefs->order == 4)
   lmin = 5*(*h);
else
   {
   fprintf(stderr,"**** Invalid order:\n");
   fprintf(stderr,"        order must be either 2 or 4\n");
   fprintf(stderr,"        exiting...\n");
   exit(-1);
   }

fmax = vmn/lmin;

fprintf(stderr,"**** Differencing accuracy:\n");
fprintf(stderr,"                    order= %d\n",coefs->order);
fprintf(stderr,"                     fmax= %f Hz\n",fmax);
fprintf(stderr,"                     vmin= %f km/s\n",vmn);
fprintf(stderr,"                   izord2= %d (vsmin= %f)\n\n",coefs->izord2,(*vsmin));
}

diff(n,a,x,b,c,d,e,fdc,ord,iflag)
struct fdcoefs *fdc;
int n, ord, iflag;
float *a, *b, *c, *d, *e, *x;
{
if(iflag == XDERIV)
   {
   ord2(1,a,x,b,c,&fdc->dtoh);
   a++; b++; c++; d++; e++; x++;

   if(ord == 2)
      ord2(n-2,a,x,b,c,&fdc->dtoh);
   else
      ord4(n-2,a,x,b,c,d,e,&fdc->c0,&fdc->c1);
   a+= (n-2); b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2); x+= (n-2);

   ord2(1,a,x,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2(n,a,x,b,c,&fdc->dtoh);
   else
      ord4(n,a,x,b,c,d,e,&fdc->c0,&fdc->c1);
   }
}

ord2(n,a,x,b,c,c1)  /* a = a + x*c1*(b-c) */
register float *c1, *a, *b, *c, *x;
int n;
{
float con1;

con1 = *c1;
while(n--)
   {
   a[0] = a[0] + x[0]*con1*(b[0]-c[0]);
   a++; b++; c++; x++;
   }
}

ord4(n,a,x,b,c,d,e,c1,c2) /* a = a + x*(c1*(b-c) + c2*(d-e)) */
register float *c1, *c2, *a, *b, *c, *d, *e, *x;
int n;
{
float con1, con2;
 
con1 = *c1;       
con2 = *c2;
while(n--)
   {
   a[0] = a[0] + x[0]*(con1*(b[0]-c[0]) + con2*(d[0]-e[0]));
   a++; b++; c++; d++; e++; x++;
   }
}

diffx(n,a1,a2,a3,x1,x2,b,c,d,e,fdc,ord,iflag)
struct fdcoefs *fdc;
int n, ord, iflag;
float *a1, *a2, *a3, *b, *c, *d, *e, *x1, *x2;
{
if(iflag == XDERIV)
   {
   ord2x(1,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   a1++; a2++; a3++;
   b++; c++; d++; e++; x1++; x2++;

   if(ord == 2)
      ord2x(n-2,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   else
      ord4x(n-2,a1,a2,a3,x1,x2,b,c,d,e,&fdc->c0,&fdc->c1);

   a1+= (n-2); a2+= (n-2); a3+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x1+= (n-2); x2+= (n-2);

   ord2x(1,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2x(n,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   else
      ord4x(n,a1,a2,a3,x1,x2,b,c,d,e,&fdc->c0,&fdc->c1);
   }
}

ord2x(n,a1,a2,a3,x1,x2,b,c,d1)
register float *d1, *a1, *a2, *a3, *b, *c, *x1, *x2;
int n;
{
float tmp;

while(n--)
   {
   tmp = (*d1)*(b[0]-c[0]);

   a1[0] = a1[0] + x1[0]*tmp;
   a2[0] = a2[0] + x2[0]*tmp;
   a3[0] = a3[0] + x2[0]*tmp;

   a1++; a2++; a3++;
   b++; c++;
   x1++; x2++;
   }
}

ord4x(n,a1,a2,a3,x1,x2,b,c,d,e,d1,d2)
register float *d1, *d2, *a1, *a2, *a3;
register float *b, *c, *d, *e, *x1, *x2;
int n;
{
int i;
float tmp;

while(n--)
   {
   tmp = (*d1)*(b[0]-c[0]) + (*d2)*(d[0]-e[0]);

   a1[0] = a1[0] + x1[0]*tmp;
   a2[0] = a2[0] + x2[0]*tmp;
   a3[0] = a3[0] + x2[0]*tmp;

   a1++; a2++; a3++;
   b++; c++; d++; e++;
   x1++; x2++;
   }
}

diff_mv(n,a,cp,x,ks,b,c,d,e,fdc,ord,iflag)
struct fdcoefs *fdc;
int n, ord, iflag;
float *a, *cp, *b, *c, *d, *e, *x, *ks;
{
if(iflag == XDERIV)
   {
   ord2_mv(1,a,cp,x,ks,b,c,&fdc->dtoh);
   a++; cp++; b++; c++; d++; e++; x++; ks++;

   if(ord == 2)
      ord2_mv(n-2,a,cp,x,ks,b,c,&fdc->dtoh);
   else
      ord4_mv(n-2,a,cp,x,ks,b,c,d,e,&fdc->c0,&fdc->c1);

   a+= (n-2); cp+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x+= (n-2); ks+= (n-2);

   ord2_mv(1,a,cp,x,ks,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2_mv(n,a,cp,x,ks,b,c,&fdc->dtoh);
   else
      ord4_mv(n,a,cp,x,ks,b,c,d,e,&fdc->c0,&fdc->c1);
   }
}

ord2_mv(n,a,cp,x,ks,b,c,d1)
register float *d1, *a, *cp, *b, *c, *x, *ks;
int n;
{
float tmp, con1;

con1 = *d1;
while(n--)
   {
/* put effective media parameter into tmp (thus into ks as well) */
   tmp = x[0]*con1*(b[0]-c[0]);

   /* OLD WAY
   a[0] = a[0] + x[0]*tmp;
   */
   a[0] = a[0] + tmp;
   cp[0] = cp[0] + ks[0]*tmp;

   a++; cp++;
   b++; c++;
   x++; ks++;
   }
}

ord4_mv(n,a,cp,x,ks,b,c,d,e,d1,d2)
register float *d1, *d2, *a, *cp, *b, *c, *d, *e, *x, *ks;
int n;
{
float tmp, con1, con2;
 
con1 = *d1;
con2 = *d2;
while(n--)
   {
/* put effective media parameter into tmp (thus into ks as well) */
   tmp = x[0]*(con1*(b[0]-c[0]) + con2*(d[0]-e[0]));

   /* OLD WAY
   a[0] = a[0] + x[0]*tmp;
   */
   a[0] = a[0] + tmp;
   cp[0] = cp[0] + ks[0]*tmp;

   a++; cp++;
   b++; c++; d++; e++;
   x++; ks++;
   }
}

diffx_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,fdc,ord,iflag)
struct fdcoefs *fdc;
int n, ord, iflag;
float *a1, *a2, *a3, *c1, *c2, *c3, *b, *c, *d, *e, *x1, *x2, *kp, *ks;
{
if(iflag == XDERIV)
   {
   ord2x_mv(1,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);

   a1++; a2++; a3++; c1++; c2++; c3++;
   b++; c++; d++; e++; x1++; x2++; kp++; ks++;

   if(ord == 2)
      ord2x_mv(n-2,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   else
      ord4x_mv(n-2,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,&fdc->c0,&fdc->c1);

   a1+= (n-2); a2+= (n-2); a3+= (n-2);
   c1+= (n-2); c2+= (n-2); c3+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x1+= (n-2); x2+= (n-2);
   kp+= (n-2); ks+= (n-2);

   ord2x_mv(1,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2x_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   else
      ord4x_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,&fdc->c0,&fdc->c1);
   }
}

ord2x_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d1)
register float *d1, *a1, *a2, *a3, *c1, *c2, *c3, *b, *c, *x1, *x2, *kp, *ks;
int n;
{
float con1;
float tmp;
float two = 2.0;

while(n--)
   {
   tmp = (*d1)*(b[0]-c[0]);
   con1 = x2[0]*tmp;

   a1[0] = a1[0] + x1[0]*tmp;
   a2[0] = a2[0] + con1;
   a3[0] = a3[0] + con1;

   /* OLD WAY
   con1 = (kp[0] - two*ks[0])*tmp;
   */

/* need to apply mu to ks coefficient, calculate from (lam+2mu)-lam,
factor of 2 carries through formula */
   con1 = (kp[0] - (x1[0]-x2[0])*ks[0])*tmp;

   c1[0] = c1[0] + kp[0]*tmp;
   c2[0] = c2[0] + con1;
   c3[0] = c3[0] + con1;

   a1++; a2++; a3++;
   c1++; c2++; c3++;
   b++; c++;
   x1++; x2++;
   kp++; ks++;
   }
}

ord4x_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,d1,d2)
float *d1, *d2, *a1, *a2, *a3, *c1, *c2, *c3;
float *b, *c, *d, *e, *x1, *x2, *kp, *ks;
int n;
{
int i;
float con1, tmp;
float two = 2.0;

while(n--)
   {
   tmp = (*d1)*(b[0]-c[0]) + (*d2)*(d[0]-e[0]);
   con1 = x2[0]*tmp;

   a1[0] = a1[0] + x1[0]*tmp;
   a2[0] = a2[0] + con1;
   a3[0] = a3[0] + con1;

   /* OLD WAY
   con1 = (kp[0] - two*ks[0])*tmp;
   */

/* need to apply mu to ks coefficient, calculate from (lam+2mu)-lam,
factor of 2 carries through formula */
   con1 = (kp[0] - (x1[0]-x2[0])*ks[0])*tmp;

   c1[0] = c1[0] + kp[0]*tmp;
   c2[0] = c2[0] + con1;
   c3[0] = c3[0] + con1;

   a1++; a2++; a3++;
   c1++; c2++; c3++;
   b++; c++; d++; e++;
   x1++; x2++;
   kp++; ks++;
   }
}

set_vminvmax_reedmod(mds,medf,fdc,nx,ny,nz,vsmin)
struct modelstorage *mds;
struct fdcoefs *fdc;
float *medf, *vsmin;
int nx, ny, nz;
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *mu, *invrho;
float a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;
float one = 1.0;

/*
   Set lowest valid shear velocity to be 10 m/sec (0.01 km/sec).
   This avoids problems in regions where mu -> 0.  The variable 'b2floor'
   is set to 0.01*0.01 = 0.0001 as a cutoff to find 'vmin'.
*/

float b2floor = 0.0001;

/*
   Find vmin and vmax for entire model grid.
*/

init_model_seek(mds,MEDIA_FIELD);
mptr = medf;

vmin = 1.0e+15;
vmax = -1.0;
for(iy=0;iy<ny;iy++)
   {
   mptr = reed_model(mds,mptr,MEDIA_FIELD);

   lam2mu = mptr;
   mu     = mptr + 2*nx*nz;
   invrho = mptr + 5*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;
         }
      }

   lam2mu = mptr;
   mu     = mptr + 2*nx*nz;
   invrho = mptr + 6*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;
         }
      }

   lam2mu = mptr;
   mu     = mptr + 2*nx*nz;
   invrho = mptr + 7*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;
         }
      }
   }

fdc->vmin = sqrt(vmin);
fdc->vmax = sqrt(vmax);

init_model_seek(mds,MEDIA_FIELD);
}

set_izord2(mds,medf,fdc,nx,ny,nz,vsmin)
struct modelstorage *mds;
struct fdcoefs *fdc;
float *medf, *vsmin;
int nx, ny, nz;
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *mu, *invrho;
float a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;
float one = 1.0;

/*
   Find 'izord2'.  This is the depth index below which all velocities are
   at least two times larger than 'vmin' (remember, 'vmin' is stored as
   velocity squared at this point).  For this part of the model,
   second-order FD operators (ord2) can be used without substantial loss
   of accuracy, but with added computational efficiency.
 
   To be tested: 05 May 1994, RWG
*/

init_model_seek(mds,MEDIA_FIELD);
mptr = medf;

if(*vsmin < 0.0)
   *vsmin = fdc->vmin;

vmin2 = 4.0*(*vsmin)*(*vsmin);   /* 2*vmin squared */

izmin = 0;
for(iy=0;iy<ny;iy++)
   {
   mptr = reed_model(mds,mptr,MEDIA_FIELD);

   lam2mu = mptr;
   mu     = mptr + 2*nx*nz;
   invrho = mptr + 5*nx*nz;

   for(ix=0;ix<nx;ix++)
      {
      iz = nz - 1;
      j = ix + iz*nx;
      b2 = mu[j]*invrho[j];
 
      while(b2 > vmin2 && iz > 0)  /* store as velocity squared for now */
	 {
	 iz--;
         j = ix + iz*nx;
         b2 = mu[j]*invrho[j];
         }
 
      if((iz+1) > izmin)
	 {
	 ixmin = ix;
	 iymin = iy;
         izmin = iz + 1;
	 }
      }
   }
 
fdc->izord2 = izmin + 1;  /* add 1 just to be sure */

init_model_seek(mds,MEDIA_FIELD);
}

mv_mult1(t,c,a,n)
register float *t, *c, *a;
int n;
{
while(n--)
   {
   t[0] = t[0] - c[0];
   c[0] = c[0]*a[0];

   t++; c++; a++;
   }
}

mv_mult1x(t,c,n)
register float *t, *c;
int n;
{
while(n--)
   {
   t[0] = t[0] - c[0];
   t++; c++;
   }
}

mv_mult1y(c,a,n)
register float *c, *a;
int n;
{
while(n--)
   {
   c[0] = c[0]*a[0];
   c++; a++;
   }
}

mv_mult2(t,c,n)
register float *t, *c;
int n;
{
while(n--)
   {
   t[0] = t[0] - c[0];
   t++; c++;
   }
}

qmult_null(qf,v0,n)
float *qf, *v0;
int n;
{
return(-1);
}

void set_vminvmax_reedmodP3(float *medf,struct fdcoefs *fdc,int nx,int ny,int nz,float *vsmin)
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *lam, *mu1, *mu2, *mu3, *invrho;
float *lraw, *mraw;
float a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;
float one = 1.0;

/*
   Set lowest valid shear velocity to be 10 m/sec (0.01 km/sec).
   This avoids problems in regions where mu -> 0.  The variable 'b2floor'
   is set to 0.01*0.01 = 0.0001 as a cutoff to find 'vmin'.
*/

float b2floor = 0.0001;

/*
   Find vmin and vmax for entire model grid.
*/

vmin = 1.0e+15;
vmax = -1.0;
for(iy=0;iy<ny;iy++)
   {
   mptr = medf + iy*N_MED_VARS*nx*nz;

   lam2mu = mptr;
   lam    = mptr +   nx*nz;
   mu1    = mptr + 2*nx*nz;
   mu2    = mptr + 3*nx*nz;
   mu3    = mptr + 4*nx*nz;
   invrho = mptr + 5*nx*nz;

   lraw = mptr + 11*nx*nz;
   mraw = mptr + 12*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=1;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;

/*
*/
if(lam2mu[j] <= 0.0 || lam[j] <= 0.0 || mu1[j] <= 0.0 || mu2[j] <= 0.0 || mu3[j] <= 0.0)
   {
   fprintf(stderr,"ERROR: %6d %6d %6d\n",ix,iy,iz);
   fprintf(stderr,"%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n",lam2mu[j],lam[j],mu1[j],mu2[j],mu3[j],lraw[j],mraw[j],1.0/invrho[j]);
   fflush(stderr);
   exit(-1);
   }

         }
      }

   lam2mu = mptr;
   mu1    = mptr + 2*nx*nz;
   invrho = mptr + 6*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;
         }
      }

   lam2mu = mptr;
   mu1    = mptr + 2*nx*nz;
   invrho = mptr + 7*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
         if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
            vmin = b2;
         if(a2 > vmax)    /* store as velocity squared for now */
            vmax = a2;
         }
      }
   }

fdc->vmin = sqrt(vmin);
fdc->vmax = sqrt(vmax);
}

void set_izord2P3(float *medf,struct fdcoefs *fdc,int nx,int ny,int nz,float *vsmin)
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *mu, *invrho;
float a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;
float one = 1.0;

/*
   Find 'izord2'.  This is the depth index below which all velocities are
   at least two times larger than 'vmin' (remember, 'vmin' is stored as
   velocity squared at this point).  For this part of the model,
   second-order FD operators (ord2) can be used without substantial loss
   of accuracy, but with added computational efficiency.
 
   To be tested: 05 May 1994, RWG
*/

if(*vsmin < 0.0)
   *vsmin = fdc->vmin;

vmin2 = 4.0*(*vsmin)*(*vsmin);   /* 2*vmin squared */

izmin = 0;
for(iy=0;iy<ny;iy++)
   {
   mptr = medf + iy*N_MED_VARS*nx*nz;

   lam2mu = mptr;
   mu     = mptr + 2*nx*nz;
   invrho = mptr + 5*nx*nz;

   for(ix=0;ix<nx;ix++)
      {
      iz = nz - 1;
      j = ix + iz*nx;
      b2 = mu[j]*invrho[j];
 
      while(b2 > vmin2 && iz > 0)  /* store as velocity squared for now */
	 {
	 iz--;
         j = ix + iz*nx;
         b2 = mu[j]*invrho[j];
         }
 
      if((iz+1) > izmin)
	 {
	 ixmin = ix;
	 iymin = iy;
         izmin = iz + 1;
	 }
      }
   }
 
fdc->izord2 = izmin + 1;  /* add 1 just to be sure */
}

void get_vel_minvmaxP3(float *medf,struct vel_minmax *vpars,int nx,int ny,int nz,struct runparamsP3 *rpars)
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *lam, *mu1, *mu2, *mu3, *invrho;
float *lraw, *mraw;
float a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;

struct nodeinfo *ni;

/*
   Set lowest valid shear velocity to be 10 m/sec (0.01 km/sec).
   This avoids problems in regions where mu -> 0.  The variable 'b2floor'
   is set to 0.01*0.01 = 0.0001 as a cutoff to find 'vmin'.
*/

float b2floor = 0.0001;

/*
   Find vpmax, vpmin, vsmax, and vsmin for entire model grid.
*/

ni = &(rpars->ni);

vpars->vpmax = -1.0;
vpars->vpmin = 1.0e+15;
vpars->vsmax = -1.0;
vpars->vsmin = 1.0e+15;
vpars->vpvsmax = -1.0;
vpars->vpvsmin = 1.0e+15;

for(iy=0;iy<ny;iy++)
   {
   mptr = medf + iy*N_MED_VARS*nx*nz;

   lam2mu = mptr;
   lam    = mptr +   nx*nz;
   mu1    = mptr + 2*nx*nz;
   mu2    = mptr + 3*nx*nz;
   mu3    = mptr + 4*nx*nz;
   invrho = mptr + 5*nx*nz;

   lraw = mptr + 11*nx*nz;
   mraw = mptr + 12*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=1;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
/* store as velocity squared for now */

/*
         if(b2 > b2floor && a2/b2 < vpars->vpvsmin)
	    {
            vpars->vpvsmin = a2/b2;
            vpars->ix_vpvsmin = ix + ni->nx1;
            vpars->iy_vpvsmin = iy + ni->ny1;
            vpars->iz_vpvsmin = iz + ni->nz1;
	    }

         if(a2/b2 > vpars->vpvsmax)
	    {
            vpars->vpvsmax = a2/b2;
            vpars->ix_vpvsmax = ix + ni->nx1;
            vpars->iy_vpvsmax = iy + ni->ny1;
            vpars->iz_vpvsmax = iz + ni->nz1;
	    }
*/

         if(b2 > b2floor && b2 < vpars->vsmin)
	    {
            vpars->vsmin = b2;
            vpars->ix_vsmin = ix + ni->nx1;
            vpars->iy_vsmin = iy + ni->ny1;
            vpars->iz_vsmin = iz + ni->nz1;
	    }

         if(b2 > vpars->vsmax)
	    {
            vpars->vsmax = b2;
            vpars->ix_vsmax = ix + ni->nx1;
            vpars->iy_vsmax = iy + ni->ny1;
            vpars->iz_vsmax = iz + ni->nz1;
	    }

         if(a2 < vpars->vpmin)
	    {
            vpars->vpmin = a2;
            vpars->ix_vpmin = ix + ni->nx1;
            vpars->iy_vpmin = iy + ni->ny1;
            vpars->iz_vpmin = iz + ni->nz1;
	    }

         if(a2 > vpars->vpmax)
	    {
            vpars->vpmax = a2;
            vpars->ix_vpmax = ix + ni->nx1;
            vpars->iy_vpmax = iy + ni->ny1;
            vpars->iz_vpmax = iz + ni->nz1;
	    }
         }
      }

   lam2mu = mptr;
   mu1    = mptr + 2*nx*nz;
   invrho = mptr + 6*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
/* store as velocity squared for now */

/*
         if(b2 > b2floor && a2/b2 < vpars->vpvsmin)
	    {
            vpars->vpvsmin = a2/b2;
            vpars->ix_vpvsmin = ix + ni->nx1;
            vpars->iy_vpvsmin = iy + ni->ny1;
            vpars->iz_vpvsmin = iz + ni->nz1;
	    }

         if(a2/b2 > vpars->vpvsmax)
	    {
            vpars->vpvsmax = a2/b2;
            vpars->ix_vpvsmax = ix + ni->nx1;
            vpars->iy_vpvsmax = iy + ni->ny1;
            vpars->iz_vpvsmax = iz + ni->nz1;
	    }
*/

         if(b2 > b2floor && b2 < vpars->vsmin)
	    {
            vpars->vsmin = b2;
            vpars->ix_vsmin = ix + ni->nx1;
            vpars->iy_vsmin = iy + ni->ny1;
            vpars->iz_vsmin = iz + ni->nz1;
	    }

         if(b2 > vpars->vsmax)
	    {
            vpars->vsmax = b2;
            vpars->ix_vsmax = ix + ni->nx1;
            vpars->iy_vsmax = iy + ni->ny1;
            vpars->iz_vsmax = iz + ni->nz1;
	    }

         if(a2 < vpars->vpmin)
	    {
            vpars->vpmin = a2;
            vpars->ix_vpmin = ix + ni->nx1;
            vpars->iy_vpmin = iy + ni->ny1;
            vpars->iz_vpmin = iz + ni->nz1;
	    }

         if(a2 > vpars->vpmax)
	    {
            vpars->vpmax = a2;
            vpars->ix_vpmax = ix + ni->nx1;
            vpars->iy_vpmax = iy + ni->ny1;
            vpars->iz_vpmax = iz + ni->nz1;
	    }
         }
      }

   lam2mu = mptr;
   mu1    = mptr + 2*nx*nz;
   invrho = mptr + 7*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;
         a2 = lam2mu[j]*invrho[j];
         b2 = mu1[j]*invrho[j];
 
/* store as velocity squared for now */

/*
         if(b2 > b2floor && a2/b2 < vpars->vpvsmin)
	    {
            vpars->vpvsmin = a2/b2;
            vpars->ix_vpvsmin = ix + ni->nx1;
            vpars->iy_vpvsmin = iy + ni->ny1;
            vpars->iz_vpvsmin = iz + ni->nz1;
	    }

         if(a2/b2 > vpars->vpvsmax)
	    {
            vpars->vpvsmax = a2/b2;
            vpars->ix_vpvsmax = ix + ni->nx1;
            vpars->iy_vpvsmax = iy + ni->ny1;
            vpars->iz_vpvsmax = iz + ni->nz1;
	    }
*/

         if(lam2mu[j] > 0.0 && lam2mu[j]/mraw[j] < vpars->vpvsmin)
	    {
            vpars->vpvsmin = lam2mu[j]/mraw[j];
            vpars->ix_vpvsmin = ix + ni->nx1;
            vpars->iy_vpvsmin = iy + ni->ny1;
            vpars->iz_vpvsmin = iz + ni->nz1;
	    }

         if(lam2mu[j] > 0.0 && lam2mu[j]/mraw[j] > vpars->vpvsmax)
	    {
            vpars->vpvsmax = lam2mu[j]/mraw[j];
            vpars->ix_vpvsmax = ix + ni->nx1;
            vpars->iy_vpvsmax = iy + ni->ny1;
            vpars->iz_vpvsmax = iz + ni->nz1;
	    }

         if(b2 > b2floor && b2 < vpars->vsmin)
	    {
            vpars->vsmin = b2;
            vpars->ix_vsmin = ix + ni->nx1;
            vpars->iy_vsmin = iy + ni->ny1;
            vpars->iz_vsmin = iz + ni->nz1;
	    }

         if(b2 > vpars->vsmax)
	    {
            vpars->vsmax = b2;
            vpars->ix_vsmax = ix + ni->nx1;
            vpars->iy_vsmax = iy + ni->ny1;
            vpars->iz_vsmax = iz + ni->nz1;
	    }

         if(a2 < vpars->vpmin)
	    {
            vpars->vpmin = a2;
            vpars->ix_vpmin = ix + ni->nx1;
            vpars->iy_vpmin = iy + ni->ny1;
            vpars->iz_vpmin = iz + ni->nz1;
	    }

         if(a2 > vpars->vpmax)
	    {
            vpars->vpmax = a2;
            vpars->ix_vpmax = ix + ni->nx1;
            vpars->iy_vpmax = iy + ni->ny1;
            vpars->iz_vpmax = iz + ni->nz1;
	    }
         }
      }
   }

vpars->vpmin = sqrt(vpars->vpmin);
vpars->vpmax = sqrt(vpars->vpmax);
vpars->vsmin = sqrt(vpars->vsmin);
vpars->vsmax = sqrt(vpars->vsmax);
vpars->vpvsmin = sqrt(vpars->vpvsmin);
vpars->vpvsmax = sqrt(vpars->vpvsmax);
}

void setcoefs_pvc(float *pvfield,struct nodeinfo *ni,int fs,int order,struct fdcoefs *coefs,float *dt,float *h,float *vsmin)
{
float sum, f, con, vmn, vmx;
float lmin, fmax;

float *pvptr, *c0ptr, *c1ptr;
int ix, iy, iz, ip, off0, off1, off2;
int loc_izord2, nx, ny, nz;

int npad_2ord = NBND_PAD_2NDORDER;

if(NBND_PAD_COPY > npad_2ord)
   npad_2ord = NBND_PAD_COPY;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

coefs->order = order;
coefs->ordx = ORDERX;      /* spatial differencing order below iz=izord2 */

sum = 0.0;
vmn = coefs->vmin;
vmx = coefs->vmax;
con = (*dt)/(*h);

coefs->dtoh = con;
coefs->c0 = 9.0/8.0;
coefs->c1 = -1.0/24.0;

if(coefs->order == 2)
   sum = 1.0;
else if(coefs->order == 4)
   sum = coefs->c0 - coefs->c1;

f = (*h)/(vmx*sum*sqrt(3.0));
if((*dt) > f)
   {
   fprintf(stderr,"**** Unstable run:\n");
   fprintf(stderr,"      vmax = %f\n",vmx);
   fprintf(stderr,"        dt = %f, must be less than %f\n",(*dt),f);
   fprintf(stderr,"        exiting...\n");
   exit(-1);
   }
else
   {
   fprintf(stderr,"**** Stable run:\n");
   fprintf(stderr,"      vmax = %f\n",vmx);
   fprintf(stderr,"        dt = %f, must be less than %f\n\n",(*dt),f);
   }

coefs->c0 = coefs->dtoh*coefs->c0;
coefs->c1 = coefs->dtoh*coefs->c1;

if(coefs->order == 2)
   lmin = 10*(*h);
else if(coefs->order == 4)
   lmin = 5*(*h);
else
   {
   fprintf(stderr,"**** Invalid order:\n");
   fprintf(stderr,"        order must be either 2 or 4\n");
   fprintf(stderr,"        exiting...\n");
   exit(-1);
   }

fmax = vmn/lmin;

fprintf(stderr,"**** Differencing accuracy:\n");
fprintf(stderr,"                    order= %d\n",coefs->order);
fprintf(stderr,"                     fmax= %f Hz\n",fmax);
fprintf(stderr,"                     vmin= %f km/s\n",vmn);
fprintf(stderr,"                   izord2= %d (vsmin= %f)\n\n",coefs->izord2,(*vsmin));

/*
   RWG 2017-10-19
   set grid point specific FD coefficients,
   2nd order in zone NBND_PAD near absorbing boundaries
*/

/* first set coefs for entire sub-model */

loc_izord2 = coefs->izord2 - ni->nz1;
if(loc_izord2 > nz)
   loc_izord2 = nz;

fprintf(stderr,"FD coefs: entire sub-model\n");
fflush(stderr);

for(iy=0;iy<ny;iy++)
   {
   pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

   c0ptr = pvptr + 15*nx*nz;
   c1ptr = pvptr + 16*nx*nz;

   if(coefs->order == 2)
      {
      for(ip=0;ip<nx*nz;ip++)
         {
         c0ptr[ip] = coefs->dtoh;
         c1ptr[ip] = 0.0;
         }
      }
   else
      {
      for(ip=0;ip<nx*loc_izord2;ip++)
         {
         c0ptr[ip] = coefs->c0;
         c1ptr[ip] = coefs->c1;
         }
      for(ip=nx*loc_izord2;ip<nx*nz;ip++)
         {
         c0ptr[ip] = coefs->dtoh;
         c1ptr[ip] = 0.0;
         }
      }
   }

if(coefs->order == 2 || npad_2ord <= 0)
   return;

/* set coefs to 2nd order near abosrbing boundaries */

/*XXXXXXXXXXXXXXX
XXXXXXXXXXXXXXX */
fprintf(stderr,"FD coefs: now checking boundaries\n");
fflush(stderr);

if(ni->minusId_x < 0)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(iz=0;iz<nz;iz++)
         {
         for(ix=0;ix<npad_2ord;ix++)
            {
	    ip = ix + iz*nx;

            c0ptr[ip] = coefs->dtoh;
            c1ptr[ip] = 0.0;
	    }
	 }
      }
   }

if(ni->plusId_x < 0)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(iz=0;iz<nz;iz++)
         {
         for(ix=nx-npad_2ord;ix<nx;ix++)
            {
	    ip = ix + iz*nx;

            c0ptr[ip] = coefs->dtoh;
            c1ptr[ip] = 0.0;
	    }
	 }
      }
   }

if(!fs && ni->minusId_z < 0)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(iz=0;iz<npad_2ord;iz++)
         {
         for(ix=0;ix<nx;ix++)
            {
	    ip = ix + iz*nx;

            c0ptr[ip] = coefs->dtoh;
            c1ptr[ip] = 0.0;
	    }
	 }
      }
   }

if(ni->plusId_z < 0)
   {
   for(iy=0;iy<ny;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(iz=nz-npad_2ord;iz<nz;iz++)
         {
         for(ix=0;ix<nx;ix++)
            {
	    ip = ix + iz*nx;

            c0ptr[ip] = coefs->dtoh;
            c1ptr[ip] = 0.0;
	    }
	 }
      }
   }

if(ni->minusId_y < 0)
   {
   for(iy=0;iy<npad_2ord;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(ip=0;ip<nx*nz;ip++)
         {
         c0ptr[ip] = coefs->dtoh;
         c1ptr[ip] = 0.0;
	 }
      }
   }

if(ni->plusId_y < 0)
   {
   for(iy=ny-npad_2ord;iy<ny;iy++)
      {
      pvptr = pvfield + N_WAVE_VARS*iy*nx*nz;

      c0ptr = pvptr + 15*nx*nz;
      c1ptr = pvptr + 16*nx*nz;

      for(ip=0;ip<nx*nz;ip++)
         {
         c0ptr[ip] = coefs->dtoh;
         c1ptr[ip] = 0.0;
	 }
      }
   }
}

void diff_pvc(int n,float *a,float *c0,float *c1,float *x,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag)
{
if(iflag == XDERIV)
   {
   ord2(1,a,x,b,c,&fdc->dtoh);
   a++; b++; c++; d++; e++; x++;
   c0++; c1++;

   if(ord == 2)
      ord2(n-2,a,x,b,c,&fdc->dtoh);
   else
      ord4_pvc(n-2,a,x,b,c,d,e,c0,c1);

   a+= (n-2); b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2); x+= (n-2);
   c0 = c0 + (n-2); c1 = c1 + (n-2);

   ord2(1,a,x,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2(n,a,x,b,c,&fdc->dtoh);
   else
      ord4_pvc(n,a,x,b,c,d,e,c0,c1);
   }
}

void ord4_pvc(int n,float *a,float *x,float *b,float *c,float *d,float *e,float *c0,float *c1)
{
int i;
 
for(i=0;i<n;i++)
   a[i] = a[i] + x[i]*(c0[i]*(b[i]-c[i]) + c1[i]*(d[i]-e[i]));
}

void diffx_pvc(int n,float *a1,float *a2,float *a3,float *c0,float *c1,float *x1,float *x2,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag)
{
if(iflag == XDERIV)
   {
   ord2x(1,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   a1++; a2++; a3++;
   b++; c++; d++; e++; x1++; x2++;
   c0++; c1++;

   if(ord == 2)
      ord2x(n-2,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   else
      ord4x_pvc(n-2,a1,a2,a3,x1,x2,b,c,d,e,c0,c1);

   a1+= (n-2); a2+= (n-2); a3+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x1+= (n-2); x2+= (n-2);
   c0+= (n-2); c1+= (n-2);

   ord2x(1,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2x(n,a1,a2,a3,x1,x2,b,c,&fdc->dtoh);
   else
      ord4x_pvc(n,a1,a2,a3,x1,x2,b,c,d,e,c0,c1);
   }
}

void ord4x_pvc(int n,float *a1,float *a2,float *a3,float *x1,float *x2,float *b,float *c,float *d,float *e,float *c0,float *c1)
{
int i;
float tmp;

for(i=0;i<n;i++)
   {
   tmp = c0[i]*(b[i]-c[i]) + c1[i]*(d[i]-e[i]);

   a1[i] = a1[i] + x1[i]*tmp;
   a2[i] = a2[i] + x2[i]*tmp;
   a3[i] = a3[i] + x2[i]*tmp;
   }
}

void diffx_mv_pvc(int n,float *a1,float *a2,float *a3,float *d0,float *d1,float *c1,float *c2,float *c3,float *x1,float *x2,float *kp,float *ks,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag)
{
if(iflag == XDERIV)
   {
   ord2x_mv(1,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);

   a1++; a2++; a3++; c1++; c2++; c3++;
   b++; c++; d++; e++; x1++; x2++; kp++; ks++;
   d0++; d1++;

   if(ord == 2)
      ord2x_mv(n-2,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   else
      ord4x_mv_pvc(n-2,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,d0,d1);

   a1+= (n-2); a2+= (n-2); a3+= (n-2);
   c1+= (n-2); c2+= (n-2); c3+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x1+= (n-2); x2+= (n-2);
   kp+= (n-2); ks+= (n-2);
   d0+= (n-2); d1+= (n-2);

   ord2x_mv(1,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2x_mv(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,&fdc->dtoh);
   else
      ord4x_mv_pvc(n,a1,a2,a3,c1,c2,c3,x1,x2,kp,ks,b,c,d,e,d0,d1);
   }
}

void ord4x_mv_pvc(int n,float *a1,float *a2,float *a3,float *c1,float *c2,float *c3,float *x1,float *x2,float *kp,float *ks,float *b,float *c,float *d,float *e,float *d0,float *d1)
{
int i;
float con1, con2, tmp;

for(i=0;i<n;i++)
   {
   tmp = d0[i]*(b[i]-c[i]) + d1[i]*(d[i]-e[i]);
   con1 = x2[i]*tmp;

   a1[i] = a1[i] + x1[i]*tmp;
   a2[i] = a2[i] + con1;
   a3[i] = a3[i] + con1;

   con2 = (kp[i] - (x1[i]-x2[i])*ks[i])*tmp;

   c1[i] = c1[i] + kp[i]*tmp;
   c2[i] = c2[i] + con2;
   c3[i] = c3[i] + con2;
   }
}

void diff_mv_pvc(int n,float *a,float *c0,float *c1,float *cp,float *x,float *ks,float *b,float *c,float *d,float *e,struct fdcoefs *fdc,int ord,int iflag)
{
if(iflag == XDERIV)
   {
   ord2_mv(1,a,cp,x,ks,b,c,&fdc->dtoh);
   a++; cp++; b++; c++; d++; e++; x++; ks++;
   c0++; c1++;

   if(ord == 2)
      ord2_mv(n-2,a,cp,x,ks,b,c,&fdc->dtoh);
   else
      ord4_mv_pvc(n-2,a,cp,x,ks,b,c,d,e,c0,c1);

   a+= (n-2); cp+= (n-2);
   b+= (n-2); c+= (n-2); d+= (n-2); e+= (n-2);
   x+= (n-2); ks+= (n-2);
   c0+= (n-2); c1+= (n-2);

   ord2_mv(1,a,cp,x,ks,b,c,&fdc->dtoh);
   }
else if(iflag == YDERIV || iflag == ZDERIV)
   {
   if(ord == 2)
      ord2_mv(n,a,cp,x,ks,b,c,&fdc->dtoh);
   else
      ord4_mv_pvc(n,a,cp,x,ks,b,c,d,e,c0,c1);
   }
}

void ord4_mv_pvc(int n,float *a,float *cp,float *x,float *ks,float *b,float *c,float *d,float *e,float *c0,float *c1)
{
int i;
float tmp;
 
for(i=0;i<n;i++)
   {
   tmp = x[i]*(c0[i]*(b[i]-c[i]) + c1[i]*(d[i]-e[i]));

   a[i] = a[i] + tmp;
   cp[i] = cp[i] + ks[i]*tmp;
   }
}

void get_vel_minmaxP3(float *medf,struct vel_minmax *vpars,int nx,int ny,int nz,struct runparamsP3 *rpars)
{
int ix, iy, iz, j;
float *mptr, *lam2mu, *lam, *mu1, *mu2, *mu3, *rho;
float *lraw, *mraw, *rraw;
float l2mraw, muavg, a2, b2, den, vmin, vmin2, vmax;
int ixmin, iymin, izmin;

struct nodeinfo *ni;

/*
   Set lowest valid shear velocity to be 10 m/sec (0.01 km/sec).
   This avoids problems in regions where mu -> 0.  The variable 'b2floor'
   is set to 0.01*0.01 = 0.0001 as a cutoff to find 'vmin'.
*/

float b2floor = 0.0001;

/*
   Find vpmax, vpmin, vsmax, and vsmin for entire model grid.
*/

ni = &(rpars->ni);

vpars->vpmax = -1.0;
vpars->vpmin = 1.0e+15;
vpars->vsmax = -1.0;
vpars->vsmin = 1.0e+15;
vpars->vpvsmax = -1.0;
vpars->vpvsmin = 1.0e+15;

for(iy=0;iy<ny;iy++)
   {
   mptr = medf + iy*N_MED_VARS*nx*nz;

   lam2mu = mptr;
   lam    = mptr +   nx*nz;
   mu1    = mptr + 2*nx*nz;
   mu2    = mptr + 3*nx*nz;
   mu3    = mptr + 4*nx*nz;

   lraw = mptr + 11*nx*nz;
   mraw = mptr + 12*nx*nz;
   rho    = mptr + 13*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 j = ix + iz*nx;

	 muavg = 3.0/(1.0/mu1[j] + 1.0/mu2[j] + 1.0/mu3[j]);
	 muavg = (mu1[j] + mu2[j] + mu3[j])/3.0;
	 muavg = 0.5*(lam2mu[j] - lam[j]);

         a2 = lam2mu[j]/rho[j];
         b2 = muavg/rho[j];
 
/* store as velocity squared for now */

         if(lam2mu[j] > 0.0 && lam2mu[j]/muavg < vpars->vpvsmin)
	    {
            vpars->vpvsmin = lam2mu[j]/muavg;
            vpars->ix_vpvsmin = ix + ni->nx1;
            vpars->iy_vpvsmin = iy + ni->ny1;
            vpars->iz_vpvsmin = iz + ni->nz1;
	    }

         if(lam2mu[j] > 0.0 && lam2mu[j]/muavg > vpars->vpvsmax)
	    {
            vpars->vpvsmax = lam2mu[j]/muavg;
            vpars->ix_vpvsmax = ix + ni->nx1;
            vpars->iy_vpvsmax = iy + ni->ny1;
            vpars->iz_vpvsmax = iz + ni->nz1;
	    }

         if(b2 > b2floor && b2 < vpars->vsmin)
	    {
            vpars->vsmin = b2;
            vpars->ix_vsmin = ix + ni->nx1;
            vpars->iy_vsmin = iy + ni->ny1;
            vpars->iz_vsmin = iz + ni->nz1;
	    }

         if(b2 > vpars->vsmax)
	    {
            vpars->vsmax = b2;
            vpars->ix_vsmax = ix + ni->nx1;
            vpars->iy_vsmax = iy + ni->ny1;
            vpars->iz_vsmax = iz + ni->nz1;
	    }

         if(a2 < vpars->vpmin)
	    {
            vpars->vpmin = a2;
            vpars->ix_vpmin = ix + ni->nx1;
            vpars->iy_vpmin = iy + ni->ny1;
            vpars->iz_vpmin = iz + ni->nz1;
	    }

         if(a2 > vpars->vpmax)
	    {
            vpars->vpmax = a2;
            vpars->ix_vpmax = ix + ni->nx1;
            vpars->iy_vpmax = iy + ni->ny1;
            vpars->iz_vpmax = iz + ni->nz1;
	    }
         }
      }
   }

vpars->vpmin = sqrt(vpars->vpmin);
vpars->vpmax = sqrt(vpars->vpmax);
vpars->vsmin = sqrt(vpars->vsmin);
vpars->vsmax = sqrt(vpars->vsmax);
vpars->vpvsmin = sqrt(vpars->vpvsmin);
vpars->vpvsmax = sqrt(vpars->vpvsmax);
}
