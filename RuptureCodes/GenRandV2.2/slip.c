#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void init_slip(struct complex *sc,int nx,int ny,float *st,float *bt)
{
float xdamp, ydamp;
int i;

for(i=0;i<nx*ny;i++)
   {
   sc[i].re = 0.5;
   sc[i].im = 0.0;
   }

/*
taper_slip(sc,nx,ny,st,bt);
*/
}

void taper_slip(struct complex *sc,int nx,int ny,float *st,float *bt)
{
float xdamp, ydamp;
int ix, iy, xb, yb;

xb = (int)((*st)*nx + 0.5);
if(xb < 0)
   xb = 0;

yb = (int)((*bt)*ny + 0.5);
if(yb < 0)
   yb = 0;

for(iy=0;iy<yb;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yb);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*ydamp*sc[ix+iy*nx].re;
      sc[ix+(ny-1-iy)*nx].re = xdamp*ydamp*sc[ix+(ny-1-iy)*nx].re;
      }
   }

for(iy=yb;iy<ny-yb;iy++)
   {
   for(ix=0;ix<xb;ix++)
      {
      xdamp = (float)(ix+1)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*sc[ix+iy*nx].re;
      sc[(nx-1-ix)+iy*nx].re = xdamp*sc[(nx-1-ix)+iy*nx].re;
      }
   }
}

void scale_slip(struct pointsource *ps,struct complex *cs,int nx,int ny,int nys,float *dx,float *dy,float *dtop,float *dip,float *mom,struct velmodel *vm,float *savg,float *smax)
{
float sum, fac, sinD, area;
float zz;
int i, j, k;

float rperd = 0.017453293;

if(*savg < 0.0)
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      zz = (*dtop) + sinD*(j + 0.5)*(*dy);

      k = 0;
      while(zz > vm->dep[k] && k < (vm->nlay)-1)
         k++;

      fac = area*vm->mu[k];
      for(i=0;i<nx;i++)
         sum = sum + fac*cs[i + (j+nys)*nx].re;
      }

   fac = (*mom)/sum;

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*cs[i + (j+nys)*nx].re;

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);
   }
else
   {
   sinD = sin((*dip)*rperd);
   area = (*dx)*(*dy)*1.0e+10;    /* in CMS units */

   sum = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         sum = sum + cs[i + (j+nys)*nx].re;
         }
      }

   fac = (*savg)*(float)(nx*ny)/(sum);

   *smax = 0.0;
   *savg = 0.0;
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         ps[i + j*nx].slip = fac*cs[i + (j+nys)*nx].re;

         *savg = *savg + ps[i + j*nx].slip;
         if(ps[i + j*nx].slip > *smax)
            *smax = ps[i + j*nx].slip;
         }
      }
   *savg = (*savg)/(float)(nx*ny);
   }
}

void kfilt(struct complex *s0,int nx0,int ny0,float *dkx,float *dky,float *xl,float *yl,long *seed,int kflag)
{
int i, j, ip;
float kx, ky, fac, amp, amp0, phs, xl2, yl2;
float phs1, fac1, wtS, wtD;
float xp, k2, invkc2;
int ndkc;

float pi = 3.14159265;
float hcoef;

hcoef = 1.80; /* H=0.80, hcoef = H + 1 */
hcoef = 1.75; /* H=0.75, hcoef = H + 1 */

amp0 = sqrt(s0[0].re*s0[0].re + s0[0].im*s0[0].im);

xl2 = (*xl)*(*xl);
yl2 = (*yl)*(*yl);

/*

   Transition between deterministic and stochastic parts of spectrum
   are given by

       F = wtS*stoch + wtD*deter

   with
    
       wtD = {1 + k2/Kc2}^-(xp)     (kind of a butterworth filter)
       wtS = 1 - wtD

   and

       k2 = kx*kx + ky*ky          (k-squared)
       Kc2 = (N*dky)*(N*dkx)       (corner wavenumber of transition)

   The parameter N specifies the number of dk's in the corner (somewhat
   like a fraction of the total wavenumber space).  The exponent (xp)
   gives the sharpness of the transition.  Based on very limited
   testing, I came up with

       N = 4
       xp = 2.0

*/

xp = 2.0;
ndkc = 4;
invkc2 = ndkc*(*dky)*ndkc*(*dkx);
invkc2 = 1.0/(invkc2);

for(j=0;j<=ny0/2;j++)  /* only do positive half, then use symmetry */
   {
   if(j <= ny0/2)
      ky = j*(*dky);
   else
      ky = (j - ny0)*(*dky);

   for(i=0;i<nx0;i++)
      {
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      ip = i + j*nx0;

      amp = kx*kx*xl2 + ky*ky*yl2;

      /* default is somerville scaling */
      fac = amp0/sqrt(1.0 + amp*amp);

      if(kflag == MAI_FLAG) /* mai scaling */
         {
         fac = exp((hcoef)*log(1.0+amp));
         fac = amp0/sqrt(fac);
	 }

      if(kflag == SOMERVILLE_FLAG)      /* somerville scaling */
         fac = amp0/sqrt(1.0 + amp*amp);

      phs = pi*sfrand(seed);

      fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

      phs1 = 0.5*pi;
      if(s0[ip].re != 0.0)
         {
         phs1 = atan(s0[ip].im/s0[ip].re);
         if(s0[ip].re < 0.0)
            phs1 = phs1 + pi;
         }
      else if(s0[ip].im < 0.0)
         phs1 = -0.5*pi;

      while(phs1 > pi)
         phs1 = phs1 - pi;
      while(phs1 < -pi)
         phs1 = phs1 + pi;

      k2 = (kx*kx + ky*ky)*invkc2;
      wtD = exp(-xp*log(1.0 + k2));
      wtS = 1.0 - wtD;

      s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
      s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);

/*

   OLD STUFF that I tried FOLLOWS from HERE

*/

/*
   Do not alter phase of lowest (k=0) and 2nd lowest (k=dk) wavenumbers
   so average slip and edge taper are not significantly modified
      if(i > 1 && j > 1)
*/

/* 
   Do not alter phase of lowest (k=0) wavenumbers and use average of
   random and determinstic phase for 2nd lowest (k=dk)
   so average slip and edge taper are not significantly modified
*/

/*
      if(i > 1 && j > 1)
         {
         s0[ip].re = fac*cos(phs);
         s0[ip].im = fac*sin(phs);
	 }
      else if((i == 1 && j > 0) || (i > 0 && j == 1))
         {
	 if(i == 1)
	    wtS = sqrt(ky*ky*invkym2);
	 if(j == 1)
	    wtS = sqrt(kx*kx*invkxm2);

	 if(wtS < 0.5)
	    wtS = 0.5;

	 wtD = 1.0 - wtS;

	 fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

	 phs1 = 0.5*pi;
	 if(s0[ip].re != 0.0)
	    {
	    phs1 = atan(s0[ip].im/s0[ip].re);
	    if(s0[ip].re < 0.0)
	       phs1 = phs1 + pi;
	    }
	 else if(s0[ip].im < 0.0)
	    phs1 = -0.5*pi;

         while(phs1 > pi)
            phs1 = phs1 - pi;
         while(phs1 < -pi)
            phs1 = phs1 + pi;

         s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
         s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);
	 }
*/
      }
   }

/* 
   Enforce Hermitian symmetry to make slip real valued
*/

for(j=1;j<=(ny0-1)/2;j++)
   {
   s0[(ny0-j)*nx0].re = s0[j*nx0].re;
   s0[(ny0-j)*nx0].im = -s0[j*nx0].im;
   }

for(i=1;i<=(nx0-1)/2;i++)
   {
   s0[nx0-i].re = s0[i].re;
   s0[nx0-i].im = -s0[i].im;
   }

for(j=1;j<=ny0/2;j++)
   {
   for(i=1;i<=nx0/2;i++)
      {
      s0[(nx0-i)+(ny0-j)*nx0].re = s0[i+j*nx0].re;
      s0[(nx0-i)+(ny0-j)*nx0].im = -s0[i+j*nx0].im;

      s0[i+(ny0-j)*nx0].re = s0[(nx0-i)+j*nx0].re;
      s0[i+(ny0-j)*nx0].im = -s0[(nx0-i)+j*nx0].im;
      }
   }
}

void fft2d(struct complex *xc,int n1,int n2,int isgn,float *d1,float *d2)
{
int i, j, ip;
float *space;
struct complex *xtc;
float normf;

normf = (*d1)*(*d2);

space = (float *) check_malloc (2*(n1+n2)*sizeof(float));

for(j=0;j<n2;j++)
   fourg_(xc+j*n1,&n1,&isgn,space);

xtc = (struct complex *) check_malloc (n2*sizeof(struct complex));

for(i=0;i<n1;i++)
   {
   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xtc[j].re = xc[ip].re;
      xtc[j].im = xc[ip].im;
      }

   fourg_(xtc,&n2,&isgn,space);

   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xc[ip].re = normf*xtc[j].re;
      xc[ip].im = normf*xtc[j].im;
      }
   }

free(xtc);
free(space);
}
