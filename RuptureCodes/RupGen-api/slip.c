#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "fourg.h"
#include "misc.h"

void init_slip_IO(struct complex *sc,int nx,int ny,float *dx,float *dy,int flip,char *file)
{
FILE *fpr;
char str[1024];
float slip, x0, y0, x1, y1;
int ix0, iy0, ix1, iy1, iystart;
int ix, iy, ip, ip1;

if(file[0] == '\0')
   {
   for(ip=0;ip<nx*ny;ip++)
      {
      sc[ip].re = 0.5;
      sc[ip].im = 0.0;
      }
   }
else
   {
   fpr = _fopfile(file,"r");

   fgets(str,1024,fpr);
   sscanf(str,"%f",&slip);

   for(ip=0;ip<nx*ny;ip++)
      {
      sc[ip].re = slip;
      sc[ip].im = 0.0;
      }

   while(fgets(str,1024,fpr) != NULL)
      {
      sscanf(str,"%f %f %f %f %f",&slip,&x0,&y0,&x1,&y1);

      iystart = 0;
      if(flip == 1)
         iystart = ny/2;

iystart = 0;

      ix0 = (int)((0.5*nx*(*dx) + x0)/(*dx) + 0.5);
      if(ix0 < 0)
         ix0 = 0;
      if(ix0 > nx)
         ix0 = nx;

      iy0 = iystart + (int)(y0/(*dy) + 0.5);
      if(iy0 < iystart)
         iy0 = iystart;
      if(iy0 > ny)
         iy0 = ny;

      ix1 = (int)((0.5*nx*(*dx) + x1)/(*dx) + 0.5);
      if(ix1 < 0)
         ix1 = 0;
      if(ix1 > nx)
         ix1 = nx;

      iy1 = iystart + (int)(y1/(*dy) + 0.5) + 1;
      if(iy1 < iystart)
         iy1 = iystart;
      if(iy1 > ny)
         iy1 = ny;

fprintf(stderr,"ix0= %d ix1= %d iy0= %d iy1= %d\n",ix0,ix1,iy0,iy1);

      for(iy=iy0;iy<iy1;iy++)
         {
         for(ix=ix0;ix<ix1;ix++)
            {
	    ip = ix+iy*nx;
	    sc[ip].re = slip;
	    }
	 }

      /*
      if(flip == 1)
         {
         for(iy=2*iystart-(iy0+1);iy>2*iystart-(iy1+1);iy--)
            {
            for(ix=ix0;ix<ix1;ix++)
               {
	       ip = ix+iy*nx;
	       sc[ip].re = slip;
	       }
	    }
	 }
      */
      }

   if(flip == 1)
      {
      for(iy=0;iy<ny/2;iy++)
         {
         for(ix=0;ix<nx;ix++)
            {
            ip = ix+iy*nx;
            ip1 = ix+((ny-1)-iy)*nx;

            sc[ip1].re = sc[ip].re;
            }
         }
      }

   }
}

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

void taper_slip_all(struct complex *sc,int nx,int ny,float *st,float *bt,float *tt)
{
float xdamp, ydamp;
int ix, iy, xb, yt, yb;

xb = (int)((*st)*nx + 0.5);
if(xb < 0)
   xb = 0;

yt = (int)((*tt)*ny + 0.5);
if(yt < 0)
   yt = 0;

yb = (int)((*bt)*ny + 0.5);
if(yb < 0)
   yb = 0;

for(iy=0;iy<yt;iy++)
   {
   ydamp = (float)(iy+1)/(float)(yt);
   for(ix=0;ix<nx;ix++)
      {
      xdamp = 1.0;
      if(ix < xb)
         xdamp = (float)(ix+1)/(float)(xb);
      if(ix > nx-xb)
         xdamp = (float)(nx-ix)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*ydamp*sc[ix+iy*nx].re;
      }
   }

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

      sc[ix+(ny-1-iy)*nx].re = xdamp*ydamp*sc[ix+(ny-1-iy)*nx].re;
      }
   }

for(iy=yt;iy<ny-yb;iy++)
   {
   for(ix=0;ix<xb;ix++)
      {
      xdamp = (float)(ix+1)/(float)(xb);

      sc[ix+iy*nx].re = xdamp*sc[ix+iy*nx].re;
      sc[(nx-1-ix)+iy*nx].re = xdamp*sc[(nx-1-ix)+iy*nx].re;
      }
   }
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
float fre, fim;
float phsb;
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

xp = 1;
ndkc = 2;

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

      phs = pi*_sfrand(seed);

      fac1 = sqrt(s0[ip].re*s0[ip].re + s0[ip].im*s0[ip].im);

      phs1 = 0.5*pi;
      if(s0[ip].re != (float)(0.0))
         {
         phs1 = atan(s0[ip].im/s0[ip].re);
         if(s0[ip].re < 0.0)
            phs1 = phs1 + pi;
         }
      else if(s0[ip].im < 0.0)
         phs1 = -0.5*pi;

/* 09/24/2009 I have no idea why I put these conditions here, now they're removed
      while(phs1 > pi)
         phs1 = phs1 - pi;
      while(phs1 < -pi)
         phs1 = phs1 + pi;
*/

      k2 = (kx*kx + ky*ky)*invkc2;
      wtD = exp(-xp*log(1.0 + k2));

      k2 = exp(xp*log(kx*kx/(ndkc*(*dkx)*ndkc*(*dkx)))) + exp(xp*log(ky*ky/(ndkc*(*dky)*ndkc*(*dky))));
      wtD = 1.0/(1.0 + k2);
      k2 = kx*kx/(ndkc*(*dkx)*ndkc*(*dkx)) + ky*ky/(ndkc*(*dky)*ndkc*(*dky));
      wtD = 1.0/(1.0 + exp(xp*log(k2)));

      wtS = 1.0 - wtD;

      s0[ip].re = wtS*fac*cos(phs) + wtD*fac1*cos(phs1);
      s0[ip].im = wtS*fac*sin(phs) + wtD*fac1*sin(phs1);

/*
wtD = 0.0;
wtS = 20.0;
s0[ip].re = fac*gaus_rand(&wtS,&wtD,seed)/sqrt(2.0*wtS*wtS);
s0[ip].im = fac*gaus_rand(&wtS,&wtD,seed)/sqrt(2.0*wtS*wtS);
*/

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

/*
for(j=0;j<ny0;j++)
   {
   for(i=0;i<nx0;i++)
      {
      ip = i + j*nx0;

      if(j <= ny0/2)
         ky = j*(*dky);
      else
         ky = (j - ny0)*(*dky);
      if(i <= nx0/2)
         kx = i*(*dkx);
      else
         kx = (i - nx0)*(*dkx);

      phs = -2.0*pi*kx*10.0;

      fre = s0[ip].re*cos(phs) - s0[ip].im*sin(phs);
      fim = s0[ip].re*sin(phs) + s0[ip].im*cos(phs);

      phs = -2.0*pi*ky*10.0;

      s0[ip].re = fre*cos(phs) - fim*sin(phs);
      s0[ip].im = fre*sin(phs) + fim*cos(phs);
      }
   }
*/
}

void fft2d(struct complex *xc,int n1,int n2,int isgn,float *d1,float *d2)
{
int i, j, ip;
float *space;
struct complex *xtc;
float normf;

normf = (*d1)*(*d2);

space = (float *) _check_malloc (2*(n1+n2)*sizeof(float));

for(j=0;j<n2;j++)
  fourg__(&((xc+j*n1)->re),&n1,&isgn,space);

xtc = (struct complex *) _check_malloc (n2*sizeof(struct complex));

for(i=0;i<n1;i++)
   {
   for(j=0;j<n2;j++)
      {
      ip = i + j*n1;

      xtc[j].re = xc[ip].re;
      xtc[j].im = xc[ip].im;
      }

   fourg__(&(xtc->re),&n2,&isgn,space);

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
