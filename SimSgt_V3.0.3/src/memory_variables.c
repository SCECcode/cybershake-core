#include "include.h"

mvar_coefs(mds,medf,nx,ny,nz,dt,nt,fs,vf0,qf0,fmin,fmax,eflg,nodeType,ny1)
struct modelstorage *mds;
float *medf, *dt, *vf0, *qf0, *fmin, *fmax;
int nx, ny, nz, nt, fs, eflg, nodeType, ny1;
{
float *tau, *qpwt, *qswt, *qpfunc, *qsfunc;
float *mdfp[2], *l2mp, *lamp, *mxyp, *mxzp, *myzp, *akp, *pkp, *skp;
float *qp0, *qp1, *qs0, *qs1, qpval, qsval;
float taumin, taumax, ltau0, dtau;
float mubar, denom;
float den2, xre, xim;
float gamma, dfac, uppr, lowr;
float sump, sums, ap, as;
int ix, iy, iz, iyend, izstart, i, j, k, mm;
int ixp, iyp, izp;
int ip, ipx, ipz, ipxz;
int print_tol, ptol, iyinc, iymin, iymax, bigN;
float ptfac;
int spng_zone, constM, xbndflag, ybndflag, zbndflag;

float vw0, qw0, pi2, pre, pim, sre, sim, pmodf, smodf;

float half = 0.5;
float one = 1.0;
float two = 2.0;
float three = 3.0;
float eight = 8.0;
float pi = 3.141592654;

pi2 = two*pi;

izstart = 0;
if(fs)
   izstart = 1;

bigN = 8;

tau = (float *) check_malloc (bigN*sizeof(float));
qpwt = (float *) check_malloc (bigN*sizeof(float));
qswt = (float *) check_malloc (bigN*sizeof(float));

qpfunc = (float *) check_malloc (bigN*sizeof(float));
qsfunc = (float *) check_malloc (bigN*sizeof(float));

if(*fmax < 0.0)
   {
   taumin = (*dt)/pi;   /* inverse of nyquist */
   *fmax = 1.0/(pi2*taumin);
   }
else
   taumin = 1.0/(pi2*(*fmax));

if(*fmin < 0.0)
   {
   taumax = 1.0e+04*taumin;
   *fmin = 1.0/(pi2*taumax);
   }
else
   taumax = 1.0/(pi2*(*fmin));

ltau0 = log(taumin);
dtau = (log(taumax) - log(taumin))/(2.0*bigN);
for(i=0;i<bigN;i++)
   tau[i] = exp(ltau0 + (2*i + 1)*dtau);

qw0 = pi2*(*qf0);
if((*qf0) <= 0.0)
   {
   qw0 = 1.0/sqrt(taumax*taumin);
   *qf0 = qw0/pi2;
   }

vw0 = pi2*(*vf0);
if((*vf0) <= 0.0)
   {
   vw0 = 1.0/sqrt(taumax*taumin);
   *vf0 = vw0/pi2;
   }

fprintf(stderr,"**** vf0= %13.5e qf0= %13.5e fmin= %13.5e fmax= %13.5e\n",*vf0,*qf0,*fmin,*fmax);

init_model_seek(mds,MEDIA_FIELD);
mdfp[0] = medf;
mdfp[1] = medf + N_MED_VARS*nx*nz;

/*
   If ny1 is odd, then periodicity of qpwt, qswt & tau's need
   to be shifted by one grid to be consistent with adjacent
   CPU's.  This is accomplished by incrementing the media
   pointers by one, so that the local iy%2 agrees with the global
   iy%2 (i.e., they are either both even or odd).  The iyend=1 shifts
   the end of the loop back one so that we don't go off the end of
   the model array.
*/

iyend = 0;
if(ny1%2)   /* TRUE if ny1 is odd */
   {
   mdfp[0] = reed_model(mds,mdfp[0],MEDIA_FIELD);
   rite_model(mds,mdfp[0],MEDIA_FIELD);
   iyend = 1;
   }

print_tol = 25;
ptol = print_tol;

iymin = 0;
iymax = ny - iyend;
iyinc = 2;
ptfac = 1.0/(float)(iyinc*(int)((iymax-1-iymin)/iyinc));

for(iy=0;iy<iymax;iy=iy+iyinc)
   {

/*SN:
  not right boundary
*/

   ybndflag = 0;
   if((nodeType == PARL_HEAD) && (iy <= (NBND_PAD+1)))
      ybndflag = 1;
   if((nodeType == PARL_TAIL) && (iy >= ((ny-1-iyend)-(NBND_PAD+1))))
      ybndflag = 1;

   mdfp[0] = reed_model(mds,mdfp[0],MEDIA_FIELD);

   iyp = 1;
   if(iy == ny-1-iyend)
      iyp = 0;

   if(iyp == 1)
      mdfp[1] = reed_model(mds,mdfp[1],MEDIA_FIELD);
   else
      mdfp[1] = mdfp[0];     /* if ny is odd, just re-use last plane */

   qp0 = mdfp[0] + 9*nx*nz;
   qp1 = mdfp[1] + 9*nx*nz;
   qs0 = mdfp[0] + 10*nx*nz;
   qs1 = mdfp[1] + 10*nx*nz;

   for(iz=izstart;iz<nz;iz=iz+2)
      {
      zbndflag = 0;
      if((fs == 0 && iz <= (NBND_PAD+1)) || iz >= ((nz-1)-(NBND_PAD+1)))
         zbndflag = 1;

      for(ix=0;ix<nx;ix=ix+2)
	 {
         xbndflag = 0;
         if(ix <= (NBND_PAD+1) || ix >= ((nx-1)-(NBND_PAD+1)))
            xbndflag = 1;

         ip = ix + iz*nx;

	 ixp = 1;
         if(ix == nx-1)
	    ixp = 0;
	 izp = 1;
         if(iz == nz-1)
	    izp = 0;

	 ipx = ip + ixp;
	 ipz = ip + izp*nx;
	 ipxz = ip + izp*nx + ixp;

/* 
   Compute average Q over 8 grid cells -> use harmonic average (1/Q).  If
   input Q is negative, then sponge zone is true.
*/

         spng_zone = 0;

         qpfunc[0] = qp0[ip];
         qpfunc[1] = qp0[ipx];
         qpfunc[2] = qp0[ipz];
         qpfunc[3] = qp0[ipxz];
         qpfunc[4] = qp1[ip];
         qpfunc[5] = qp1[ipx];
         qpfunc[6] = qp1[ipz];
         qpfunc[7] = qp1[ipxz];

         qsfunc[0] = qs0[ip];
         qsfunc[1] = qs0[ipx];
         qsfunc[2] = qs0[ipz];
         qsfunc[3] = qs0[ipxz];
         qsfunc[4] = qs1[ip];
         qsfunc[5] = qs1[ipx];
         qsfunc[6] = qs1[ipz];
         qsfunc[7] = qs1[ipxz];

         qpval = 0.0;
         qsval = 0.0;

         for(i=0;i<bigN;i++)
            {
            if(qpfunc[i] < 0.0)
               {
               qpfunc[i] = -qpfunc[i];
               spng_zone = 1;
               }
            qpval = qpval + one/qpfunc[i];

            if(qsfunc[i] < 0.0)
               {
               qsfunc[i] = -qsfunc[i];
               spng_zone = 1;
               }
            qsval = qsval + one/qsfunc[i];
            }

         qpval = eight/qpval;
         qsval = eight/qsval;

constM = 0;
/* for constantM formulation
      constM = 1;
      */

/* The following sets qwt coefficients to Steve Day's values
*/

	 ap = 2.*log(taumax/taumin)/(pi*qpval-2.*log(qw0*taumin));
	 as = 2.*log(taumax/taumin)/(pi*qsval-2.*log(qw0*taumin));
         for(i=0;i<bigN;i++)
	    {
            tau[i] = exp(ltau0 + (2*i + 1)*dtau);
            qpwt[i] = ap;
            qswt[i] = as;
	    }

/* more appropriate form for different bandwidths
*/

	 uppr = 0.0;
	 lowr = 0.0;
         for(i=0;i<bigN;i++)
	    {
	    denom = 1.0/(qw0*qw0*tau[i]*tau[i] + 1.0);
	    lowr = lowr + denom/(float)(bigN);
	    uppr = uppr + qw0*tau[i]*denom/(float)(bigN);
	    }

	 ap = 1.0/(lowr + qpval*uppr);
	 as = 1.0/(lowr + qsval*uppr);
         for(i=0;i<bigN;i++)
	    {
            tau[i] = exp(ltau0 + (2*i + 1)*dtau);
            qpwt[i] = ap;
            qswt[i] = as;
	    }

/*
   The following sets qwt coefficients to Robertsson 1994 values with
   one relaxation mechanism -> use this near boundaries to ensure stability
*/

      if(xbndflag || ybndflag || zbndflag)
         {
         for(i=0;i<bigN;i++)
	    {
	    tau[i] = 1.0/qw0;
            qpwt[i] = 2.0/(qpval + 1.0);
            qswt[i] = 2.0/(qsval + 1.0);
	    }
         }

/*
   Check for stability -> qpwt,qswt must be < 1,
		       -> 2*tau > dt.
*/

      for(i=0;i<bigN;i++)
         {
         if(qpwt[i] >= 1.0 || qswt[i] >= 1.0)
	    {
	    fprintf(stderr,"**** ix= %d iy= %d iz= %d\n",ix,iy,iz);
	    fprintf(stderr,"     qpval= %13.5e qsval= %13.5e\n",qpval,qsval);
	    fprintf(stderr,"     qpwt= %13.5e qswt= %13.5e -> unstable, exiting...\n",qpwt[i],qswt[i]);
            exit(-10);
	    }

         if(tau[i] <= 0.5*(*dt))
	    {
	    fprintf(stderr,"**** tau[%d]= %13.5e is too small\n",i,tau[i]);
	    fprintf(stderr,"     must be greater than %13.5e\n",0.5*(*dt));
	    fprintf(stderr,"     possible instability, exiting...\n");
            exit(-10);
	    }
         }

/*
   Use RWG's form for constant Mu using effective M when near boundaries
   in "sponge zone".
*/

	 if(spng_zone)
	    {
            pre = 0.0;
            pim = 0.0;
            sre = 0.0;
            sim = 0.0;
            for(i=0;i<bigN;i++)
               {
               denom = one/(vw0*vw0*tau[i]*tau[i] + 1);
               xre = 1.0 - qpwt[i]*denom;
               xim = qpwt[i]*vw0*tau[i]*denom;

	       den2 = 1.0/(xre*xre + xim*xim);
               pre = pre + xre*den2;
               pim = pim + xim*den2;

               xre = 1.0 - qswt[i]*denom;
               xim = qswt[i]*vw0*tau[i]*denom;

	       den2 = 1.0/(xre*xre + xim*xim);
               sre = sre + xre*den2;
               sim = sim + xim*den2;
               }

	    denom = 1.0/(pre*pre + pim*pim);
	    pre = bigN*pre*denom;
	    pim = bigN*pim*denom;

	    gamma = pim/pre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    pmodf = dfac*0.5*(1.0 + dfac)/pre;

	    denom = 1.0/(sre*sre + sim*sim);
	    sre = bigN*sre*denom;
	    sim = bigN*sim*denom;

	    gamma = sim/sre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    smodf = dfac*0.5*(1.0 + dfac)/sre;
	    }

	 if(constM)  /* Steve Day's formulation for Mu  */
	    {
            pre = 1.0;
            pim = 0.0;
            sre = 1.0;
            sim = 0.0;
            for(i=0;i<bigN;i++)
               {
               denom = one/(bigN*(vw0*vw0*tau[i]*tau[i] + 1));

               pre = pre - qpwt[i]*denom;
               pim = pim + qpwt[i]*vw0*tau[i]*denom;

               sre = sre - qswt[i]*denom;
               sim = sim + qswt[i]*vw0*tau[i]*denom;
               }

	    gamma = pim/pre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    pmodf = dfac*0.5*(1.0 + dfac)/pre;

	    gamma = sim/sre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    smodf = dfac*0.5*(1.0 + dfac)/sre;
	    }

         for(j=0;j<=iyp;j++)
            {
            l2mp = mdfp[j];
            lamp = mdfp[j] +   nx*nz;
            mxyp = mdfp[j] + 2*nx*nz;
            mxzp = mdfp[j] + 3*nx*nz;
            myzp = mdfp[j] + 4*nx*nz;
            akp  = mdfp[j] + 8*nx*nz;
            pkp  = mdfp[j] + 9*nx*nz;
            skp  = mdfp[j] + 10*nx*nz;

            for(k=0;k<=izp;k++)
               {
               for(i=0;i<=ixp;i++)
                  {
		  ip = (ix+i) + (iz+k)*nx;
		  mm = (i%2) + 2*(j%2) + 4*(k%2);

/*
   Use RWG 2001 element specific modulus correction factors for interior
   of model
*/

	          if(spng_zone == 0 && constM == 0)
	             {
                     denom = one/(vw0*vw0*tau[mm]*tau[mm] + 1);

                     pre = 1.0 - qpwt[mm]*denom;
                     pim = qpwt[mm]*vw0*tau[mm]*denom;
                     gamma = pim/pre;
                     dfac = 1.0/(sqrt(1.0 + gamma*gamma));
                     pmodf = dfac*0.5*(1.0 + dfac)/pre;

                     sre = 1.0 - qswt[mm]*denom;
                     sim = qswt[mm]*vw0*tau[mm]*denom;
                     gamma = sim/sre;
                     dfac = 1.0/(sqrt(1.0 + gamma*gamma));
                     smodf = dfac*0.5*(1.0 + dfac)/sre;
		     }

if(pmodf < 0.0 || smodf < 0.0)
   fprintf(stderr,"%d %d) pm=%13.5e sm=%13.5e\n",ix+i,iz+k,pmodf,smodf);

                  mubar = half*(l2mp[ip] - lamp[ip]);
                  l2mp[ip] = pmodf*l2mp[ip];
                  mubar = smodf*mubar;
                  lamp[ip] = l2mp[ip] - two*mubar;

/*
   If eflg=0, then media averaging has not yet been done.  Only one value
   of rigidity is stored (as 1/mu) in the field myzp = mdfp[j] + 4*nx*nz.
   This is the only field that needs to be modified.  Media averaging will
   be done after returning to main().
*/

		  if(eflg == 0)
                     myzp[ip] = myzp[ip]/smodf;
                  else
		     {
                     mxyp[ip] = smodf*mxyp[ip];
                     mxzp[ip] = smodf*mxzp[ip];
                     myzp[ip] = smodf*myzp[ip];
		     }

		  denom = one/(two*tau[mm] + (*dt));

		  akp[ip] = (two*tau[mm] - (*dt))*denom;
		  pkp[ip] = (*dt)*qpwt[mm]*l2mp[ip]*denom;

/*
OLD WAY-> put mubar into sk coefficient.  This creates instability
when sharp velocity contrasts are present in model.
		  skp[ip] = (*dt)*qswt[mm]*mubar*denom;

NEW WAY 07/03/02-> Apply effective mu to sk coefficient in the
routines diff_mv() and diffx_mv().  This seems to prevent instability
and is compatible with the effective mu's that are used for the stress
updates (txy, txz, tyz).
*/

		  skp[ip] = (*dt)*qswt[mm]*denom;
                  }
               }
            }
	 }
      }

   for(j=0;j<=iyp;j++)
      rite_model(mds,mdfp[j],MEDIA_FIELD);

   if((float)(100.0*(iy-iymin))*ptfac >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy-iymin,iymax-1-iymin);
      ptol = ptol + print_tol;
      }
   }

free(tau);
free(qpwt);
free(qswt);
free(qpfunc);
free(qsfunc);

init_model_seek(mds,MEDIA_FIELD);
}

smv_coefs(mds,medf,nx,ny,nz,dt,nt,fs,vf0,qf0,fmin,fmax,eflg)
struct modelstorage *mds;
float *medf, *dt, *vf0, *qf0, *fmin, *fmax;
int nx, ny, nz, nt, fs, eflg;
{
float tau0, qpwt, qswt;
float *mdfp[1], *l2mp, *lamp, *mxyp, *mxzp, *myzp, *akp, *pkp, *skp;
float *qp0, *qs0, qpval, qsval;
float mubar, denom, taumin, taumax;
float xre, xim;
float gamma, dfac, uppr, lowr;
int ix, iy, iz, i, j, k;
int ixp, iyp, izp;
int ip;
int print_tol;

float vw0, qw0, pi2, pre, pim, sre, sim, pmodf, smodf;

float half = 0.5;
float one = 1.0;
float two = 2.0;
float pi = 3.141592654;

pi2 = two*pi;

if(*fmax < 0.0)
   {
   taumin = (*dt)/pi;   /* inverse of nyquist */
   *fmax = 1.0/(pi2*taumin);
   }
else
   taumin = 1.0/(pi2*(*fmax));

if(*fmin < 0.0)
   {
   taumax = 1.0e+04*taumin;
   *fmin = 1.0/(pi2*taumax);
   }
else
   taumax = 1.0/(pi2*(*fmin));

qw0 = pi2*(*qf0);
if((*qf0) <= 0.0)
   {
   qw0 = 1.0/sqrt(taumax*taumin);
   *qf0 = qw0/pi2;
   }

tau0 = 1.0/qw0;

fprintf(stderr,"**** SMV: f0= %13.5e\n",*qf0);

init_model_seek(mds,MEDIA_FIELD);
mdfp[0] = medf;

print_tol = 25;
for(iy=0;iy<ny;iy++)
   {
   if(100.0*(float)(iy+1)/(float)(ny-1) >= (float)(print_tol))
      {
      fprintf(stderr,"\t %3d percent done\n",print_tol);
      print_tol = print_tol + 25;
      }

   mdfp[0] = reed_model(mds,mdfp[0],MEDIA_FIELD);

   l2mp = mdfp[0];
   lamp = mdfp[0] +   nx*nz;
   mxyp = mdfp[0] + 2*nx*nz;
   mxzp = mdfp[0] + 3*nx*nz;
   myzp = mdfp[0] + 4*nx*nz;
   akp  = mdfp[0] + 8*nx*nz;
   pkp  = mdfp[0] + 9*nx*nz;
   skp  = mdfp[0] + 10*nx*nz;

   qp0 = mdfp[0] + 9*nx*nz;
   qs0 = mdfp[0] + 10*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
	 {
         ip = ix + iz*nx;

         qpval = qp0[ip];
         qsval = qs0[ip];

	 if(qpval < 0.0)
	    qpval = -qpval;
	 if(qsval < 0.0)
	    qsval = -qsval;

/*  Rob's SMV weights   */
         qpwt = two/(qpval + one);
         qswt = two/(qsval + one);

         pre = one - half*qpwt;
         pim = half*qpwt;
         gamma = pim/pre;
         dfac = one/(sqrt(one + gamma*gamma));
         pmodf = dfac*half*(one + dfac)/pre;

         sre = one - half*qswt;
         sim = half*qswt;
         gamma = sim/sre;
         dfac = one/(sqrt(one + gamma*gamma));
         smodf = dfac*half*(one + dfac)/sre;

/*  Shawn's (Robertsson's) SMV weights   */
         qpwt = two/(qpval + two);
         qswt = two/(qsval + two);

	 /* strictly Blanch's */
         pmodf = sqrt((two + qpval)/qpval);
         smodf = sqrt((two + qsval)/qsval);

	 /* strictly Shawn's */
         pmodf = one;
         smodf = one;

if(ix == 0 && iy == 0 && iz == 0)
   {
   fprintf(stderr,"**** qpwt= %13.5e qswt= %13.5e\n",qpwt,qswt);
   fprintf(stderr,"     pmod= %13.5e smod= %13.5e\n",pmodf,smodf);
   }

/*
   Check for stability -> qpwt,qswt must be < 1,
		       -> 2*tau > dt.
*/

         if(qpwt >= 1.0 || qswt >= 1.0)
	    {
	    fprintf(stderr,"**** ix= %d iy= %d iz= %d\n",ix,iy,iz);
	    fprintf(stderr,"     qpval= %13.5e qsval= %13.5e\n",qpval,qsval);
	    fprintf(stderr,"     qpwt= %13.5e qswt= %13.5e -> unstable, exiting...\n",qpwt,qswt);
            exit(-10);
	    }

         if(tau0 <= 0.5*(*dt))
	    {
	    fprintf(stderr,"**** tau0= %13.5e is too small\n",tau0);
	    fprintf(stderr,"     must be greater than %13.5e\n",0.5*(*dt));
	    fprintf(stderr,"     possible instability, exiting...\n");
            exit(-10);
	    }

         mubar = half*(l2mp[ip] - lamp[ip]);
         l2mp[ip] = pmodf*l2mp[ip];
         mubar = smodf*mubar;
         lamp[ip] = l2mp[ip] - two*mubar;

/*
   If eflg=0, then media averaging has not yet been done.  Only one value
   of rigidity is stored (as 1/mu) in the field myzp = mdfp[j] + 4*nx*nz.
   This is the only field that needs to be modified.  Media averaging will
   be done after returning to main().
*/

         if(eflg == 0)
            myzp[ip] = myzp[ip]/smodf;
         else
            {
            mxyp[ip] = smodf*mxyp[ip];
            mxzp[ip] = smodf*mxzp[ip];
            myzp[ip] = smodf*myzp[ip];
            }

         denom = one/(two*tau0 + (*dt));

         akp[ip] = (two*tau0 - (*dt))*denom;
         pkp[ip] = (*dt)*qpwt*l2mp[ip]*denom;
         skp[ip] = (*dt)*qswt*denom;
	 }
      }

   rite_model(mds,mdfp[0],MEDIA_FIELD);
   }
init_model_seek(mds,MEDIA_FIELD);
}

void mvar_coefsP3(float *medf,struct runparamsP3 *rpars,float *vf0,float *qf0,float *fmin,float *fmax,int eflg)
{
float *tau, *qpwt, *qswt, *qpfunc, *qsfunc;
float *mdfp[2], *l2mp, *lamp, *mxyp, *mxzp, *myzp, *akp, *pkp, *skp;
float *qp0, *qp1, *qs0, *qs1, qpval, qsval;
float taumin, taumax, ltau0, dtau;
float mubar, denom;
float den2, xre, xim;
float gamma, dfac, uppr, lowr;
float sump, sums, ap, as;
int ix, iy, iz, ixshft, iyshft, izshft, i, j, k, mm;
int ixp, iyp, izp;
int ip, ipx, ipz, ipxz;
int print_tol, ptol, iyinc, iymin, iymax, bigN;
float ptfac;
int spng_zone, constM, xbndflag, ybndflag, zbndflag;

float vw0, qw0, pi2, pre, pim, sre, sim, pmodf, smodf;

float half = 0.5;
float one = 1.0;
float two = 2.0;
float three = 3.0;
float eight = 8.0;
float pi = 3.141592654;

float dt;
int nx, ny, nz, nx1, ny1, nz1;
struct nodeinfo *ni;

ni = &(rpars->ni);
dt = rpars->dt;

nx = ni->loc_nx;
nx1 = ni->nx1;
ny = ni->loc_ny;
ny1 = ni->ny1;
nz = ni->loc_nz;
nz1 = ni->nz1;

pi2 = two*pi;
bigN = 8;

tau = (float *) check_malloc (bigN*sizeof(float));
qpwt = (float *) check_malloc (bigN*sizeof(float));
qswt = (float *) check_malloc (bigN*sizeof(float));

qpfunc = (float *) check_malloc (bigN*sizeof(float));
qsfunc = (float *) check_malloc (bigN*sizeof(float));

if(*fmax < 0.0)
   {
   taumin = dt/pi;   /* inverse of nyquist */
   *fmax = 1.0/(pi2*taumin);
   }
else
   taumin = 1.0/(pi2*(*fmax));

if(*fmin < 0.0)
   {
   taumax = 1.0e+04*taumin;
   *fmin = 1.0/(pi2*taumax);
   }
else
   taumax = 1.0/(pi2*(*fmin));

ltau0 = log(taumin);
dtau = (log(taumax) - log(taumin))/(2.0*bigN);
for(i=0;i<bigN;i++)
   tau[i] = exp(ltau0 + (2*i + 1)*dtau);

qw0 = pi2*(*qf0);
if((*qf0) <= 0.0)
   {
   qw0 = 1.0/sqrt(taumax*taumin);
   *qf0 = qw0/pi2;
   }

vw0 = pi2*(*vf0);
if((*vf0) <= 0.0)
   {
   vw0 = 1.0/sqrt(taumax*taumin);
   *vf0 = vw0/pi2;
   }

fprintf(stderr,"**** vf0= %13.5e qf0= %13.5e fmin= %13.5e fmax= %13.5e\n",*vf0,*qf0,*fmin,*fmax);

/*
   If nx1, ny1, nz1 is odd then periodicity of qpwt, qswt & tau's need
   to be shifted by one grid along each respective axis to be consistent with adjacent
   CPU's.  This is accomplished by incrementing the media
   pointers by one, so that the local ix%2, iy%2, iz%2 agrees with the global
   ix%2, iy%2, iz%2, respectively (i.e., they are either both even or odd).
   The ixshft=1, iyshft=1, izshft=1 shifts the index by one.

   For nz1, this criteria is reversed if the top of the model has a free-surfurace
   because the entire model is then shifted down by 1 grid.
*/

ixshft = 0;
if(nx1%2)   /* TRUE if nx1 is odd */
   ixshft = 1;

iyshft = 0;
if(ny1%2)   /* TRUE if ny1 is odd */
   iyshft = 1;

izshft = 0;
if((nz1%2 && rpars->freesurf==0) || ((nz1%2)==0 && rpars->freesurf))
   izshft = 1;

print_tol = 25;
ptol = print_tol;

iyinc = 2;
ptfac = 1.0/(float)(iyinc*(int)(((ny - iyshft)-1)/iyinc));

for(iy=iyshft;iy<ny;iy=iy+iyinc)
   {
   ybndflag = 0;
   if((ni->minusId_y == -1) && (iy <= (NBND_PAD+1)))
      ybndflag = 1;
   if((ni->plusId_y == -1) && (iy >= ((ny-1)-(NBND_PAD+1))))
      ybndflag = 1;

   iyp = 1;
   if(iy == ny-1)
      iyp = 0;

   mdfp[0] = medf + iy*N_MED_VARS*nx*nz;

   if(iyp == 1)
      mdfp[1] = medf + (iy+1)*N_MED_VARS*nx*nz;
   else
      mdfp[1] = mdfp[0];     /* if ny is odd, just re-use last plane */

   qp0 = mdfp[0] + 9*nx*nz;
   qp1 = mdfp[1] + 9*nx*nz;
   qs0 = mdfp[0] + 10*nx*nz;
   qs1 = mdfp[1] + 10*nx*nz;

   for(iz=izshft;iz<nz;iz=iz+2)
      {
      zbndflag = 0;
      if((ni->minusId_z == -1) && (rpars->freesurf == 0) && (iz <= (NBND_PAD+1)))
         zbndflag = 1;
      if((ni->plusId_z == -1) && (iz >= ((nz-1)-(NBND_PAD+1))))
         zbndflag = 1;

      izp = 1;
      if(iz == nz-1)
         izp = 0;

      for(ix=ixshft;ix<nx;ix=ix+2)
	 {
         xbndflag = 0;
         if((ni->minusId_x == -1) && (ix <= (NBND_PAD+1)))
            xbndflag = 1;
         if((ni->plusId_x == -1) && (ix >= ((nx-1)-(NBND_PAD+1))))
            xbndflag = 1;
	 ixp = 1;
         if(ix == nx-1)
	    ixp = 0;

         ip = ix + iz*nx;
	 ipx = ip + ixp;
	 ipz = ip + izp*nx;
	 ipxz = ip + izp*nx + ixp;

/* 
   Compute average Q over 8 grid cells -> use harmonic average (1/Q).  If
   input Q is negative, then sponge zone is true.
*/

         spng_zone = 0;

         qpfunc[0] = qp0[ip];
         qpfunc[1] = qp0[ipx];
         qpfunc[2] = qp0[ipz];
         qpfunc[3] = qp0[ipxz];
         qpfunc[4] = qp1[ip];
         qpfunc[5] = qp1[ipx];
         qpfunc[6] = qp1[ipz];
         qpfunc[7] = qp1[ipxz];

         qsfunc[0] = qs0[ip];
         qsfunc[1] = qs0[ipx];
         qsfunc[2] = qs0[ipz];
         qsfunc[3] = qs0[ipxz];
         qsfunc[4] = qs1[ip];
         qsfunc[5] = qs1[ipx];
         qsfunc[6] = qs1[ipz];
         qsfunc[7] = qs1[ipxz];

         qpval = 0.0;
         qsval = 0.0;

         for(i=0;i<bigN;i++)
            {
            if(qpfunc[i] < 0.0)
               {
               qpfunc[i] = -qpfunc[i];
               spng_zone = 1;
               }
            qpval = qpval + one/qpfunc[i];

            if(qsfunc[i] < 0.0)
               {
               qsfunc[i] = -qsfunc[i];
               spng_zone = 1;
               }
            qsval = qsval + one/qsfunc[i];
            }

         qpval = eight/qpval;
         qsval = eight/qsval;

constM = 0;
/* for constantM formulation
      constM = 1;
      */

/* The following sets qwt coefficients to Steve Day's values
*/

	 ap = 2.*log(taumax/taumin)/(pi*qpval-2.*log(qw0*taumin));
	 as = 2.*log(taumax/taumin)/(pi*qsval-2.*log(qw0*taumin));
         for(i=0;i<bigN;i++)
	    {
            tau[i] = exp(ltau0 + (2*i + 1)*dtau);
            qpwt[i] = ap;
            qswt[i] = as;
	    }

/* more appropriate form for different bandwidths
*/

	 uppr = 0.0;
	 lowr = 0.0;
         for(i=0;i<bigN;i++)
	    {
	    denom = 1.0/(qw0*qw0*tau[i]*tau[i] + 1.0);
	    lowr = lowr + denom/(float)(bigN);
	    uppr = uppr + qw0*tau[i]*denom/(float)(bigN);
	    }

	 ap = 1.0/(lowr + qpval*uppr);
	 as = 1.0/(lowr + qsval*uppr);
         for(i=0;i<bigN;i++)
	    {
            tau[i] = exp(ltau0 + (2*i + 1)*dtau);
            qpwt[i] = ap;
            qswt[i] = as;
	    }

/*
   The following sets qwt coefficients to Robertsson 1994 values with
   one relaxation mechanism -> use this near boundaries to ensure stability
*/

      if(xbndflag || ybndflag || zbndflag)
         {
         for(i=0;i<bigN;i++)
	    {
	    tau[i] = 1.0/qw0;
            qpwt[i] = 2.0/(qpval + 1.0);
            qswt[i] = 2.0/(qsval + 1.0);
	    }
         }

/*
   Check for stability -> qpwt,qswt must be < 1,
		       -> 2*tau > dt.
*/

      for(i=0;i<bigN;i++)
         {
         if(qpwt[i] >= 1.0 || qswt[i] >= 1.0)
	    {
	    fprintf(stderr,"**** ix= %d iy= %d iz= %d\n",ix,iy,iz);
	    fprintf(stderr,"     qpval= %13.5e qsval= %13.5e\n",qpval,qsval);
	    fprintf(stderr,"     qpwt= %13.5e qswt= %13.5e -> unstable, exiting...\n",qpwt[i],qswt[i]);
            exit(-10);
	    }

         if(tau[i] <= 0.5*dt)
	    {
	    fprintf(stderr,"**** tau[%d]= %13.5e is too small\n",i,tau[i]);
	    fprintf(stderr,"     must be greater than %13.5e\n",0.5*dt);
	    fprintf(stderr,"     possible instability, exiting...\n");
            exit(-10);
	    }
         }

/*
   Use RWG's form for constant Mu using effective M when near boundaries
   in "sponge zone".
*/

	 if(spng_zone)
	    {
            pre = 0.0;
            pim = 0.0;
            sre = 0.0;
            sim = 0.0;
            for(i=0;i<bigN;i++)
               {
               denom = one/(vw0*vw0*tau[i]*tau[i] + 1);
               xre = 1.0 - qpwt[i]*denom;
               xim = qpwt[i]*vw0*tau[i]*denom;

	       den2 = 1.0/(xre*xre + xim*xim);
               pre = pre + xre*den2;
               pim = pim + xim*den2;

               xre = 1.0 - qswt[i]*denom;
               xim = qswt[i]*vw0*tau[i]*denom;

	       den2 = 1.0/(xre*xre + xim*xim);
               sre = sre + xre*den2;
               sim = sim + xim*den2;
               }

	    denom = 1.0/(pre*pre + pim*pim);
	    pre = bigN*pre*denom;
	    pim = bigN*pim*denom;

	    gamma = pim/pre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    pmodf = dfac*0.5*(1.0 + dfac)/pre;

	    denom = 1.0/(sre*sre + sim*sim);
	    sre = bigN*sre*denom;
	    sim = bigN*sim*denom;

	    gamma = sim/sre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    smodf = dfac*0.5*(1.0 + dfac)/sre;
	    }

	 if(constM)  /* Steve Day's formulation for Mu  */
	    {
            pre = 1.0;
            pim = 0.0;
            sre = 1.0;
            sim = 0.0;
            for(i=0;i<bigN;i++)
               {
               denom = one/(bigN*(vw0*vw0*tau[i]*tau[i] + 1));

               pre = pre - qpwt[i]*denom;
               pim = pim + qpwt[i]*vw0*tau[i]*denom;

               sre = sre - qswt[i]*denom;
               sim = sim + qswt[i]*vw0*tau[i]*denom;
               }

	    gamma = pim/pre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    pmodf = dfac*0.5*(1.0 + dfac)/pre;

	    gamma = sim/sre;
	    dfac = 1.0/(sqrt(1.0 + gamma*gamma));
	    smodf = dfac*0.5*(1.0 + dfac)/sre;
	    }

         for(j=0;j<=iyp;j++)
            {
            l2mp = mdfp[j];
            lamp = mdfp[j] +   nx*nz;
            mxyp = mdfp[j] + 2*nx*nz;
            mxzp = mdfp[j] + 3*nx*nz;
            myzp = mdfp[j] + 4*nx*nz;
            akp  = mdfp[j] + 8*nx*nz;
            pkp  = mdfp[j] + 9*nx*nz;
            skp  = mdfp[j] + 10*nx*nz;

            for(k=0;k<=izp;k++)
               {
               for(i=0;i<=ixp;i++)
                  {
		  ip = (ix+i) + (iz+k)*nx;
		  mm = (i%2) + 2*(j%2) + 4*(k%2);

/*
   Use RWG 2001 element specific modulus correction factors for interior
   of model
*/

	          if(spng_zone == 0 && constM == 0)
	             {
                     denom = one/(vw0*vw0*tau[mm]*tau[mm] + 1);

                     pre = 1.0 - qpwt[mm]*denom;
                     pim = qpwt[mm]*vw0*tau[mm]*denom;
                     gamma = pim/pre;
                     dfac = 1.0/(sqrt(1.0 + gamma*gamma));
                     pmodf = dfac*0.5*(1.0 + dfac)/pre;

                     sre = 1.0 - qswt[mm]*denom;
                     sim = qswt[mm]*vw0*tau[mm]*denom;
                     gamma = sim/sre;
                     dfac = 1.0/(sqrt(1.0 + gamma*gamma));
                     smodf = dfac*0.5*(1.0 + dfac)/sre;
		     }

if(pmodf < 0.0 || smodf < 0.0)
   fprintf(stderr,"%d %d) pm=%13.5e sm=%13.5e\n",ix+i,iz+k,pmodf,smodf);

                  mubar = half*(l2mp[ip] - lamp[ip]);
                  l2mp[ip] = pmodf*l2mp[ip];
                  mubar = smodf*mubar;
                  lamp[ip] = l2mp[ip] - two*mubar;

/*
   If eflg=0, then media averaging has not yet been done.  Only one value
   of rigidity is stored (as 1/mu) in the field myzp = mdfp[j] + 4*nx*nz.
   This is the only field that needs to be modified.  Media averaging will
   be done after returning to main().
*/

		  if(eflg == 0)
                     myzp[ip] = myzp[ip]/smodf;
                  else
		     {
                     mxyp[ip] = smodf*mxyp[ip];
                     mxzp[ip] = smodf*mxzp[ip];
                     myzp[ip] = smodf*myzp[ip];
		     }

		  denom = one/(two*tau[mm] + dt);

		  akp[ip] = (two*tau[mm] - dt)*denom;
		  pkp[ip] = dt*qpwt[mm]*l2mp[ip]*denom;

		  skp[ip] = dt*qswt[mm]*denom;
                  }
               }
            }
	 }
      }

   if((float)(100.0*(iy-iyshft))*ptfac >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy-iyshft,(ny - iyshft)-1);
      ptol = ptol + print_tol;
      }
   }

free(tau);
free(qpwt);
free(qswt);
free(qpfunc);
free(qsfunc);
}
