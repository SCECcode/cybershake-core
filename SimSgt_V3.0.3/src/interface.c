/*
   interface_cfac.c:

   GF response is interpolated at the central location of the grid indicies
   of vx, vy, and vz.  Thus only a single array is needed for cfac[] weights,
   and only one  GF start time (itgfst) is required.  The location of the
   central point is given as (ix+h/4,iy+h/4,iz+h/4).  These parameters are
   stored in the structure "gfindexpar".  See the file 'structure_cfac.h'.
*/

/*
   interface.c contains the following functions:

      add_infc_bnd()
      calc_infc_pnt()
      calc_infc_pnt1D()
      get_infc_par()
      init_pointsrc()
      norm_stf()
      partition()
      sub_infc_bnd()
*/

#include "include.h"

/*
   sub_infc_bnd() removes the interfacing wave from the grid points along
   the outer and inner edges of the model boundaries.

      v1  = velocity field for previous time step along outer boundary edge,
          -> stored on xz plane
      v0  = velocity field for current time step along interior boundary edge,
          -> stored on xz plane
      vir = velocity field for previous time step along interior boundary edge,
          -> stored on xz plane when iflag=YZERO,YN
             stored in x and z strips otherwise [see store_iring() in absorb.c]
*/

sub_infc_bnd(nx,nz,vx1,vx0,vxir,vy1,vy0,vyir,vz1,vz0,vzir,kflag,med,intfac,iflag)
struct interface *intfac;
int nx, nz, kflag, iflag;
float *vx0, *vx1, *vxir, *vy0, *vy1, *vyir, *vz0, *vz1, *vzir, *med;
{
switch(intfac->mode)
   {
   case FAULT:
      get_infc_fault(-1,nx,nz,vx1,vx0,vxir,vy1,vy0,vyir,vz1,vz0,vzir,kflag,med,intfac,iflag);
      break;
   }
}

/*
   add_infc_bnd() adds in the interfacing wave to the grid points along
   the outer and inner edges of the model boundaries.

      v1 = velocity field for current time step along outer boundary edge,
         -> stored on xz plane
      v0 = velocity field for current time step along interior boundary edge,
         -> stored on xz plane
*/

add_infc_bnd(nx,nz,vx1,vx0,vy1,vy0,vz1,vz0,kflag,med,intfac,iflag)
struct interface *intfac;
int nx, nz, kflag, iflag;
float *vx0, *vx1, *vy0, *vy1, *vz0, *vz1, *med;
{
int go = 0;

switch(intfac->mode)
   {
   case FAULT:
      get_infc_fault(1,nx,nz,vx1,vx0,vx0,vy1,vy0,vy0,vz1,vz0,vz0,kflag,med,intfac,iflag);
      break;
   }
return(go);
}

/*
   get_infc_fault() retrieves the values of the interfacing wave field along
   the outer and inner edges of the model boundaries and adds or subtracts
   these values from the FD wave field according to the value of sgn.
*/

get_infc_fault(sgn,nx,nz,vx1,vx0,vxir,vy1,vy0,vyir,vz1,vz0,vzir,kflag,med,intfac,rflag)
struct interface *intfac;
int nx, nz, kflag, sgn, rflag;
float *vx0, *vx1, *vxir, *vy0, *vy1, *vyir, *vz0, *vz1, *vzir, *med;
{
struct gfsourcepar *gfsrcp;
struct gfindexpar *gfip;
struct gftable *gftab;
float *vp, *vpi, *vpx, *vpxi, *vpy, *vpyi, *vpz, *vpzi;
float av;
int it, ic, iyspt, ip, j, jstart, jend, ibnd;
int *nsrc, isrc, ofst, *indx1, *indx2, nbyte;
int xndx, yndx, zndx, itspt, itspt1, irngflag;
int ixs, ixe, izs, ize, xinc, zinc, ix, iz;
int go, *gopln, *nwgo, fdr, fdr_int, fdr_ext;
int k, vxshft, vzshft;
char ychar[8], fname_int[256], fname_ext[256], *bufptr;

/*
int *igf;
*/
unsigned short *jsrc, *igf, govflag;
float *cfac;

int size_int = sizeof(int);
int size_gfsrcp = sizeof(struct gfsourcepar);
int size_gfip = sizeof(struct gfindexpar);

int ibndst = 0;

iyspt = intfac->iy - intfac->yshft;

/*
   First, check to see if this plane expects a non-zero response for the
   interfacing wave field at this time step. This is only done if sgn>0
   [call from add_infc_pnt()]:
      
      -If no response is expected (intfac->goplane[iyspt]=0), return without
       computing anything.

      -If a non-zero response is expected along this plane of the model
       (intfac->goplane[iyspt]=1), then the interfacing wave field is
       read in from the temporary files.  If the interfacing wave field is
       active (non-zero), reset intfac->newgo[iyspt]=1 before returning
       from function.  If the interfacing wave field is no longer active,
       function returns with intfac->newgo[iyspt]=0.
*/

gopln = intfac->goplane;
nwgo = intfac->newgo;

if(!gopln[iyspt])
   go = 0;
else
   go = 1;

/*
ibndst = 1;
*/

if(go)
   {
   gftab = intfac->gftable;

   if(iyspt < 10)
      sprintf(ychar,"00%1d",iyspt);
   else if(iyspt < 100)
      sprintf(ychar,"0%2d",iyspt);
   else if(iyspt < 1000)
      sprintf(ychar,"%3d",iyspt);

   for(ibnd=ibndst;ibnd<2;ibnd++)
      {
      irngflag = 0;

      switch(ibnd)
         {
/*
   case 0: interior boundary edge - always at current time step for v0,
				  - previous time step for vir when sgn=-1
   case 1: exterior boundary edge - if sgn=-1 (sub), use previous time step
		                  - if sgn=1 (add), use current time step
*/
         case 0:
            sprintf(fname_int,"%sint.%3s",intfac->gfname,ychar);
	    vpx = vx0;
	    vpy = vy0;
	    vpz = vz0;
	    ip = 1;
            it = intfac->it;

            if(sgn < 0)
               irngflag = 1;
	    break;
         case 1:
            sprintf(fname_ext,"%sext.%3s",intfac->gfname,ychar);
	    vpx = vx1;
	    vpy = vy1;
	    vpz = vz1;
	    ip = 0;
            it = intfac->it;
            if(it && sgn < 0)
               it = intfac->it - 1;
	    break;
         }

      if(iyspt == 0 || iyspt == intfac->ny-1)
         {
         jstart = 4;
         jend = 5;

	 if(ibnd == 0)
            fdr = opfile(fname_int);
	 else
            fdr = opfile(fname_ext);
         }
      else
         {
         jstart = 0;
         jend = 4;

	   /* read all gf information in one block when iy != 0,ny-1 */

	 if(rflag == OPEN_READ)
	    {
	    if(ibnd == 0)
	       {
               fdr_int = opfile(fname_int);
               reed(fdr_int,&nbyte,size_int);
               reed(fdr_int,intfac->buf_int,nbyte);
	       close(fdr_int);
	       }
	    else
	       {
               fdr_ext = opfile(fname_ext);
               reed(fdr_ext,&nbyte,size_int);
               reed(fdr_ext,intfac->buf_ext,nbyte);
	       close(fdr_ext);
	       }
	    }

	 if(ibnd == 0)
            bufptr = intfac->buf_int;
	 else
            bufptr = intfac->buf_ext;
         }

      for(j=jstart;j<jend;j++)
	 {
	 vxshft = vzshft = 0;

         switch(j)
            {
/*
   case 0: iz=0         boundary
   case 1: ix=nx-1      boundary
   case 2: iz=nz-1      boundary
   case 3: ix=0         boundary
   case 4: iy=0,ny-1    boundary
*/
            case 0:
	       ofst = 2;
               ixs = 1 + ip;      ixe = nx - ip; xinc = 1;
               izs = ip;          ize = 1 + ip;  zinc = 1;
	       break;
            case 1:
	       ofst = nx + 2;
               ixs = nx - 1 - ip; ixe = nx - ip; xinc = 1;
               izs = 1 + ip;      ize = nz - ip; zinc = 1;

	       vxshft = -1;
	       break;
            case 2:
	       ofst = nx + nz + 2;
               ixs = nx - 2 - ip; ixe = ip - 1;  xinc = -1;
               izs = nz - 1 - ip; ize = nz - ip; zinc = 1;

	       vzshft = -1;
	       break;
            case 3:
	       ofst = 2*nx + nz + 2;
               ixs = ip;          ixe = 1 + ip;  xinc = 1;
               izs = nz - 2 - ip; ize = ip - 1;  zinc = -1;
	       break;
            case 4:
	       ofst = 0;
               ixs = ip;          ixe = nx - ip; xinc = 1;
               izs = ip;          ize = nz - ip; zinc = 1;
	       break;
            }

/*
   set-up parameters for subtracting interfacing wave field from inner
   ring variables at previous time step
*/
	 if(irngflag)
	    {
            vpxi = vxir + ofst;
            vpyi = vyir + ofst;
            vpzi = vzir + ofst;
	    }

	 ic = 0; /* counter for inner ring arrays */
         for(iz=izs;iz!=ize;iz=iz+zinc)
            {

	    if(j == 4) /* read gf information in x-strips when iy = 0,ny-1 */
	       {
               reed(fdr,&nbyte,size_int);
               reed(fdr,intfac->buf_ext,nbyte);
               bufptr = intfac->buf_ext;
	       }

            for(ix=ixs;ix!=ixe;ix=ix+xinc)
               {
	       xndx = ix+vxshft + iz*nx;
	       yndx = ix + iz*nx;
	       zndx = ix + (iz+vzshft)*nx;

	       nsrc = (int *) bufptr;
	       bufptr = bufptr + size_int;

               if(*nsrc)
		  {
		  gfip = (struct gfindexpar *) bufptr;
	          bufptr = bufptr + size_gfip;

                  for(isrc=0;isrc<(*nsrc);isrc++)
                     {
		     gfsrcp = (struct gfsourcepar *) bufptr;
	             bufptr = bufptr + size_gfsrcp;

		     for(k=0;k<3;k=k+kflag)
		        {
		        switch(k)
			   {
                           case 0:    /* VX update */
			      vp = vpx;
			      vpi = vpxi;
                              govflag = gfsrcp->vflag & ONE_MASK;

                              indx1 = &xndx;
	                      break;
                           case 1:    /* VY update */
			      vp = vpy;
			      vpi = vpyi;
                              govflag = (gfsrcp->vflag >> 1) & ONE_MASK;

                              indx1 = &yndx;
	                      break;
                           case 2:    /* VZ update */
			      vp = vpz;
                              vpi = vpzi;
                              govflag = (gfsrcp->vflag >> 2) & ONE_MASK;

                              indx1 = &zndx;
	                      break;
			   }

                        itspt = it - gfsrcp->itgfst;

                        if(itspt < 0)  /* response hasn't arrived yet */
                           nwgo[iyspt] = 1;

                        else if(itspt <= (int)(gfsrcp->npts) && govflag)
		           {
                           nwgo[iyspt] = 1;
			   jsrc = &gfsrcp->isrc;
			   igf = &gfip->igf[0];
                           cfac = &gfip->cfac[0];

			   if(itspt < (int)(gfsrcp->npts))
                              {
                              get_gfamp(gftab,jsrc,itspt,k,igf,cfac,&av);
                              vp[*indx1] = vp[*indx1] + sgn*av;
                              }

			   itspt1 = itspt - 1;
			   if(irngflag && itspt1 >= 0)
                              {
		              indx2 = indx1;
                              if(j < 4)
                                 indx2 = &ic;

                              get_gfamp(gftab,jsrc,itspt1,k,igf,cfac,&av);

                              vpi[*indx2] = vpi[*indx2] - av;
                              }
		           }
		        }
		     }
		  }
	       ic++;  /* increment counter for inner ring arrays */
	       }
	    }
	 }

      if(iyspt != 0 && iyspt != intfac->ny-1)
	 {
	 if(nwgo[iyspt])    /* fill in corners of interior ring values */
	    {
            if(irngflag)
	       {
	       vxir[1] = vxir[2*(nx+nz)-2];
	       vxir[nx+1] = vxir[nx-2-1];
	       vxir[nx+nz+1] = vxir[nx+nz-2];
	       vxir[2*nx+nz+1] = vxir[2*nx+nz-2];

	       if(kflag == 1)
		  {
	          vyir[1] = vyir[2*(nx+nz)-2];
	          vyir[nx+1] = vyir[nx-2];
	          vyir[nx+nz+1] = vyir[nx+nz-2];
	          vyir[2*nx+nz+1] = vyir[2*nx+nz-2];
		  }

	       vzir[1] = vzir[2*(nx+nz)-2];
	       vzir[nx+1] = vzir[nx-2];
	       vzir[nx+nz+1] = vzir[nx+nz-2-1];
	       vzir[2*nx+nz+1] = vzir[2*nx+nz-2];
	       }
	    }
	 }
      else
         close(fdr);
      }
   }
}

/*
   partition_ff() calculates the analytic response at the current grid
   location using a ray tracing algorithm for individual subfaults of
   a finite-fault

      *** far-field S-wave response for a point double-couple source

   The response is partitioned into vx, vy, and vz components.
*/

struct vector
   {
   float x;
   float y;
   float z;
   };

partition_ff(ff,ps,srcp,gftab,isrc,gf,itend,beta)
struct finitefault *ff;
struct pointsource *ps;
struct sourceparams *srcp;
struct gftable *gftab;
float *beta, *gf;
int itend, isrc;
{
float *gfvx, *gfvy, *gfvz;
float xsrc, ysrc, zsrc, rng, rayp;
float ttime, rad, fsv, fsh, samp;
float xs, ys, zs, invrad;
float cosI, sinI, cos2I, sin2I;
float cosP0, sinP0, cos2P0, sin2P0;
float cosP, sinP;
float invdt, src_const;
float qfac, qconst;
int ix, iy, iz, nstf2;
int i, j, it, its, ite;

float fac;
int it0;

struct vector uv, vv, iv, pv, tv;
float at, bt, ct, dt;

float half = 0.5;
float one = 1.0;
float two = 2.0;
float four = 4.0;
float eps = 1.0e-06;

float transh = 1.0;
float transv = 1.0;

int nstf = 8;  /* time length (T=nstf*dt) of point source time function */
nstf = 8;

ix = gftab->ix;
iy = gftab->iy;
iz = gftab->iz;

qconst = -PI*0.5/200.0;
qconst = 0.0;

rayp = 0.5*eps;
invdt = one/(ps->dt);

src_const = two*invdt/nstf; /* stf normalized amplitude */
src_const = four*invdt*invdt/(nstf*nstf); /* stf normalized amplitude */

xs = ix*ps->h;
ys = iy*ps->h;
zs = iz*ps->h;

for(j=0;j<ff->npnt2;j++)
   {
   samp = src_const*ps->srcamp[j]; /* source amplitude in homogeneous region */

   xsrc = xs - ps->xs[j];
   ysrc = ys - ps->ys[j];
   zsrc = zs - ps->zs[j];

   rng = sqrt(xsrc*xsrc + ysrc*ysrc);

/*
   do earthquake-> Aki & Richards, pg 115

   first, calculate amplitude:
   find ray parameter, ray length, and travel time for source location
*/

   get_rayp(ff,&ps->h,&ps->zs[j],&zs,&rng,&rayp,&rad,&ttime,&eps);

   if(rayp >= 0.0) /* layered media */
      {
      invrad = one/rad;

      sinI = rayp*ps->vsrc[j];        /* I is take-off angle at source */
      if(zsrc < 0.0)
         cosI = -sqrt(one - sinI*sinI);
      else
         cosI = sqrt(one - sinI*sinI);

/*
   Get SS transmission coefficients for ray propagating through layered media.
   These are multiplied together and stored in transh, transv.
*/

      get_transmission(ff,&ps->zs[j],&zs,&rayp,&transh,&transv);
      }
   else /* homogeneous media */
      {
      rad = sqrt(xsrc*xsrc + ysrc*ysrc + zsrc*zsrc);
      ttime = rad/ps->vsrc[j];
      invrad = one/rad;

      sinI = rng*invrad;
      cosI = zsrc*invrad;
      }

   if(rng > 0.0)
      {
      cosP = -ysrc/rng;
      sinP = xsrc/rng;
      }
   else
      {
      cosP = one;
      sinP = 0.0;
      }

/* compute Q factor for anelastic attenuation */

   qfac = exp(qconst*rad/(ps->vsrc[j]));

/* compute SV and SH radiation pattern coefficients */

/* XXXXXXXXXXXXXXXXXXXXXXX

uv.x = ps->cosL*ps->sinA - ps->cosD*ps->sinL*ps->cosA;
uv.y = -(ps->cosL*ps->cosA + ps->cosD*ps->sinL*ps->sinA);
uv.z = -ps->sinL*ps->sinD;

vv.x = ps->sinD*ps->cosA;
vv.y = ps->sinD*ps->sinA;
vv.z = -ps->cosD;

iv.x = sinI*sinP;
iv.y = -sinI*cosP;
iv.z = cosI;

pv.x = cosI*sinP;
pv.y = -cosI*cosP;
pv.z = -sinI;

tv.x = cosP;
tv.y = sinP;
tv.z = 0.0;

dotp(&at,&iv,&vv);
dotp(&bt,&uv,&pv);
dotp(&ct,&iv,&uv);
dotp(&dt,&vv,&pv);

fsv = at*bt + ct*dt;
*/

/* XXXXXXXXXXXXXXXXXXXXXXX */

   cos2I = cosI*cosI - sinI*sinI;
   sin2I = two*sinI*cosI;

   cosP0 = cosP*ps->cosA + sinP*ps->sinA;
   sinP0 = sinP*ps->cosA - cosP*ps->sinA;

   cos2P0 = cosP0*cosP0 - sinP0*sinP0;
   sin2P0 = two*sinP0*cosP0;

   fsv = ps->sinL*ps->cos2D*cos2I*sinP0
             - ps->cosL*ps->cosD*cos2I*cosP0
             + half*ps->cosL*ps->sinD*sin2I*sin2P0
             - half*ps->sinL*ps->sin2D*sin2I*(one + sinP0*sinP0);

/*
if(transv > transh)
   fprintf(stderr,"SH=%10.3f  SV=%10.3f\n",transh,transv);
*/

   fsv = fsv*transv;

/*
   fsv = fsv*transv;

   07/15/94:  SV transmission coefficient seems to be very small.  Need to
   check accuracy of computation formulas in get_transmission().  For now,
   just use SH transmission coefficient for SV.
*/

   fsh = ps->cosL*ps->cosD*cosI*sinP0
             + ps->cosL*ps->sinD*sinI*cos2P0
             + ps->sinL*ps->cos2D*cosI*cosP0
             - half*ps->sinL*ps->sin2D*sinI*sin2P0;

   fsh = fsh*transh;

/*
   Next, calculate tstart for this point source.  Also set gftab[isrc]->its
   to be the earliest start time for this subfault.
*/

/* NEW TSMIN CALCULATION
   srcp[j].ts = (ps->rtime[j] + ttime) - (ps->itshft)*(ps->dt);
*/
   srcp[j].ts = (ps->rtime[j] + ttime);

   if(srcp[j].ts < ff->tsmin)
      ff->tsmin = srcp[j].ts;

   if(srcp[j].ts < gftab->ts[isrc])
      gftab->ts[isrc] = srcp[j].ts;

/*
   Partition into x, y, and z velocity components.  The amplitude
   used in the calculations is further modified by moment scaling,
   impedance, radiation pattern, and subfault weight, which have all
   been precomputed in the variable 'samp'. In addition, I is changed
   to incidence angle at grid location.
*/

   sinI = sinI*(*beta)/ps->vsrc[j];
   if(zsrc < 0.0)
      cosI = -sqrt(one - sinI*sinI);
   else
      cosI = sqrt(one - sinI*sinI);

   srcp[j].vxamp = samp*(cosI*sinP*fsv + cosP*fsh)*invrad*qfac;
   srcp[j].vyamp = samp*(-cosI*cosP*fsv + sinP*fsh)*invrad*qfac;
   srcp[j].vzamp = samp*(-sinI*fsv)*invrad*qfac;
   }

/*
   Use srcp array to build up Green's function.  Do 2-sided boxcar for
   'nstf' timesteps.
   
   While source pulse is active, the shape of the velocity
   waveform is given by the time derivative of the far-field displacement
   pulse, which we have assumed to be a triangle of unit area.  The velocity
   waveform is thus a two-sided boxcar with amplitude A given by:
 
                   A =  4/(T*T);        0    < it <= nstf/2
                   A = -4/(T*T);      nstf/2 < it <= nstf
 
   where T = nstf*dt is the time width of the pulse.  The amplitude
   scaling has already been applied to each point source using the variable
   'samp' (also 'src_const'), above.
*/

gfvx = gf;
gfvy = gf + itend;
gfvz = gf + 2*itend;
it0 = 1;
ite = -1;
nstf2 = nstf/2;
for(j=0;j<ff->npnt2;j++)
   {
   its = (srcp[j].ts - gftab->ts[isrc])*invdt;

/*
fprintf(stderr,"isrc=%d j=%d its=%d\n",isrc,j,its);
*/

   /*     boxcar    */
   for(it=0;it<nstf2;it++)
      {
      i = it + its + it0;
      gfvx[i] = gfvx[i] + srcp[j].vxamp;
      gfvy[i] = gfvy[i] + srcp[j].vyamp;
      gfvz[i] = gfvz[i] + srcp[j].vzamp;
      }
   for(it=nstf2;it<nstf;it++)
      {
      i = it + its + it0;
      gfvx[i] = gfvx[i] - srcp[j].vxamp;
      gfvy[i] = gfvy[i] - srcp[j].vyamp;
      gfvz[i] = gfvz[i] - srcp[j].vzamp;
      }

   /*     triangle
   for(it=0;it<=nstf2;it++)
      {
      i = it + its + it0;
      fac = (float)(it-it0)/(float)(nstf2);

      gfvx[i] = gfvx[i] + fac*srcp[j].vxamp;
      gfvy[i] = gfvy[i] + fac*srcp[j].vyamp;
      gfvz[i] = gfvz[i] + fac*srcp[j].vzamp;
      }
   for(it=nstf2+1;it<nstf;it++)
      {
      i = it + its + it0;
      fac = (float)(nstf+it0-it-1)/(float)(nstf2);

      gfvx[i] = gfvx[i] + fac*srcp[j].vxamp;
      gfvy[i] = gfvy[i] + fac*srcp[j].vyamp;
      gfvz[i] = gfvz[i] + fac*srcp[j].vzamp;
      }
      */

   if(i > ite)
      ite = i;
   }

gftab->npts[isrc] = ite + 1;
}

get_infc_par(intfac)
struct interface *intfac;
{
char srctype[32];

mstpar("srctype","s",srctype);
mstpar("gfdz","f",&intfac->dz);
mstpar("gfdh","f",&intfac->dh);

if(strncmp(srctype,"fault",5) == 0)
   intfac->mode = FAULT;
else
   {
   fprintf(stderr,"*** invalid source type for interfacing, exiting...\n");
   exit(-1);
   }
}

get_fault_par(intfac)
struct interface *intfac;
{
struct finitefault *ff;

ff = &intfac->ffault;

mstpar("xs","f",&ff->xf);
mstpar("ys","f",&ff->yf);
mstpar("zs","f",&ff->zf);
mstpar("moment","f",&ff->moment);
mstpar("strike","f",&ff->strike);
mstpar("dip","f",&ff->dip);
mstpar("length","f",&ff->length);
mstpar("width","f",&ff->width);
mstpar("nstrk","d",&ff->nstrk);
mstpar("ndown","d",&ff->ndown);
mstpar("npntsrc","d",&ff->npntsrc);
mstpar("shypo","f",&ff->shypo);
mstpar("dhypo","f",&ff->dhypo);
mstpar("rvel","f",&ff->rvel);
mstpar("slipmodel","s",ff->slipmodel);
mstpar("srcmodel","s",ff->srcmodel);
}

init_fault(nx,ny,nz,nt,ntp2,s,h,dt,intfac,ms,mdp,nprof,med,fs)
struct interface *intfac;
struct medstencil *ms;
struct medprof *mdp;
float *h, *dt, *s, *med;
int nx, ny, nz, nprof, fs, ntp2;
{
FILE *fpr, *fopfile();
struct finitefault *ff;
struct boundflags *bf;
struct pointsource *ps;
struct sourceparams *srcparams;
struct gfsourcepar gfsrcpar, *gfsrcp;
struct gfindexpar gfindxpar, *gfindxp;
struct gftable *gftab;
float x, y, z, x2, y2, z2, lmin, lngth;
float arg, srcfac, slipfac;
float u, v, w, bigl2, len, wid, nrmslip, dep;
float cosA, sinA, cosD, sinD, xhypo, yhypo, zhypo;
float u0, v0;
int ix, iy, iz, ibnd, ip, fdw, iyspt, jstart, jend, j;
int go, i, ixs, ixe, iys, iye, izs, ize, xinc, zinc, nsrcnt, nbyte;
int isrc, istrk, idwn, jj, offset;
int vxf, vyf, vzf;
float tsmin, vsgrid;
int jsrc, is1, id1;
float invdt, dlen, dwid;
char ychar[8], fname[256];

float *gftmp, dx, dy, dz;
float cfx, cfy, cfz;
float tstx, h2, inv4pi;
float memalloc;
int idz, idh, gfindx, npts, ir, ixmod, iymod, izmod, ixst, iyst, izst;
int *gfgotab, goflag, fdgf;
/*
int igftmp[8];
*/
unsigned short igftmp[8];

float normf = 1.0e-20;
int size_int = sizeof(int);
int size_fl = sizeof(float);

int *tpts;
int maxbyte = 0;

float zap = 0.0;
float quart = 0.25;
float half = 0.5;
float one = 1.0;

int size_gfsrcp = sizeof(struct gfsourcepar);
int size_gfindxp = sizeof(struct gfindexpar);
gfsrcp = &gfsrcpar;
gfindxp = &gfindxpar;

inv4pi = one/(4.0*PI);
invdt = one/(*dt);

/*
   calculate number of subfaults and allocate memory
*/

ff = &intfac->ffault;
ff->nsrc = ff->nstrk*ff->ndown;
ff->npnt2 = ff->npntsrc*ff->npntsrc;

intfac->bflags = (struct boundflags *) check_malloc ((ff->nsrc)*sizeof(struct boundflags));

ff->pntsrc = (struct pointsource *) check_malloc ((ff->nsrc)*sizeof(struct pointsource));

srcparams = (struct sourceparams *) check_malloc (ff->npnt2*sizeof(struct sourceparams));

ff->slip = (float *) check_malloc ((ff->nsrc)*size_fl);
ff->rake = (float *) check_malloc ((ff->nsrc)*size_fl);
ff->goslip = (int *) check_malloc ((ff->nsrc)*size_fl);

/*
   get slip model, ordering is as follows:


      1   2   3   4   5   6   --> strike
      7   8   9  10  11  12
     13  14  15  16  17  18

*/

fpr = fopfile(ff->slipmodel,"r");
for(isrc=0;isrc<ff->nsrc;isrc++)
   {
   fscanf(fpr,"%f %f %d",&ff->slip[isrc],&ff->rake[isrc],&ff->goslip[isrc]);

   if(ff->goslip[isrc] != 0)
      ff->goslip[isrc] = 1;

   if(ff->slip[isrc] == 0.0)
      ff->goslip[isrc] = 0;

fprintf(stderr,"%3d) w=%10.5f r=%8.1f f=%d\n",isrc,ff->slip[isrc],ff->rake[isrc],ff->goslip[isrc]);
   }
fclose(fpr);

/*
   get source velocity model
*/

fpr = fopfile(ff->srcmodel,"r");

fscanf(fpr,"%d",&ff->nlay);
ff->vp = (float *) check_malloc ((ff->nlay)*size_fl);
ff->inva2 = (float *) check_malloc ((ff->nlay)*size_fl);
ff->vs = (float *) check_malloc ((ff->nlay)*size_fl);
ff->invb2 = (float *) check_malloc ((ff->nlay)*size_fl);
ff->vs2 = (float *) check_malloc ((ff->nlay)*size_fl);
ff->den = (float *) check_malloc ((ff->nlay)*size_fl);
ff->th = (float *) check_malloc ((ff->nlay)*size_fl);

for(i=0;i<ff->nlay;i++)
   {
   fscanf(fpr,"%f %f %f %f",&ff->vp[i],&ff->vs[i],&ff->den[i],&ff->th[i]);
   ff->inva2[i] = 1.0/(ff->vp[i]*ff->vp[i]);
   ff->vs2[i] = ff->vs[i]*ff->vs[i];
   ff->invb2[i] = 1.0/ff->vs2[i];
   }
fclose(fpr);

/*
   set-up source time function parameters
*/

make_rndm_stf(s,ntp2,dt);     /* random component of stf */

intfac->dt = (*dt);
intfac->nt = nt;
intfac->stf = s;

/*
   calculate subfault dimensions
*/

len = ff->length/ff->nstrk;
wid = ff->width/ff->ndown;

bigl2 = 0.5*ff->length;

/*
   calculate hypocenter location in (x,y,z) space
*/

arg = RPERD*ff->dip;
cosD = cos(arg);
sinD = sin(arg);

arg = RPERD*ff->strike;
cosA = cos(arg);
sinA = sin(arg);

u = -bigl2 + ff->shypo;
v = ff->dhypo;

xhypo = ff->xf + v*cosD*cosA + u*sinA;
yhypo = ff->yf + v*cosD*sinA - u*cosA;
zhypo = ff->zf + v*sinD;

/*
   Set time shift for source (tshft): Time shift is determined from the
   hypocenter to the nearest point on the model grid, using the velocity
   of the deepest layer for the source velocity model.  The following code
   is not elegant, but it works.
*/

/* OLD TSMIN CALCULATION */

lmin = 1.0e+15;

for(i=0;i<6;i++)
   {
   switch(i)
      {
      case 0:        /* x0 boundary */
         ixs = 0;    ixe = 1;  iys = 0;    iye = ny; izs = 0;    ize = nz;
         break;
      case 1:        /* x0 boundary */
         ixs = nx-1; ixe = nx; iys = 0;    iye = ny; izs = 0;    ize = nz;
         break;
      case 2:        /* y0 boundary */
         ixs = 0;    ixe = nx; iys = 0;    iye = 1;  izs = 0;    ize = nz;
         break;
      case 3:        /* yn boundary */
         ixs = 0;    ixe = nx; iys = ny-1; iye = ny; izs = 0;    ize = nz;
         break;
      case 4:        /* z0 boundary */
         ixs = 0;    ixe = nx; iys = 0;    iye = ny; izs = 0;    ize = 1;
         break;
      case 5:        /* zn boundary */
         ixs = 0;    ixe = nx; iys = 0;    iye = ny; izs = nz-1; ize = nz;
         break;
      }
   for(iz=izs;iz<ize;iz++)
      {
      z = iz*(*h) - zhypo;
      z2 = z*z;
      for(iy=iys;iy<iye;iy++)
         {
         y = iy*(*h) - yhypo;
         y2 = y*y;
         for(ix=ixs;ix<ixe;ix++)
            {
            x = ix*(*h) - xhypo;
            x2 = x*x;
            
            lngth = sqrt(x2 + y2 + z2);
            if(lngth < lmin)
               lmin = lngth;
            }
         }
      }
   }

/*
tsmin = 0.85*lmin/(ff->vs[(ff->nlay)-1]);
tsmin = lmin/(ff->vs[(ff->nlay)-1]);
fprintf(stderr,"**** Interfacing time shift = %f sec.\n\n",tsmin);
*/

/* NEW TSMIN CALCULATION */
tsmin = 0.0;
ff->tsmin = 1.0e+15;

/*
   calculate subfault parameters
*/

dlen = len/ff->npntsrc;
dwid = wid/ff->npntsrc;

nrmslip = 0.0;
isrc = 0;
for(idwn=0;idwn<ff->ndown;idwn++)
   {
   for(istrk=0;istrk<ff->nstrk;istrk++)
      {
      bf = &intfac->bflags[isrc];
      ps = &ff->pntsrc[isrc];

      ps->xs = (float *) check_malloc (ff->npnt2*size_fl);
      ps->ys = (float *) check_malloc (ff->npnt2*size_fl);
      ps->zs = (float *) check_malloc (ff->npnt2*size_fl);
      ps->rtime = (float *) check_malloc (ff->npnt2*size_fl);
      ps->vsrc = (float *) check_malloc (ff->npnt2*size_fl);
      ps->den = (float *) check_malloc (ff->npnt2*size_fl);
      ps->srcamp = (float *) check_malloc (ff->npnt2*size_fl);

      ps->itshft = tsmin/(*dt);    /*   set time shift for each subfault */

      ps->expl = 0;

      ps->nsrc = ff->nsrc;
      ps->isrc = isrc;

/* FD parameters */

      ps->h = (*h);
      ps->dt = (*dt);

/* fault parameters */

      ps->rake = ff->rake[isrc];
      arg = RPERD*ps->rake;
      ps->cosL = cos(arg);
      ps->sinL = sin(arg);

      arg = RPERD*ff->dip;
      ps->cosD = cos(arg);
      ps->sinD = sin(arg);

      arg = 2.0*arg;
      ps->cos2D = cos(arg);
      ps->sin2D = sin(arg);

      arg = RPERD*ff->strike;
      ps->cosA = cos(arg);
      ps->sinA = sin(arg);

/*
   calculate subfault initiation and termination times - coordinates
   measured with respect to hypocenter location
      u0 = along strike coordinate on left of subfault
      v0 = down-dip coordinate at top of subfault
*/

      jsrc = 0;
      for(id1=0;id1<ff->npntsrc;id1++)
	 {
         v0 = idwn*wid + (id1+0.5)*dwid;
         for(is1=0;is1<ff->npntsrc;is1++)
	    {
            u0 = istrk*len + (is1+0.5)*dlen;

            get_ruptime(&ps->rtime[jsrc],ff,&u0,&v0);
	    jsrc++;
	    }
	 }

/*
   find subfault corner locations relative to model origin

      u = along strike coordinate
      v = down-dip coordinate projected onto horizontal axis
      w = down-dip coordinate projected onto vertical axis
*/

      jsrc = 0;
      for(id1=0;id1<ff->npntsrc;id1++)
	 {
         v = (idwn*wid + (id1+0.5)*dwid)*ps->cosD;
         w = (idwn*wid + (id1+0.5)*dwid)*ps->sinD;
         for(is1=0;is1<ff->npntsrc;is1++)
	    {
            u = -bigl2 + istrk*len + (is1+0.5)*dlen;

            ps->xs[jsrc] = ff->xf + u*ps->sinA + v*ps->cosA;
            ps->ys[jsrc] = ff->yf - u*ps->cosA + v*ps->sinA;
            ps->zs[jsrc] = ff->zf + w;
	    jsrc++;
	    }
	 }

/*
   find grid boundaries illuminated by source -> use center of subfault
*/

      jsrc = ff->npnt2/2;

                         /* x boundaries */
      if(ps->xs[jsrc] < 0.0)
         {
         bf->x0bnd = 1; bf->xnbnd = 0;
         }
      else if(ps->xs[jsrc] >= 0.0 && ps->xs[jsrc] <= nx*(*h))
         {
         bf->x0bnd = 0; bf->xnbnd = 0;
         }
      else
         {
         bf->x0bnd = 0; bf->xnbnd = 1;
         }

                         /* y boundaries */
      if(ps->ys[jsrc] < 0.0)
         {
         bf->y0bnd = 1; bf->ynbnd = 0;
         }
      else if(ps->ys[jsrc] >= 0.0 && ps->ys[jsrc] <= ny*(*h))
         {
         bf->y0bnd = 0; bf->ynbnd = 0;
         }
      else
         {
         bf->y0bnd = 0; bf->ynbnd = 1;
         }

                         /* z boundaries */
      if(ps->zs[jsrc] < 0.0)
         {
         bf->z0bnd = 1; bf->znbnd = 0;
         }
      else if(ps->zs[jsrc] >= 0.0 && ps->zs[jsrc] <= nz*(*h))
         {
         bf->z0bnd = 0; bf->znbnd = 0;
         }
      else
         {
         bf->z0bnd = 0; bf->znbnd = 1;
         }

      if(bf->x0bnd == 0 && bf->xnbnd == 0 && bf->y0bnd == 0 && bf->ynbnd == 0 && bf->z0bnd == 0 && bf->znbnd == 0)
         {
         fprintf(stderr,"*** interface source inside FD grid, exiting...\n");
         exit(-1);
         }

/* 
   Source parameters -Source time function is normalized so that the it has
		      an area of unity, it is assumed that the input time
		      function is the far-field displacement pulse.  This
		      will be differentiated to velocity after normalization.
                     -Units of parameters are as follows:

		     moment  -> dyne-cm
		     density -> g/(cm*cm*cm)
		     vsrc    -> km/sec
		     rad     -> km [this term applied in calc_infc_pnt()]

		     normf = 1.0e-20 for unit normalization
*/

      for(jsrc=0;jsrc<ff->npnt2;jsrc++)
	 {
         jj = 0;
         dep = ff->th[jj];
         while(ps->zs[jsrc] > dep)
            {
            jj++;
            dep = dep + ff->th[jj];
            }

         ps->vsrc[jsrc] = ff->vs[jj];
         ps->den[jsrc] = ff->den[jj];

         nrmslip = nrmslip + (ps->vsrc[jsrc])*(ps->vsrc[jsrc])*(ps->den[jsrc])*(ff->slip[isrc]);
	 }

      ps->st = s;

      isrc++;
      }
   }

/*
   normalize point source amplitudes using rigidity, relative slip, and
   total moment.
*/

for(isrc=0;isrc<ff->nsrc;isrc++)
   {
   ps = &ff->pntsrc[isrc];

   for(jsrc=0;jsrc<ff->npnt2;jsrc++)
      {
      srcfac = (normf*ff->moment)*inv4pi/((ps->vsrc[jsrc])*(ps->vsrc[jsrc])*(ps->vsrc[jsrc])*(ps->den[jsrc]));

      slipfac = (ps->vsrc[jsrc])*(ps->vsrc[jsrc])*(ps->den[jsrc])*ff->slip[isrc]/nrmslip;

      ps->srcamp[jsrc] = ff->goslip[isrc]*slipfac*srcfac;
      }
   }

/*
   Calculate regular grid of source Green's functions.  These will be
   stored in the array intfac->gfvbuf.
*/

idz = (0.5 + intfac->dz/(*h));
idh = (0.5 + intfac->dh/(*h));
h2 = 0.5*(*h);

if(idh < 3)
   idh = 3;
if(idz < 3)
   idz = 3;

/*   determine number of gf's to calculate */

ixst = 0;
if(nx%idh == 1)
   ixst = -1;

iyst = 0;
if(ny%idh == 1)
   iyst = -1;

izst = 0;
if(nz%idz == 1)
   izst = -1;

gfindx = 0;
for(iy=iyst;(iy-idh)<ny;iy=iy+idh)
   {
   for(iz=izst;(iz-idz)<nz;iz=iz+idz)
      {
      for(ix=ixst;(ix-idh)<nx;ix=ix+idh)
         {
	 goflag = 0;
	 if(iy > iyst+idh && iy < ny-idh)        /* inside y boundaries */
	    {
	    if(iz > izst+idz && iz < nz-idz)     /* inside z boundaries */
	       {
	       if(ix > ixst+idh && ix < nx-idh)  /* inside x boundaries */
	          {
		  while(ix < nx - 2*idh)
		     ix = ix + idh;
	          }
	       else
		  goflag = 1;
	       }
	    else
	       goflag = 1;
	    }
	 else
	    goflag = 1;

         if(goflag)
            gfindx++;
         }
      }
   }

intfac->ngfsrc = gfindx;

/*   allocate memory for gf's and associated parameters */

gftmp = (float *) check_malloc (3*nt*size_fl);
gfgotab = (int *) check_malloc (intfac->ngfsrc*size_int);
intfac->gftable = (struct gftable *) check_malloc (intfac->ngfsrc*sizeof(struct gftable));

intfac->gfvbuf = (float **) check_malloc (ff->nsrc*sizeof(float *));
tpts = (int *) check_malloc (ff->nsrc*size_int);

sprintf(fname,"%s_tmpfile",intfac->gfname);
fdgf = croptrfile(fname);

/* calculate gf's */

gfindx = 0;
for(iy=iyst;(iy-idh)<ny;iy=iy+idh)
   {
   iymod = iy;
   if(iymod < 0)
      iymod = 0;
   if(iymod >= ny)
      iymod = ny - 1;

   for(iz=izst;(iz-idz)<nz;iz=iz+idz)
      {
      izmod = iz;
      if(izmod < 0)
         izmod = 0;
      if(izmod >= nz)
	 izmod = nz - 1;

      for(ix=ixst;(ix-idh)<nx;ix=ix+idh)
         {
	 ixmod = ix;
         if(ixmod < 0)
            ixmod = 0;
	 if(ixmod >= nx)
	    ixmod = nx - 1;

	 goflag = 0;
	 if(iy > iyst+idh && iy < ny-idh)        /* true-> inside y border */
	    {
	    if(iz > izst+idz && iz < nz-idz)     /* true-> inside z border */
	       {
	       if(ix > ixst+idh && ix < nx-idh)  /* true-> inside x border */
	          {
		  while(ix < (nx-2*idh))  /* model interior-> skip to border */
		     ix = ix + idh;
	          }
	       else
		  goflag = 1;
	       }
	    else
	       goflag = 1;
	    }
	 else
	    goflag = 1;

         if(goflag) /* true-> grid location is somewhere along model border */
	    {
            gfgotab[gfindx] = 0; /* reset later if gf is used */

	    jj = 0;
            dep = ff->th[jj];
	    while(izmod*(*h) > dep)
	       {
               jj++;
               dep = dep + ff->th[jj];
               }
            vsgrid = ff->vs[jj];

	    gftab = intfac->gftable + gfindx;
            gftab->ix = ix;
	    gftab->iy = iy;
	    gftab->iz = iz;

            gftab->ts = (float *) check_malloc (ff->nsrc*size_fl);
            gftab->npts = (int *) check_malloc (ff->nsrc*size_int);
	    gftab->gfv = (float **) check_malloc (ff->nsrc*sizeof(float *));

            for(isrc=0;isrc<ff->nsrc;isrc++)
	       {
	       zero(gftmp,3*nt);

               ps = &ff->pntsrc[isrc];
	       gftab->ts[isrc] = 9999999.9;

	       partition_ff(ff,ps,srcparams,gftab,isrc,gftmp,nt,&vsgrid);

	       if(gftab->npts[isrc] > nt)
	          gftab->npts[isrc] = nt;

/*
   gftmp + 0*nt:	gfvx
   gftmp + 1*nt:	gfvy
   gftmp + 2*nt:	gfvz
*/

               rite(fdgf,gftmp,gftab->npts[isrc]*size_fl);
               rite(fdgf,gftmp+nt,gftab->npts[isrc]*size_fl);
               rite(fdgf,gftmp+2*nt,gftab->npts[isrc]*size_fl);
	       }

	    gfindx++;
	    }
         }
      }
   }

/* NEW TSMIN CALCULATION */

ff->tsmin = ff->tsmin - 20*(*dt);
fprintf(stderr,"**** Interfacing time shift = %f sec.\n\n",ff->tsmin);

/*
   loop over model boundaries and calculate:
		itgs   - gf start time relative to ps->itshft
		npts   - # of time points that gf is active
		igf[j] - array indices for 8 surrounding gf's
		cf[j]  - averaging factors for 8 surrounding gf's
		vxf    - flag for vx response at grid point
		vyf    - flag for vy response at grid point
		vzf    - flag for vz response at grid point
   for each grid point
*/

for(iy=0;iy<ny;iy++)
   {
   if(iy < 10)
      sprintf(ychar,"00%1d",iy);
   else if(iy < 100)
      sprintf(ychar,"0%2d",iy);
   else if(iy < 1000)
      sprintf(ychar,"%3d",iy);

   for(ibnd=0;ibnd<2;ibnd++)
      {
      switch(ibnd)
         {
/*
   case 0: interior boundary edge
   case 1: exterior boundary edge
*/
         case 0:
            sprintf(fname,"%sint.%3s",intfac->gfname,ychar);
            ip = 1;
            break;
         case 1:
            sprintf(fname,"%sext.%3s",intfac->gfname,ychar);
            ip = 0;
            break;
         }
 
      fdw = croptrfile(fname);

      if(iy == 0 || iy == ny-1)
         {
	 iyspt = iy;
	 if(ibnd == 0)
	    {
	    if(iy == 0)
	       iyspt = 1;
	    if(iy == ny-1)
	       iyspt = ny-2;
	    }

         jstart = 4;
         jend = 5;
         }
      else
         {
	 iyspt = iy;
         jstart = 0;
         jend = 4;

	   /* store all gf information in one block when iy != 0,ny-1 */
         nbyte = 0;
         rite(fdw,&nbyte,size_int); /* replace later */
         }

      for(j=jstart;j<jend;j++)
         {
         switch(j)
            {
/*
   case 0: iz=0         boundary
   case 1: ix=nx-1      boundary
   case 2: iz=nz-1      boundary
   case 3: ix=0         boundary
   case 4: iy=0,ny-1    boundary
*/
            case 0:
               ixs = 1 + ip;      ixe = nx - ip; xinc = 1;
               izs = ip;          ize = 1 + ip;  zinc = 1;
               break;
            case 1:
               ixs = nx - 1 - ip; ixe = nx - ip; xinc = 1;
               izs = 1 + ip;      ize = nz - ip; zinc = 1;
               break;
            case 2:
               ixs = nx - 2 - ip; ixe = ip - 1;  xinc = -1;
               izs = nz - 1 - ip; ize = nz - ip; zinc = 1;
               break;
            case 3:
               ixs = ip;          ixe = 1 + ip;  xinc = 1;
               izs = nz - 2 - ip; ize = ip - 1;  zinc = -1;
               break;
            case 4:
               ixs = ip;          ixe = nx - ip; xinc = 1;
               izs = ip;          ize = nz - ip; zinc = 1;
               break;
            }

         for(iz=izs;iz!=ize;iz=iz+zinc)
            {
	    if(j == 4) /* store gf information in x-strips when iy = 0,ny-1 */
	       {
               nbyte = 0;
               rite(fdw,&nbyte,size_int); /* replace later */
	       }

            for(ix=ixs;ix!=ixe;ix=ix+xinc)
               {
               nsrcnt = 0;
               rite(fdw,&nsrcnt,size_int); /* replace later */
               nbyte = nbyte + size_int;

               for(isrc=0;isrc<ff->nsrc;isrc++)
                  {
                  vxf = vyf = vzf = 1;
                  ps = &ff->pntsrc[isrc];

                     /* determine if this source illuminates this boundary */

                  bf = &intfac->bflags[isrc];
                  go = 0;

                  if(j == 0)
                     {
                     if(bf->z0bnd)
		        {
                        go = 1;
	                if(ix+xinc == ixe)
	                   {
		           vxf = 0;
			   vyf = vzf = 1;
			   }
                        }
		     else if(bf->x0bnd && ix == ixs && ibnd == 1)
			go = 1;
                     else if(bf->xnbnd && ix+3*xinc == ixe && ibnd == 1)
			{
			go = 1;
			vxf = 1;
			vyf = vzf = 0;
			}
                     else if(bf->xnbnd && ix+2*xinc == ixe)
			{
			go = 1;
			if(ibnd == 0)
			   {
			   vxf = 1;
			   vyf = vzf = 0;
			   }
			}
                     else if(bf->xnbnd && ix+xinc == ixe)
			{
			go = 1;
			vxf = 0;
			vyf = vzf = 1;
			}
                     }
                  else if(j == 1)
                     {
                     if(bf->xnbnd)
			{
                        go = 1;
			if(iz+zinc == ize)
			   {
			   vzf = 0;
			   vxf = vyf = 1;
			   }
			}
                     else if(bf->z0bnd && iz == izs && ibnd == 1)
                        go = 1;
                     else if(bf->znbnd && iz+3*zinc == ize && ibnd == 1)
			{
			go = 1;
			vzf = 1;
			vxf = vyf = 0;
			}
                     else if(bf->znbnd && iz+2*zinc == ize)
			{
			go = 1;
			if(ibnd == 0)
			   {
			   vzf = 1;
			   vxf = vyf = 0;
			   }
			}
                     else if(bf->znbnd && iz+zinc == ize)
			{
                        go = 1;
			vzf = 0;
			vxf = vyf = 1;
			}
                     }
                  else if(j == 2)
                     {
                     if(bf->znbnd)
			{
                        go = 1;
			if(ix == ixs)
			   {
			   vxf = 0;
			   vyf = vzf = 1;
			   }
			}
                     else if(bf->xnbnd && ix == ixs)
			{
                        go = 1;
			vxf = 0;
			vyf = vzf = 1;
			}
                     else if(bf->xnbnd && ix == ixs+xinc && ibnd == 1)
			{
                        go = 1;
			vxf = 1;
			vyf = vzf = 0;
			}
                     else if(bf->x0bnd && ix+2*xinc == ixe && ibnd == 1)
                        go = 1;
                     else if(bf->x0bnd && ix+xinc == ixe)
                        go = 1;
                     }
                  else if(j == 3)
                     {
                     if(bf->x0bnd)
			{
                        go = 1;
			if(iz == izs)
			   {
			   vzf = 0;
			   vxf = vyf = 1;
			   }
			}
                     else if(bf->znbnd && iz == izs)
			{
                        go = 1;
			vzf = 0;
			vxf = vyf = 1;
			}
                     else if(bf->znbnd && iz == izs+zinc && ibnd == 1)
			{
                        go = 1;
			vzf = 1;
			vxf = vyf = 0;
			}
                     else if(bf->z0bnd && iz+2*zinc == ize && ibnd == 1)
                        go = 1;
                     else if(bf->z0bnd && iz+zinc == ize)
                        go = 1;
                     }
	          else if(j == 4)  /* iy=0 or iy=ny-1 */
		     {
	             if(bf->y0bnd && iy == 0)
	                go = 1;
	             else if(bf->ynbnd && iy == ny-1)
	                go = 1;
	             else
			{
			if(bf->z0bnd)
			   {
	                   if(iz == izs)
	                      go = 1;
	                   else if(iz == izs+zinc && ibnd == 1)
	                      go = 1;
			   }
	                if(bf->xnbnd)
			   {
	                   if(ix+xinc == ixe || ix+2*xinc == ixe)
	                      go = 1;
	                   else if(ix+3*xinc == ixe && ibnd == 1)
	                      go = 1;
			   }
	                if(bf->znbnd)
			   {
	                   if(iz+zinc == ize || iz+2*zinc == ize)
	                      go = 1;
	                   else if(iz+3*zinc == ize && ibnd == 1)
	                      go = 1;
			   }
	                if(bf->x0bnd)
			   {
	                   if(ix == ixs)
	                      go = 1;
	                   else if(ix == ixs+xinc && ibnd == 1)
	                      go = 1;
			   }
			}
		     }

		  /* if source has zero slip, don't do */

                  if(ff->goslip[isrc] == 0)
                     go = 0;

                  if(go)    /* write GF parameters */
		     {
		     nsrcnt++;
		     gfsrcp->isrc = (unsigned short)(isrc);

/*
   Set the bits in gfsrcp->vflag using the following table:

             vxf=0 vyf=0 vzf=0 -> vflag(bits)=0000000 [vflag=0]
             vxf=1 vyf=0 vzf=0 -> vflag(bits)=0000001 [vflag=1]
             vxf=0 vyf=1 vzf=0 -> vflag(bits)=0000010 [vflag=2]
             vxf=1 vyf=1 vzf=0 -> vflag(bits)=0000011 [vflag=3]
             vxf=0 vyf=0 vzf=1 -> vflag(bits)=0000100 [vflag=4]
             vxf=1 vyf=0 vzf=1 -> vflag(bits)=0000101 [vflag=5]
             vxf=0 vyf=1 vzf=1 -> vflag(bits)=0000110 [vflag=6]
             vxf=1 vyf=1 vzf=1 -> vflag(bits)=0000111 [vflag=7]

   Bits will be decoded later using ONE_MASK and bit shifting.
*/

                     gfsrcp->vflag = 0;

                     if(vxf)
                        gfsrcp->vflag = gfsrcp->vflag + 1;
                     if(vyf)
                        gfsrcp->vflag = gfsrcp->vflag + 2;
                     if(vzf)
                        gfsrcp->vflag = gfsrcp->vflag + 4;

/*
   Determine 8 surrounding gf's.  Do this only when 'nsrcnt=1', that is only
   for the first subfault source that is active at this grid location.  All
   subsequent subfault sources will have the same set of 8 surrounding gf
   locations.

   The 8 surrounding gf locations are found by "bracketing" the grid location
   sequentially along the x-, y-, and z-axes.  For each of the 8 gf locations,
   the distance between the grid location and the gf location is then
   determined and stored in the array 'rad'.
*/

		     if(nsrcnt == 1)
			{
			ir = 0;
                        for(gfindx=0;gfindx<intfac->ngfsrc;gfindx++)
			   {
	                   gftab = intfac->gftable + gfindx;

			   dx = (ix - gftab->ix);
			   if((dx >= 0.0 && dx < idh) || (dx < 0.0 && -dx <= idh))
			      {
			      dy = (iyspt - gftab->iy);
			      if((dy >= 0.0 && dy < idh) || (dy < 0.0 && -dy <= idh))
			         {
			         dz = (iz - gftab->iz);
			         if((dz >= 0.0 && dz < idz) || (dz < 0.0 && -dz <= idz))
				    {
			            gfindxp->igf[ir] = (unsigned short)(gfindx);
				    gfgotab[gfindx] = 1;

				    ir++;
				    if(ir > 8)
				       {
				       fprintf(stderr,"ir > 8, exiting ...\n");
				       exit(-88);
				       }
				    }
			         }
			      }
			   }

                        if(ir < 8)
                           {
                           fprintf(stderr,"ir=%d (< 8), exiting ...\n",ir);
                           exit(-88);
                           }

/*
   Now, sort the gf indicies acording to grid location.  The ordering is
   as follows (ix0 < ix1), (iy0 < iy1), (iz0 < iz1):

	gfindxp->igf[0]		ix0,iy0,iz0		dx+,dy+,dz+
	gfindxp->igf[1]		ix1,iy0,iz0		dx-,dy+,dz+
	gfindxp->igf[2]		ix0,iy1,iz0		dx+,dy-,dz+
	gfindxp->igf[3]		ix1,iy1,iz0		dx-,dy-,dz+
	gfindxp->igf[4]		ix0,iy0,iz1		dx+,dy+,dz-
	gfindxp->igf[5]		ix1,iy0,iz1		dx-,dy+,dz-
	gfindxp->igf[6]		ix0,iy1,iz1		dx+,dy-,dz-
	gfindxp->igf[7]		ix1,iy1,iz1		dx-,dy-,dz-

   This will allow for easy calculation of dx/idh, dy/idh, and dz/idz factors
   later on.
*/

                        for(ir=0;ir<8;ir++)
			   {
	                   gftab = intfac->gftable + gfindxp->igf[ir];

			   dx = (ix - gftab->ix);
			   dy = (iyspt - gftab->iy);
			   dz = (iz - gftab->iz);

			   if(dz >= 0.0)
			      {
			      if(dy >= 0.0)
				 {
			         if(dx >= 0.0)
				    igftmp[0] = gfindxp->igf[ir];
                                 else
				    igftmp[1] = gfindxp->igf[ir];
				 }
                              else
				 {
			         if(dx >= 0.0)
				    igftmp[2] = gfindxp->igf[ir];
                                 else
				    igftmp[3] = gfindxp->igf[ir];
				 }
			      }
                           else
			      {
			      if(dy >= 0.0)
				 {
			         if(dx >= 0.0)
				    igftmp[4] = gfindxp->igf[ir];
                                 else
				    igftmp[5] = gfindxp->igf[ir];
				 }
                              else
				 {
			         if(dx >= 0.0)
				    igftmp[6] = gfindxp->igf[ir];
                                 else
				    igftmp[7] = gfindxp->igf[ir];
				 }
			      }
			   }

                        for(ir=0;ir<8;ir++)
                           gfindxp->igf[ir] = igftmp[ir];

/*
   The interpolated value of the gf at x=a,y=b,z=c is given by

	g(a,b,c) = g0 + a*dg/dx + b*dg/dy + c*dgdz.

   Calculate the weighting factors 'cfx', 'cfy', 'cfz', ....
   The weights are based on the distance between the central grid location and
   each of the gf locations.

			cfx = (dx+0.25)/idh
			cfy = (dy+0.25)/idh
			cfz = (dz+0.25)/idz

*/

                        gftab = intfac->gftable + gfindxp->igf[0];

			dx = (ix - gftab->ix);
			dy = (iyspt - gftab->iy);
			dz = (iz - gftab->iz);

                        cfx = quart*((float)(dx+quart)/(float)(idh));
                        cfy = quart*((float)(dy+quart)/(float)(idh));
                        cfz = quart*((float)(dz+quart)/(float)(idz));

/*
   Finally, combine the individual weighting factors so that they can
   be used directly on the value of the GF at each location.  The result
   is just a linear sum of these factors applied to the appropriate GF.
*/

			gfindxp->cfac[0] = -cfx - cfy  - cfz;
			gfindxp->cfac[1] =  cfx - cfy  - cfz;
			gfindxp->cfac[2] = -cfx + cfy  - cfz;
			gfindxp->cfac[3] =  cfx + cfy  - cfz;
			gfindxp->cfac[4] = -cfx - cfy  + cfz;
			gfindxp->cfac[5] =  cfx - cfy  + cfz;
			gfindxp->cfac[6] = -cfx + cfy  + cfz;
			gfindxp->cfac[7] =  cfx + cfy  + cfz;

	                rite(fdw,gfindxp,size_gfindxp);
	                nbyte = nbyte + size_gfindxp;
			}

/*
   Calculate itgfst & npts for this 'isrc' using the 'cf' weights.  These
   will be different for each subfault source that is active at this grid
   location.
*/

/* Start with the value at the first interpolation point */

                     gftab = intfac->gftable + gfindxp->igf[0];
		     tstx = gftab->ts[isrc];

                     gfsrcp->npts = 0;
                     for(ir=0;ir<8;ir++)
                        {
                        gftab = intfac->gftable + gfindxp->igf[ir];

			tstx = tstx + gfindxp->cfac[ir]*gftab->ts[isrc];

                        if(gftab->npts[isrc] > (int)(gfsrcp->npts))
                           gfsrcp->npts = (unsigned short)(gftab->npts[isrc]);
                        }

/* NEW TSMIN CALCULATION */
                     gfsrcp->itgfst = (unsigned short)((tstx - ff->tsmin)*invdt);

	             rite(fdw,gfsrcp,size_gfsrcp);
	             nbyte = nbyte + size_gfsrcp;
		     }
	          }

               if(nsrcnt)    /* replace with correct value */
                  {
		  offset = nsrcnt*size_gfsrcp + size_gfindxp;
                  lseek(fdw,-(offset+size_int),SEEK_CUR);
                  rite(fdw,&nsrcnt,size_int);
                  lseek(fdw,offset,SEEK_CUR);
                  }
               }

            if(j == 4)    /* replace with correct value when iy = 0,ny-1*/
               {
               lseek(fdw,-(nbyte+size_int),SEEK_CUR);
               rite(fdw,&nbyte,size_int);
               lseek(fdw,nbyte,SEEK_CUR);

	       if(nbyte > maxbyte)
		  maxbyte = nbyte;
               }
	    }
         }

      if(jend == 4)    /* replace with correct value when iy != 0,ny-1*/
         {
         lseek(fdw,0,SEEK_SET);
         rite(fdw,&nbyte,size_int);

	 if(nbyte > maxbyte)
	    maxbyte = nbyte;
         }
      close(fdw);
      }
   }

/*
   Read in GF's from temporary file that will be used in simulation,
   ie., gfgotab[gfindx] = 1.
*/

/* First, find out how much internal memory is needed for GF's. */

for(isrc=0;isrc<ff->nsrc;isrc++)
   tpts[isrc] = 0;

intfac->ngfmem = 0;
for(gfindx=0;gfindx<intfac->ngfsrc;gfindx++)
   {
   if(gfgotab[gfindx])
      {
      gftab = intfac->gftable + gfindx;

      for(isrc=0;isrc<ff->nsrc;isrc++)
	 tpts[isrc] = tpts[isrc] + 3*gftab->npts[isrc];

      intfac->ngfmem = intfac->ngfmem + 1;
      }
   }

/* Allocate memory for GF's. */

memalloc = 0.0;
for(isrc=0;isrc<ff->nsrc;isrc++)
   {
   intfac->gfvbuf[isrc] = (float *) check_malloc (tpts[isrc]*size_fl);
   memalloc = memalloc + tpts[isrc]*size_fl;
   tpts[isrc] = 0;
   }

/* Now, set memory pointers and read in GF's */

lseek(fdgf,0,SEEK_SET);
nbyte = 0;
for(gfindx=0;gfindx<intfac->ngfsrc;gfindx++)
   {
   gftab = intfac->gftable + gfindx;

   if(gfgotab[gfindx])
      {
      for(isrc=0;isrc<ff->nsrc;isrc++)
	 {
	 npts = gftab->npts[isrc];
	 gftab->gfv[isrc] = intfac->gfvbuf[isrc] + tpts[isrc];

	 reed(fdgf,gftab->gfv[isrc],npts*size_fl);
	 reed(fdgf,gftab->gfv[isrc]+npts,npts*size_fl);
	 reed(fdgf,gftab->gfv[isrc]+2*npts,npts*size_fl);

	 tpts[isrc] = tpts[isrc] + 3*npts;
	 nbyte = nbyte + 3*npts*size_fl;

/*
fprintf(stderr,"XXXX gf[%d].ts[%d]= %f\n",gfindx,isrc,gftab->ts[isrc]);
*/
	 }
      }
   else
      {
      offset = 0;
      for(isrc=0;isrc<ff->nsrc;isrc++)
         offset = offset + 3*gftab->npts[isrc];

      lseek(fdgf,offset*size_fl,SEEK_CUR);
      nbyte = nbyte + offset*size_fl;
      }
   }

/* Remove temporary GF file */

close(fdgf);
sprintf(fname,"%s_tmpfile",intfac->gfname);
unlink(fname);

/*
   Free memory for temporary buffers that will no longer be used.
*/

for(isrc=0;isrc<ff->nsrc;isrc++)
   {
   ps = &ff->pntsrc[isrc];

   free(ps->xs);
   free(ps->ys);
   free(ps->zs);
   free(ps->rtime);
   free(ps->vsrc);
   free(ps->den);
   free(ps->srcamp);
   }
free(ff->pntsrc);

free(gftmp);
free(intfac->bflags);
free(srcparams);

/*
   Allocate memory for temporary buffer when reading in interfacing wave
   field parameters.  Use value of 'maxbyte' determined above.
*/

intfac->buf_int = (char *) check_malloc (maxbyte);
intfac->buf_ext = (char *) check_malloc (maxbyte);

fprintf(stderr,"***** maxbyte= %d\n",maxbyte);
fprintf(stderr,"  ngfsrc(tot)= %d\n",intfac->ngfsrc);
fprintf(stderr,"  ngfsrc(mem)= %d\n",intfac->ngfmem);
fprintf(stderr,"  GF memalloc= %8.1f Mb\n",memalloc/1.0e+06);
fprintf(stderr,"\n");
}

get_ruptime(rtime,ff,u0,v0)
struct finitefault *ff;
float *rtime, *u0, *v0;
{
float u, v;

u = *u0 - ff->shypo;
v = *v0 - ff->dhypo;

*rtime = sqrt(u*u + v*v)/(ff->rvel);
}

norm_stf(s,nt,dt)
float *s, *dt;
int nt;
{
float norm, sum, samp2, s2;
int go, ite, it;

/* 
   normalize source time function to have area of unity and differentiate
*/

sum = 0.0;
for(it=0;it<nt;it++)
   sum = sum + (*dt)*s[it];

norm = 1.0/sum;
for(it=0;it<nt;it++)
   s[it] = s[it]*norm;

/* this form actually shifts STF by 0.5*dt
norm = 1.0/(*dt);
s[0] = 0.0;
for(it=1;it<nt;it++)
   s[it] = (s[it+1] - s[it])*norm;
*/

/* find max. amp. of source time function and set cutoff at 1% of max. amp. */
samp2 = 0.0;
for(it=0;it<nt;it++)
   {
   s2 = s[it]*s[it];
   if(s2 > samp2)
      samp2 = s2;
   }

ite = nt;
go = 0;
samp2 = 0.0001*samp2;
for(it=0;it<nt;it++)
   {
   s2 = s[it]*s[it];

   if(s2 > samp2)
      go = 1;
   if(go && s2 < samp2)
      {
      ite = it;
      go = 0;
      }
   }
return(ite);
}

rm_infc_files(intfac)
struct interface *intfac;
{
int iy;
char ychar[8], fname[512];

for(iy=0;iy<intfac->ny;iy++)
   {
   if(iy < 10)
      sprintf(ychar,"00%1d",iy);
   else if(iy < 100)
      sprintf(ychar,"0%2d",iy);
   else if(iy < 1000)
      sprintf(ychar,"%3d",iy);

   sprintf(fname,"%sint.%3s",intfac->gfname,ychar);
   unlink(fname);
   sprintf(fname,"%sext.%3s",intfac->gfname,ychar);
   unlink(fname);
   }
}

check_infc_goflag(intfac) 
struct interface *intfac;
{
int *goptr, iy;

/* switch go flag pointers */

goptr = intfac->goplane;
intfac->goplane = intfac->newgo;
intfac->newgo = goptr;

intfac->go = 0;
for(iy=0;iy<intfac->ny;iy++)
   {
   if(intfac->goplane[iy])
      {
      intfac->go = 1;
      break;
      }
   }

for(iy=0;iy<intfac->ny;iy++)
   intfac->newgo[iy] = 0;
}

get_rayp(ff,h,srcd,recd,rng,p,rad,tt,eps)
struct finitefault *ff;
float *h, *rng, *recd, *p, *rad, *tt, *srcd, *eps;
{
float sdep, rdep, delr, sth, rth;
float pp, pm, tol;
float r0, r1, r2;
int k, slay, rlay, linc;

float zap = 0.0;
float half = 0.5;
float one = 1.0;
float tenth = 0.1;

tol = tenth*(*h);

k = 0;
sdep = ff->th[0];
while((*srcd) > sdep)
   {
   k++;
   sdep = sdep + ff->th[k];
   }
slay = k;

k = 0;
rdep = ff->th[0];
while((*recd) > rdep)
   {
   k++;
   rdep = rdep + ff->th[k];
   }
rlay = k;

if(slay != rlay)
   {
   if(sdep > rdep)
      {
      sth = ff->th[slay] - (sdep - *srcd);
      rth = rdep - *recd;
      linc = -1;
      }
   else
      {
      sth = sdep - *srcd;
      rth = ff->th[rlay] - (rdep - *recd);
      linc = 1;
      }

/*
   bisection method
   if(*p < *eps)
*/
   *p = 0.5*(*eps);
      bisect_p(ff,slay,&sth,rlay,&rth,p,eps,&tol,rng,linc);

/*
   Newton's method
   else
      newton_p(ff,slay,&sth,rlay,&rth,p,eps,&tol,rng,linc);
*/

         /* get path length and travel time for correct ray parameter */

   get_radtime(ff,slay,&sth,rlay,&rth,p,rad,tt,linc);
   }
else
   *p = -1.0;
}

/*
   07/15/94:  SV transmission coefficient seems to be very small.  Need to
   check accuracy of computation formulas in get_transmission().  For now,
   just use SH transmission coefficient for SV.
*/

get_transmission(ff,srcd,recd,p,tsh,tsv)
struct finitefault *ff;
float *recd, *p, *tsh, *tsv, *srcd;
{
float sdep, rdep;
float p2, gam0, gam1, arg;
int i, k, slay, rlay, linc;
float t0, t1, u0, u1, a0, b0, c0, d0, e0, f0, g0, h0, bigd;
float r0, i0, r1, i1, mag;
struct complex ec, gc, hc, bdc, tmp, inv;

float zap = 0.0;
float one = 1.0;
float two = 2.0;

k = 0;
sdep = ff->th[0];
while((*srcd) > sdep)
   {
   k++;
   sdep = sdep + ff->th[k];
   }
slay = k;

k = 0;
rdep = ff->th[0];
while((*recd) > rdep)
   {
   k++;
   rdep = rdep + ff->th[k];
   }
rlay = k;

/*
   Calculate SH transmission coefficients for layered stack.
   Aki and Richards, p. 144.
*/

*tsh = one;
if(slay != rlay)
   {
   if(sdep > rdep)
      linc = -1;
   else
      linc = 1;

   p2 = (*p)*(*p);

   for(i=slay;i!=rlay;i=i+linc)
      {
      k = i + linc;
      gam0 = ff->den[i]*ff->vs[i]*sqrt(one - ff->vs2[i]*p2);
      gam1 = ff->den[k]*ff->vs[k]*sqrt(one - ff->vs2[k]*p2);

      *tsh = (*tsh)*two*gam0/(gam0+gam1);
      }
   }

/*
   Calculate SV transmission coefficients for layered stack.
   Aki and Richards, p. 149-150.
*/

*tsv = one;
if(slay != rlay)
   {
   if(sdep > rdep)
      linc = -1;
   else
      linc = 1;

   p2 = (*p)*(*p);

   for(i=slay;i!=rlay;i=i+linc)
      {
      k = i + linc;

      u0 = two*ff->den[i]*ff->vs2[i];
      u1 = two*ff->den[k]*ff->vs2[k];
      t0 = ff->den[i] - u0*p2;
      t1 = ff->den[k] - u1*p2;

      a0 = t0 - t1;
      b0 = t0 + p2*u1;
      c0 = t1 + p2*u0;
      d0 = u0 - u1;

      arg = ff->inva2[i] - p2;
      if(arg >= zap)
	 {
         t0 = sqrt(arg);
	 r0 = one;
	 i0 = zap;
	 }
      else   /* cos(i) -> imaginary (inhomogeneous wave); set cos(i)=0.0 */
	 {
         t0 = sqrt(-arg);
	 r0 = zap;
	 i0 = one;
	 }

      arg = ff->inva2[k] - p2;
      if(arg >= zap)
	 {
         t1 = sqrt(arg);
	 r1 = one;
	 i1 = zap;
	 }
      else   /* cos(i) -> imaginary (inhomogeneous wave); set cos(i)=0.0 */
	 {
         t1 = sqrt(-arg);
	 r1 = zap;
	 i1 = one;
	 }

      u0 = sqrt(ff->invb2[i] - p2);
      u1 = sqrt(ff->invb2[k] - p2);

      ec.re = r1*b0*t1 + r0*c0*t0;
      ec.im = i1*b0*t1 + i0*c0*t0;

      f0 = b0*u1 + c0*u0;

      gc.re = a0 - r1*d0*t1*u0;
      gc.im = -i1*d0*t1*u0;

      hc.re = a0 - r0*d0*t0*u1;
      hc.im = -i0*d0*t0*u1;

      bdc.re = ec.re*f0 + (gc.re*hc.re - gc.im*hc.im)*p2;
      bdc.im = ec.im*f0 + (gc.re*hc.im + gc.im*hc.re)*p2;

      mag = bdc.re*bdc.re + bdc.im*bdc.im;
      inv.re = bdc.re/mag;
      inv.im = -bdc.im/mag;

      tmp.re = ec.re*inv.re - ec.im*inv.im;
      tmp.im = ec.re*inv.im + ec.im*inv.re;

      mag = sqrt(tmp.re*tmp.re + tmp.im*tmp.im);

      *tsv = (*tsv)*two*ff->den[i]*u0*ff->vs[i]*mag/ff->vs[k];

/*
      e0 = b0*t1 + c0*t0;
      f0 = b0*u1 + c0*u0;
      g0 = a0 - d0*t1*u0;
      h0 = a0 - d0*t0*u1;

      bigd = e0*f0 + g0*h0*p2;

      *tsv = (*tsv)*two*ff->den[i]*u0*e0*ff->vs[i]/(ff->vs[k]*bigd);
*/
      }
   }
}

get_range(ff,slay,sth,rlay,rth,p,r0,linc)
struct finitefault *ff;
float *sth, *rth, *p, *r0;
int linc, slay, rlay;
{
int i;
float denom, arg;
float invp2;
float one = 1.0;

invp2 = one/((*p)*(*p));

denom = sqrt(invp2*ff->invb2[slay] - one);
*r0 = (*sth)/denom;

for(i=slay+linc;i!=rlay;i=i+linc)
   {
   denom = sqrt(invp2*ff->invb2[i] - one);
   *r0 = *r0 + ff->th[i]/denom;
   }

denom = sqrt(invp2*ff->invb2[rlay] - one);
*r0 = *r0 + (*rth)/denom;
}

get_radtime(ff,slay,sth,rlay,rth,p,r0,tt,linc)
struct finitefault *ff;
float *sth, *rth, *p, *r0, *tt;
int linc, slay, rlay;
{
int i;
float r1, rad, arg;
float denom, invp2;
float one = 1.0;

if(*p > 0.0)
   {
   invp2 = one/((*p)*(*p));

   arg = invp2*ff->invb2[slay] - one;
   denom = sqrt(arg);
   r1 = (*sth)/denom;

   *r0 = sqrt(r1*r1 + (*sth)*(*sth));
   *tt = *r0/ff->vs[slay];

   for(i=slay+linc;i!=rlay;i=i+linc)
      {
      arg = invp2*ff->invb2[i] - one;
      denom = sqrt(arg);
      r1 = ff->th[i]/denom;

      rad = sqrt(r1*r1 + ff->th[i]*ff->th[i]);
      *r0 = *r0 + rad;
      *tt = *tt + rad/ff->vs[i];
      }

   arg = invp2*ff->invb2[rlay] - one;
   denom = sqrt(arg);
   r1 = (*rth)/denom;

   rad = sqrt(r1*r1 + (*rth)*(*rth));
   *r0 = *r0 + rad;
   *tt = *tt + rad/ff->vs[rlay];
   }
else
   {
   *r0 = *sth;
   *tt = *r0/ff->vs[slay];

   for(i=slay+linc;i!=rlay;i=i+linc)
      {
      *r0 = *r0 + ff->th[i];
      *tt = *tt + ff->th[i]/ff->vs[i];
      }

   *r0 = *r0 + *rth;
   *tt = *tt + (*rth)/ff->vs[rlay];
   }
}

newton_p(ff,slay,sth,rlay,rth,p,eps,tol,r0,linc)
struct finitefault *ff;
float *sth, *rth, *p, *r0, *eps, *tol;
int linc, slay, rlay;
{
int i;
float arg0, arg1, arg2, arg3;
float fp, fpp, tol2, p0;
float fofp;

float one = 1.0;
float two = 2.0;
int nc = 100;
int ic = 0;

p0 = one/ff->vs[slay];

tol2 = (*tol)*(*tol);
tol2 = two*(*eps);
fofp = two*tol2;

while(fofp > tol2)
   {
   ic++;
   fp = -(*r0);

   if(*p < *eps)
      *p = *eps;
   if(*p > p0)
      *p = p0*(one - *eps);

   arg1 = (*p)*ff->vs[slay];
   arg2 = one/(one - arg1*arg1);
   arg3 = sqrt(arg2);

   fp = fp + (*sth)*arg1*arg3;
   fpp = ff->vs[slay]*(*sth)*arg2*arg3;

   for(i=slay+linc;i!=rlay;i=i+linc)
      {
      arg1 = (*p)*ff->vs[i];
      arg2 = one/(one - arg1*arg1);
      arg3 = sqrt(arg2);

      fp = fp + (ff->th[i])*arg1*arg3;
      fpp = fpp + ff->vs[i]*(ff->th[i])*arg2*arg3;
      }

   arg1 = (*p)*ff->vs[rlay];
   arg2 = one/(one - arg1*arg1);
   arg3 = sqrt(arg2);

   fp = fp + (*rth)*arg1*arg3;
   fpp = fpp + ff->vs[rlay]*(*rth)*arg2*arg3;

   fofp = fp/(fpp);

   *p = *p - fofp;

   if(fofp < 0.0)
      fofp = -fofp;

   if(ic > nc)
      {
      fprintf(stderr,"(fp/fpp)= %e tol2=%e\n",fofp,tol2);
      fprintf(stderr,"p= %e p0= %e\n",*p,p0);
      fofp = 0.0;
      if(*p > p0)
         *p = p0*(one - *eps);
      }
   }
}

/*
*/
/*matherr(x)
register struct exception *x;
{
switch (x->type)
   {
   case DOMAIN:
      if (!strcmp(x->name, "sqrt"))
         {
         fprintf(stderr,"arg= %e\n",x->arg1);
         }
      break;
   }
return(0);
}*/

bisect_p(ff,slay,sth,rlay,rth,p,eps,tol,rng,linc)
struct finitefault *ff;
float *sth, *rth, *p, *rng, *eps, *tol;
int linc, slay, rlay;
{
float r0, delr, pp, pm, p0;
float tp0;

int i, ic;
int nc = 100;

float one = 1.0;
float half = 0.5;

p0 = one/ff->vs[slay];
for(i=slay+linc;i!=rlay;i=i+linc)
   {
   tp0 = one/ff->vs[i];
   if(tp0 < p0)
      p0 = tp0;
   }
tp0 = one/ff->vs[rlay];
if(tp0 < p0)
   p0 = tp0;

*p = *eps;

get_range(ff,slay,sth,rlay,rth,p,&r0,linc);

if(r0 < *rng)  /* if not, then p=0 (vertical ray) */
   {
		   /* bracket range with ray parameter extremes */

   *p = p0*(one - *eps);
   get_range(ff,slay,sth,rlay,rth,p,&r0,linc);

   pp = *p;
   pm = *eps;

   delr = r0 - *rng;

/*
   use bisection to converge to correct ray parameter
*/

   ic = 0;
   while(delr > *tol)
      {
      *p = half*(pp + pm);

      if(*p == pp || *p == pm) /* beyond single precision accuracy */
         break;

      get_range(ff,slay,sth,rlay,rth,p,&r0,linc);
      if(r0 >= *rng)
         {
         delr = r0 - *rng;
         pp = *p;
         }
      else
         {
         delr = *rng - r0;
         pm = *p;
         }

      ic++;
      if(ic > nc)
         break;
      }
   }
else
   *p = 0.0;
}

make_rndm_stf(ut,nt,dt)
float *ut, *dt;
int nt;
{
struct complex *uc;
float f0 = 1.0;
float f1 = 0.1;

uc = (struct complex *) ut;

make_stf(ut,nt);
forfft(uc,nt,1);
wfilter_stf(uc,nt,dt,&f0,&f1);
invfft(uc,nt,-1);
normal_stf(ut,nt,dt);
}

make_stf(u,nt)
float *u;
int nt;
{
int it;
double frand();

for(it=0;it<nt;it++)
   u[it] = frand();

for(it=0;it<nt;it++)            /* take derivative */
   u[it] = u[it+1] - u[it];

u[nt-1] = u[nt-2];
}

wfilter_stf(uc,nt,dt,f0,f1)
struct complex *uc;
float *f0, *f1, *dt;
int nt;
{
float df, fac, ff, freq;
int it;
float one = 1.0;

uc[0].re = 0.0;
uc[0].im = 0.0;

df = one/(nt*(*dt));
for(it=1;it<nt/2;it++)
   {
   freq = it*df;

   ff = freq/(*f0);
   ff = ff*ff*ff*ff;
   fac = one/(one + ff);

   ff = (*f1)/freq;
   ff = ff*ff*ff*ff;
   fac = fac/(one + ff);

   uc[it].re = fac*uc[it].re;
   uc[it].im = fac*uc[it].im;
   }
}

normal_stf(u,nt,dt)
float *u, *dt;
int nt;
{
FILE *fpw, *fopfile();
int nt6, j;
int it, i;
float max = 0.0;
 
for(it=0;it<nt;it++)
   {
   if(u[it] > max)
      max = u[it];
   else if(-u[it] > max)
      max = -u[it];
   }
 
max = 1.0/max;
 
for(it=0;it<nt;it++)
   u[it] = (*dt)*u[it]*max;

fpw = fopfile("stf_temp","w");

fprintf(fpw,"stf       all No Title\n");
 
fprintf(fpw,"%d %f 0 0 0.0 0.0 0.0 0.0\n",nt,*dt);

nt6 = nt/6;
for(i=0;i<nt6;i++)
   {
   for(j=0;j<6;j++)
      fprintf(fpw,"%13.5e",u[6*i + j]);
 
   fprintf(fpw,"\n");
   }
if(6*nt6 != nt)
   {
   for(i=6*nt6;i<nt;i++)
      fprintf(fpw,"%13.5e",u[i]);
 
   fprintf(fpw,"\n");
   }
fclose(fpw);
}

get_gfamp(gftab,isrc,it,k,igf,cfac,amp)
struct gftable *gftab;
float *amp, *cfac;
int it, k;
unsigned short *isrc, *igf;
{
struct gftable *gt;
int j;
float *gf;
float zap = 0.0;

/* Start with the value at the first interpolation point */

gt = gftab + igf[0];
gf = gt->gfv[*isrc] + k*gt->npts[*isrc] + it;
*amp = gf[0];

for(j=0;j<8;j++)
   {
   gt = gftab + igf[j];
   if(it < gt->npts[*isrc])
      {
      gf = gt->gfv[*isrc] + k*gt->npts[*isrc] + it;
 
      *amp = *amp + cfac[j]*gf[0];
      }
   }
}

dotp(a,v1,v2)
float *a;
struct vector *v1, *v2;
{
*a = (v1->x)*(v2->x) + (v1->y)*(v2->y) + (v1->z)*(v2->z);
}
