/*
   tstep.c contains the follow functions:

      tsteppbnd()
      tstepp()
      tstepvbnd()
      tstepv()
      fsurf()
*/

#include "include.h"

/*
	tsteppbnd() updates the stress variables txx, tyy, tzz, txy, txz
        and tyz for the new time step when iy=0 or iy=ny-1.  The varibles
        txx, tyy, tzz and txz are centered in plane 1,3 when
        iy=0,ny-1; respectively, and the variables txy and tyz are centered
        in plane 0,2 when iy=0,ny-1; respectively.  The field txz is not needed
        when iy=ny-1.
*/

tsteppbnd(nx,nz,pvptr,medptr,fdc,fs,iflag,eflag)
struct fdcoefs *fdc;
float **pvptr, **medptr;
int nx, nz, fs, iflag, eflag;
{
float *vx0, *vx1;
float *vy0, *vy1;
float *vz0, *vz1;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *taup, *vp, *med;
float *tp1, *tp2, *tp3, *med1, *med2;
float *cp1, *cp2, *cp3, *pkp, *skp;
float *vp0, *vp1;
float *l2mp, *lp, *mxyp, *mxzp, *myzp;
float *ak, *pk, *sk;
int ix, iz, ord, ordx, ip, ipm;

/*
   ip = 0 when iflag = YZERO (iy=0) -> solve all.
   ip = 2 when iflag = YN (iy=ny-1) -> solve only txy and tyz.

   ipm = 1 when iflag = YZERO
   ipm = 0 when iflag = YN
           => this shifts media pointer to be 1 plane inside
	      boundary (consistent with absorbing boundary condition)
*/

if(iflag == YZERO)
   {
   ipm = 1;
   ip = 0;
   }
else if(iflag == YN)
   {
   ipm = 0;
   ip = 2;
   }

l2mp  = medptr[ip+ipm];
lp    = medptr[ip+ipm] +   nx*nz;
mxyp  = medptr[ip+ipm] + 2*nx*nz;
mxzp  = medptr[ip+ipm] + 3*nx*nz;
myzp  = medptr[ip+ipm] + 4*nx*nz;

ak    = medptr[ip+ipm] +  8*nx*nz;
pk    = medptr[ip+ipm] +  9*nx*nz;
sk    = medptr[ip+ipm] + 10*nx*nz;

vx0 = pvptr[ip];
vx1 = pvptr[ip+1];

vy0 = pvptr[ip]   + nx*nz;
vy1 = pvptr[ip+1] + nx*nz;

vz0 = pvptr[ip]   + 2*nx*nz;
vz1 = pvptr[ip+1] + 2*nx*nz;

txx = pvptr[ip+1] + 3*nx*nz;
tyy = pvptr[ip+1] + 4*nx*nz;
tzz = pvptr[ip+1] + 5*nx*nz;
txy = pvptr[ip]   + 6*nx*nz; /* centered on plane 0,2 */
txz = pvptr[ip+1] + 7*nx*nz;
tyz = pvptr[ip]   + 8*nx*nz; /* centered on plane 0,2 */

cxx = pvptr[ip+1] +  9*nx*nz;
cyy = pvptr[ip+1] + 10*nx*nz;
czz = pvptr[ip+1] + 11*nx*nz;
cxy = pvptr[ip]   + 12*nx*nz; /* centered on plane 0,2 */
cxz = pvptr[ip+1] + 13*nx*nz;
cyz = pvptr[ip]   + 14*nx*nz; /* centered on plane 0,2 */

if(eflag == 0)
   {
   mv_mult1(txx,cxx,ak,nx*nz);
   mv_mult1(tyy,cyy,ak,nx*nz);
   mv_mult1(tzz,czz,ak,nx*nz);
   mv_mult1(txy,cxy,ak,nx*nz);
   mv_mult1(txz,cxz,ak,nx*nz);
   mv_mult1(tyz,cyz,ak,nx*nz);
   }

for(iz=1;iz<nz;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

/* x derivative */

   if(iflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = txx  + iz*nx + 1;
      tp2  = tyy  + iz*nx + 1;
      tp3  = tzz  + iz*nx + 1;
      vp   = vx1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cxx  + iz*nx + 1;
      cp2  = cyy  + iz*nx + 1;
      cp3  = czz  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
	 {
	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
	 }
      }

   if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp   = vy0  + iz*nx + 1;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sk   + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }

   if(iflag == YZERO) /* txz not needed at iy = ny-1 */
      {
      taup = txz  + (iz-1)*nx;
      vp   = vz1  + (iz-1)*nx + 1;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }
  
/* y derivative */

   if(iflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = tyy  + iz*nx + 1;
      tp2  = tzz  + iz*nx + 1;
      tp3  = txx  + iz*nx + 1;
      vp0  = vy0  + iz*nx + 1;
      vp1  = vy1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cyy  + iz*nx + 1;
      cp2  = czz  + iz*nx + 1;
      cp3  = cxx  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
	 {
	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp1,vp0,vp,vp,fdc,2,YDERIV);
	 }
      }

   taup = tyz  + (iz-1)*nx + 1;
   vp0  = vz0  + (iz-1)*nx + 1;
   vp1  = vz1  + (iz-1)*nx + 1;
   med  = myzp + (iz-1)*nx + 1;

   cp1  = cyz  + (iz-1)*nx + 1;
   skp  = sk   + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */

   if(eflag == 0)
      diff_mv(nx-2,taup,cp1,med,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
   else
      diff(nx-2,taup,med,vp1,vp0,vp,vp,fdc,2,YDERIV);

   if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp0  = vx0  + iz*nx;
      vp1  = vx1  + iz*nx;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sk   + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
      else
         diff(nx-1,taup,med,vp1,vp0,vp,vp,fdc,2,YDERIV);
      }
  
/* z derivative */

   if(iflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = tzz  + iz*nx + 1;
      tp2  = txx  + iz*nx + 1;
      tp3  = tyy  + iz*nx + 1;
      vp   = vz1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = czz  + iz*nx + 1;
      cp2  = cxx  + iz*nx + 1;
      cp3  = cyy  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1)
         {
         if(iz == 1 || iz == nz-2)
            ordx = 2;
         else
            ordx = ord;

	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
         }
      }

   if(iz == 1 || iz == nz-1)
      ordx = 2;
   else
      ordx = ord;

   if(iflag == YZERO) /* txz not needed at iy = ny-1 */
      {
      taup = txz  + (iz-1)*nx;
      vp   = vx1  +     iz*nx;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }

   taup = tyz  + (iz-1)*nx + 1;
   vp   = vy0  +     iz*nx + 1;
   med  = myzp + (iz-1)*nx + 1;

   cp1  = cyz  + (iz-1)*nx + 1;
   skp  = sk   + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
   if(eflag == 0)
      diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
   else
      diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
   }

if(eflag == 0)
   {
   mv_mult2(txx,cxx,nx*nz);
   mv_mult2(tyy,cyy,nx*nz);
   mv_mult2(tzz,czz,nx*nz);
   mv_mult2(txy,cxy,nx*nz);
   mv_mult2(txz,cxz,nx*nz);
   mv_mult2(tyz,cyz,nx*nz);
   }
}

/*
	tstepp() updates the stress variables txx, tyy, tzz, txy, txz and tyz
	for the new time step.  The varibles txx, tyy, tzz and txz are centered
	in plane 2 and the variables txy and tyz are centered in plane 1.
*/

tstepp(nx,nz,pvptr,medptr,fdc,fs,iflag,eflag)
struct fdcoefs *fdc;
float **pvptr, **medptr;
int nx, nz, fs, iflag, eflag;
{
float *vx0, *vx1, *vx2, *vx3;
float *vy0, *vy1, *vy2, *vy3;
float *vz0, *vz1, *vz2, *vz3;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *taup, *vp, *med;
float *tp1, *tp2, *tp3, *med1, *med2;
float *cp1, *cp2, *cp3, *pkp, *skp;
float *vp0, *vp1, *vp2, *vp3;
float *l2mp, *lp, *mxyp, *mxzp, *myzp;
float *ak, *pk, *sk, *aky, *sky;
int ix, iz, ord, ordx, j;

l2mp  = medptr[2];
lp    = medptr[2] +   nx*nz;
mxyp  = medptr[1] + 2*nx*nz; /* centered on plane 1 for txy */
mxzp  = medptr[2] + 3*nx*nz; /* centered on plane 2 for txz */
myzp  = medptr[1] + 4*nx*nz; /* centered on plane 1 for tyz */

ak    = medptr[2] +  8*nx*nz;
aky   = medptr[1] +  8*nx*nz; /* centered on plane 1 for txy,tyz */
pk    = medptr[2] +  9*nx*nz;
sk    = medptr[2] + 10*nx*nz;
sky   = medptr[1] + 10*nx*nz; /* centered on plane 1 for txy,tyz */
 
vx0 = pvptr[0];
vx1 = pvptr[1];
vx2 = pvptr[2];
vx3 = pvptr[3];

vy0 = pvptr[0] +   nx*nz;
vy1 = pvptr[1] +   nx*nz;
vy2 = pvptr[2] +   nx*nz;
vy3 = pvptr[3] +   nx*nz;

vz0 = pvptr[0] + 2*nx*nz;
vz1 = pvptr[1] + 2*nx*nz;
vz2 = pvptr[2] + 2*nx*nz;
vz3 = pvptr[3] + 2*nx*nz;

txx = pvptr[2] + 3*nx*nz;
tyy = pvptr[2] + 4*nx*nz;
tzz = pvptr[2] + 5*nx*nz;
txy = pvptr[1] + 6*nx*nz; /* centered on plane 1 */
txz = pvptr[2] + 7*nx*nz;
tyz = pvptr[1] + 8*nx*nz; /* centered on plane 1 */

cxx = pvptr[2] +  9*nx*nz;
cyy = pvptr[2] + 10*nx*nz;
czz = pvptr[2] + 11*nx*nz;
cxy = pvptr[1] + 12*nx*nz; /* centered on plane 1 */
cxz = pvptr[2] + 13*nx*nz;
cyz = pvptr[1] + 14*nx*nz; /* centered on plane 1 */

if(eflag == 0)
   {
   if(iflag != SKIP_PLANE_Y2)
      {
      mv_mult1(txx,cxx,ak,nx*nz);
      mv_mult1(tyy,cyy,ak,nx*nz);
      mv_mult1(tzz,czz,ak,nx*nz);
      mv_mult1(txz,cxz,ak,nx*nz);
      }
   if(iflag != SKIP_PLANE_Y1)
      {
      mv_mult1(txy,cxy,aky,nx*nz);
      mv_mult1(tyz,cyz,aky,nx*nz);
      }
   }

for(iz=1;iz<nz;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

   ordx = ord;

/* x derivative */

   if(iflag != SKIP_PLANE_Y2)
      {
      tp1  = txx  + iz*nx + 1;
      tp2  = tyy  + iz*nx + 1;
      tp3  = tzz  + iz*nx + 1;
      vp   = vx2  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cxx  + iz*nx + 1;
      cp2  = cyy  + iz*nx + 1;
      cp3  = czz  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
         {
         if(iz == 1 && fs && l2mp[nx] < 0.0)
            ordx = 2;
         else
            ordx = ord;

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-1,vp+1,vp-2,fdc,ordx,XDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-1,vp+1,vp-2,fdc,ordx,XDERIV);
         }

      taup = txz  + (iz-1)*nx;
      vp   = vz2  + (iz-1)*nx + 1;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }

   if(iz < (nz-1) && iflag != SKIP_PLANE_Y1) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp   = vy1  + iz*nx + 1;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sky  + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }
  
/* y derivative */

   if(iflag != SKIP_PLANE_Y2)
      {
      tp1  = tyy  + iz*nx + 1;
      tp2  = tzz  + iz*nx + 1;
      tp3  = txx  + iz*nx + 1;
      vp0  = vy0  + iz*nx + 1;
      vp1  = vy1  + iz*nx + 1;
      vp2  = vy2  + iz*nx + 1;
      vp3  = vy3  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cyy  + iz*nx + 1;
      cp2  = czz  + iz*nx + 1;
      cp3  = cxx  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
         {
         if((iz == 1 && fs && l2mp[nx] < 0.0) || iflag == YN)
            ordx = 2;
         else
            ordx = ord;

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp2,vp1,vp3,vp0,fdc,ordx,YDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp2,vp1,vp3,vp0,fdc,ordx,YDERIV);
         }
      }

   if(iflag != SKIP_PLANE_Y1)
      {
      taup = tyz  + (iz-1)*nx + 1;
      vp0  = vz0  + (iz-1)*nx + 1;
      vp1  = vz1  + (iz-1)*nx + 1;
      vp2  = vz2  + (iz-1)*nx + 1;
      vp3  = vz3  + (iz-1)*nx + 1;
      med  = myzp + (iz-1)*nx + 1;

      cp1  = cyz  + (iz-1)*nx + 1;
      skp  = sky  + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
      if(eflag == 0)
         diff_mv(nx-2,taup,cp1,med,skp,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
      else
         diff(nx-2,taup,med,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);

      if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
         {
         taup = txy  + iz*nx;
         vp0  = vx0  + iz*nx;
         vp1  = vx1  + iz*nx;
         vp2  = vx2  + iz*nx;
         vp3  = vx3  + iz*nx;
         med  = mxyp + iz*nx;

         cp1  = cxy  + iz*nx;
         skp  = sky  + iz*nx;

         if(eflag == 0)
            diff_mv(nx-1,taup,cp1,med,skp,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
         else
            diff(nx-1,taup,med,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
         }
      }
  
/* do x-strips: z derivative */

   if(iflag != SKIP_PLANE_Y2)
      {
      tp1  = tzz  + iz*nx + 1;
      tp2  = txx  + iz*nx + 1;
      tp3  = tyy  + iz*nx + 1;
      vp   = vz2  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = czz  + iz*nx + 1;
      cp2  = cxx  + iz*nx + 1;
      cp3  = cyy  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1)
         {
         ordx = ord;
         if(iz == 1 || iz == nz-2)
            ordx = 2;

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
         }

      ordx = ord;
      if(iz == 1 || iz == nz-1)
         ordx = 2;

      taup = txz  + (iz-1)*nx;
      vp   = vx2  +     iz*nx;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }

   if(iflag != SKIP_PLANE_Y1)
      {
      ordx = ord;
      if(iz == 1 || iz == nz-1)
         ordx = 2;

      taup = tyz  + (iz-1)*nx + 1;
      vp   = vy1  +     iz*nx + 1;
      med  = myzp + (iz-1)*nx + 1;

      cp1  = cyz  + (iz-1)*nx + 1;
      skp  = sky  + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
      if(eflag == 0)
         diff_mv(nx-2,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-2,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }
   }

if(eflag == 0)
   {
   if(iflag != SKIP_PLANE_Y2)
      {
      mv_mult2(txx,cxx,nx*nz);
      mv_mult2(tyy,cyy,nx*nz);
      mv_mult2(tzz,czz,nx*nz);
      mv_mult2(txz,cxz,nx*nz);
      }
   if(iflag != SKIP_PLANE_Y1)
      {
      mv_mult2(txy,cxy,nx*nz);
      mv_mult2(tyz,cyz,nx*nz);
      }
   }

 /* apply 2nd order boundary condition at iy=0,ny-1 */
if(iflag == YZERO || iflag == YN)
   tsteppbnd(nx,nz,pvptr,medptr,fdc,fs,iflag,eflag);

if(fs)    /* do freesurface for tzz, txz, and tyz */
   fsurf(nx,nz,pvptr,medptr,iflag,TAU);
   /*
   fsurfSHEAR(nx,nz,pvptr,medptr,iflag,TAU);
   */
}

/*
	When iflag=YZERO (near iy=0 boundary), tstepvbnd() updates the velocity
        variables vx and vz for the new time step on plane 1 using 2nd order
        y derivatives.
      
	The values of vx, vy and vz on plane 0 (iflag=YZERO) or plane 3
        (iflag=YN) are updated using an absorbing boundary condition in the
        function abs_ybnd().
*/

tstepvbnd(nx,nz,pvptr,medptr,pbnd,fdc,fs,iflag,intfac)
struct interface *intfac;
struct fdcoefs *fdc;
float **pvptr, **medptr, *pbnd;
int nx, nz, fs, iflag;
{
float *vx, *vy, *vz, *txx, *tzz, *txz;
float *txy0, *txy1;
float *tyz0, *tyz1;
float *tyy2, *tyy3;
float *txy2;
float *tyz2;
float *vxi, *vyi, *vzi;
float *taup, *vp, *med;
float *tp0, *tp1, *tp2, *tp3;
float *bxp, *byp, *bzp;
int ix, iz, ord, ordx, off1, off2, off3, ip;
 
/*
         iflag=YZERO -> do interior solution for vx and vz;
                        vy already done in tstepv().
         iflag=YN    -> go to abs_ybnd(), vx,vy and vz already done in tstepv().
*/

if(iflag == YZERO) /* set pointers for vx and vz update */
   {
/*
   vy, vyi not solved, but needed for store_iring()
*/
   bxp = medptr[1] + 5*nx*nz;
   byp = medptr[1] + 6*nx*nz;
   bzp = medptr[1] + 7*nx*nz;

   vx  = pvptr[1];
   vy  = pvptr[1] +   nx*nz;
   vz  = pvptr[1] + 2*nx*nz;
   txx = pvptr[1] + 3*nx*nz;
   tzz = pvptr[1] + 5*nx*nz;
   txz = pvptr[1] + 7*nx*nz;

   txy0 = pvptr[0] + 6*nx*nz;
   txy1 = pvptr[1] + 6*nx*nz;

   tyz0 = pvptr[0] + 8*nx*nz;
   tyz1 = pvptr[1] + 8*nx*nz;
 
   /* retain previous time step on inner ring for abs() */

   vxi = pbnd + 3*nx*nz;
   vyi = pbnd + 3*nx*nz + 2*(nx+nz);
   vzi = pbnd + 3*nx*nz + 4*(nx+nz);
   store_iring(vx,vxi,vy,vyi,vz,vzi,nx,nz);

   for(iz=1;iz<nz-1;iz++)
      {
      if(iz < fdc->izord2)
         ord = fdc->order;
      else
         ord = fdc->ordx;

/* x derivative */

      taup = txx + iz*nx + 2;
      vp   = vx  + iz*nx + 1;
      med  = bxp + iz*nx + 1;

      diff(nx-3,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

      taup = txz + iz*nx + 1;
      vp   = vz  + iz*nx + 1;
      med  = bzp + iz*nx + 1;

      if(iz < nz-2)
         diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);
   
/* 2nd order y derivative */

      tp0 = txy0 + iz*nx + 1;
      tp1 = txy1 + iz*nx + 1;
      vp  = vx   + iz*nx + 1;
      med = bxp  + iz*nx + 1;

      diff(nx-3,vp,med,tp1,tp0,taup,taup,fdc,2,YDERIV);

      tp0 = tyz0 + iz*nx + 1;
      tp1 = tyz1 + iz*nx + 1;
      vp  = vz   + iz*nx + 1;
      med = bzp  + iz*nx + 1;
   
      if(iz < nz-2)
         diff(nx-2,vp,med,tp1,tp0,taup,taup,fdc,2,YDERIV);
   
/* z derivative */

      if(iz == 1 || iz == nz-2)
         ordx = 2;
      else
         ordx = ord;

      taup = txz + iz*nx + 1;
      vp   = vx  + iz*nx + 1;
      med  = bxp + iz*nx + 1;

      diff(nx-3,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);

      taup = tzz + (iz+1)*nx + 1;
      vp   = vz  +     iz*nx + 1;
      med  = bzp +     iz*nx + 1;

      if(iz < nz-2)
         {
      /* only have 4th-order information above free-surface for tzz */
         if(iz == 1 && fs)
            ordx = ord;
         if(iz == nz-3)
            ordx = 2;

         diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
	 }
      }

   /* apply absorbing boundary along x and z edges */
   abs_xzbnd1(nx,nz,pvptr,pbnd,&fdc->dtoh,medptr,fs,intfac);
   }

/* absorb vx and vz on plane at iy=0,ny-1 */
/* absorb vy on plane at iy=0,ny-2 */
abs_ybnd(nx,nz,pvptr,medptr,pbnd,fdc,fs,iflag,intfac);
}

/*
	tstepv() updates the velocity variables vx, vy, and vz for the new
	time step.  The varibles vx and vz are centered in plane 2 and the
	variable vy is centered in plane 1.  These variables are absorbed
	at the x, y and z boundaries.

	if iflag=YZERO -> 2nd order y derivative of tyy, affects vy near iy=0;
                          also calls function tstepvbnd() to update velocities
                          at iy=0 boundary.
	if iflag=YN    -> 2nd order y derivative of tyy, txy and txz, affects
                          vx, vy and vz near iy=ny-1;  also calls function
			  tstepvbnd() to update vx and vz at iy=ny-1 boundary.
*/

tstepv(nx,nz,pvptr,medptr,pbnd,fdc,fs,iflag,intfac)
struct interface *intfac;
struct fdcoefs *fdc;
float **pvptr, **medptr, *pbnd;
int nx, nz, fs, iflag;
{
float *vx, *vy, *vz, *txx, *tzz, *txz;
float *tyy0, *tyy1, *tyy2, *tyy3;
float *txy0, *txy1, *txy2, *txy3;
float *tyz0, *tyz1, *tyz2, *tyz3;
float *vxi, *vyi, *vzi;
float *taup, *vp, *med, *med1, *med2;
float *tp0, *tp1, *tp2, *tp3;
float *bxp, *byp, *bzp;
int ix, iz, ord, ordx, off1, off2, off3;
 
/* do interior solution */

bxp  = medptr[2] + 5*nx*nz; /* center on plane 2 for vx */
byp  = medptr[1] + 6*nx*nz; /* center on plane 1 for vy */
bzp  = medptr[2] + 7*nx*nz; /* center on plane 2 for vz */

vx  = pvptr[2];
vy  = pvptr[1] +   nx*nz;   /* center on plane 1 */
vz  = pvptr[2] + 2*nx*nz;
txx = pvptr[2] + 3*nx*nz;
tzz = pvptr[2] + 5*nx*nz;
txz = pvptr[2] + 7*nx*nz;

tyy0 = pvptr[0] + 4*nx*nz;
tyy1 = pvptr[1] + 4*nx*nz;
tyy2 = pvptr[2] + 4*nx*nz;
tyy3 = pvptr[3] + 4*nx*nz;

txy0 = pvptr[0] + 6*nx*nz;
txy1 = pvptr[1] + 6*nx*nz;
txy2 = pvptr[2] + 6*nx*nz;
txy3 = pvptr[3] + 6*nx*nz;

tyz0 = pvptr[0] + 8*nx*nz;
tyz1 = pvptr[1] + 8*nx*nz;
tyz2 = pvptr[2] + 8*nx*nz;
tyz3 = pvptr[3] + 8*nx*nz;

/* retain previous time step on inner ring for abs() */

vxi = pbnd + 3*nx*nz;
vyi = pbnd + 3*nx*nz + 2*(nx+nz);
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

store_iring(vx,vxi,vy,vyi,vz,vzi,nx,nz);

/* retain previous time step on plane if near iy=0 or iy=ny-1 */

if(iflag == YZERO || iflag == YN)
   storev(nx,nz,pvptr,pbnd,iflag);

for(iz=1;iz<nz-1;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

/* x derivative */

   taup = txx + iz*nx + 2;
   vp   = vx  + iz*nx + 1;
   med  = bxp + iz*nx + 1;

   if(iflag != SKIP_PLANE_Y2)
      diff(nx-3,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

   taup = txy1 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;

   if(iflag != SKIP_PLANE_Y1)
      diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);
   
   taup = txz + iz*nx + 1;
   vp   = vz  + iz*nx + 1;
   med  = bzp + iz*nx + 1;

   if(iz < nz-2 && iflag != SKIP_PLANE_Y2)
      diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

/* y derivative */

   tp0 = txy0 + iz*nx + 1;
   tp1 = txy1 + iz*nx + 1;
   tp2 = txy2 + iz*nx + 1;
   tp3 = txy3 + iz*nx + 1;
   vp  = vx   + iz*nx + 1;
   med = bxp  + iz*nx + 1;

   ordx = ord;
   if(iflag == YN) /* true -> near y=ny boundary */
      ordx = 2;

   if(iflag != SKIP_PLANE_Y2)
      diff(nx-3,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);
   
   tp0  = tyy0 + iz*nx + 1;
   tp1  = tyy1 + iz*nx + 1;
   tp2  = tyy2 + iz*nx + 1;
   tp3  = tyy3 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;
   
   ordx = ord;
   if(iflag == YZERO || iflag == YN) /* true -> near y=0 or y=ny boundary */
      ordx = 2;

   if(iflag != SKIP_PLANE_Y1)
      diff(nx-2,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);

   tp0 = tyz0 + iz*nx + 1;
   tp1 = tyz1 + iz*nx + 1;
   tp2 = tyz2 + iz*nx + 1;
   tp3 = tyz3 + iz*nx + 1;
   vp  = vz   + iz*nx + 1;
   med = bzp  + iz*nx + 1;
   
   ordx = ord;
   if(iz < nz-2 && iflag != SKIP_PLANE_Y2)
      {
      if(iflag == YN) /* true -> near y=ny boundary */
         ordx = 2;

      diff(nx-2,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);
      }
   
/* z derivative */

   if(iz == 1 || iz == nz-2)
      ordx = 2;
   else
      ordx = ord;

   taup = txz + iz*nx + 1;
   vp   = vx  + iz*nx + 1;
   med  = bxp + iz*nx + 1;
   
   if(iflag != SKIP_PLANE_Y2)
      diff(nx-3,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);

   taup = tyz1 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;
   
   if(iflag != SKIP_PLANE_Y1)
      diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);

   taup = tzz + (iz+1)*nx + 1;
   vp   = vz  +     iz*nx + 1;
   med  = bzp +     iz*nx + 1;

   if(iz < nz-2 && iflag != SKIP_PLANE_Y2)
      {
      if(iz == 1 && fs)
         ordx = ord;
      if(iz == nz-3)
         ordx = 2;

      diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
      }
   }

/* apply absorbing boundary along x and z edges */
abs_xzbnd(nx,nz,pvptr,pbnd,&fdc->dtoh,medptr,iflag,fs,intfac);

/* apply boundary condition at iy=0,ny-1 */
if(iflag == YZERO || iflag == YN)
   tstepvbnd(nx,nz,pvptr,medptr,pbnd,fdc,fs,iflag,intfac);

if(fs)    /* do freesurface for vz; vx, vy not needed */
   fsurf(nx,nz,pvptr,medptr,iflag,VEL);
   /*
   fsurfSHEAR(nx,nz,pvptr,medptr,iflag,VEL);
   */
}

/*
      fsurf() calculates velocities and stresses at and just above the free-
      surface, acording to the stress free boundary condition txz=tyz=tzz=0.

      fsurf() with iflag=VEL calculates vz just above the free-surface using
      vx and vy with the condition tzz=0 at the free-surface.  The free-surface
      is set at iz=1, coincident with the node at which txx, tyy, and tzz are
      defined.  The value for vz just above the free-surface is not needed at
      iy=0, but is needed for the calculation of txx and tyy at iy=ny-1.  Vx
      and vy are not needed above the free-surface.

      fsurf() with iflag=TAU sets tzz=0 at the free-surface (iz=1) and sets
      tzz, txz and tyz anti-symmetric above the boundary (iz=0).  The values
      for txx, tyy, and txy are not needed above the free-surface.
*/

fsurf(nx,nz,pvptr,medp,bndflag,iflag)
float **medp, **pvptr;
int nx, nz, bndflag, iflag;
{
float *vx, *vxp, *vy, *vyp, *vz, *vzp;
float *tzz, *tzzp, *txz, *txzp, *tyz, *tyzp;
float *lp, *l2mp, *mup;
int n1, ip, pend;
float a0, a1, sum, ff;

float two = 2.0;

a0 = 0.5;
a1 = 1.0;

a0 = 0.0;
a1 = 1.0;

sum = 1.0/(a0 + a1);

if(iflag == VEL)
   {
   pend = 2;        /* do plane 1 for all iy */

   if(bndflag == YN)  /* do planes 1, 2 and 3 at iy=ny-1 boundary */
      pend = 4;

   for(ip=1;ip<pend;ip++)
      {
      l2mp = medp[ip] +            nx + 1;
      lp   = medp[ip] +    nx*nz + nx + 1;
      /*
      lp   = medp[ip] + 11*nx*nz + nx + 1;
      mup  = medp[ip] + 12*nx*nz + nx + 1;
      */

      vx  = pvptr[ip]   +           nx;
      vxp = pvptr[ip]   +           nx + 1;
      vy  = pvptr[ip-1] +   nx*nz + nx + 1;
      vyp = pvptr[ip]   +   nx*nz + nx + 1;
      vz  = pvptr[ip]   + 2*nx*nz +      1;
      vzp = pvptr[ip]   + 2*nx*nz + nx + 1;

      n1 = nx - 1;
      while(n1--)
         {
         ff = sum*(a0 + a1*lp[0]/l2mp[0]);
	 /*
         ff = sum*(a0 + a1*lp[0]/(lp[0]+two*mup[0]));
	 */

         vz[0] = vzp[0] + ff*(vxp[0] - vx[0] + vyp[0] - vy[0]);

         vx++; vxp++;
         vy++; vyp++;
         vz++; vzp++;
         lp++; l2mp++; mup++;
         }

/* NEW 11/20/96 */
/* let vx, vy =0.0 above free-surface -> seems to work just as well
      vx  = pvptr[ip];
      vxp = pvptr[ip]   +           2*nx;
      vz  = pvptr[ip]   + 2*nx*nz +        1;
      vzp = pvptr[ip]   + 2*nx*nz +   nx + 1;

      n1 = nx - 1;
      while(n1--)
         {
         vx[0] = vxp[0] + vzp[0] - vzp[-1] + vz[0] - vz[-1];

         vx++; vxp++;
         vz++; vzp++;
         }

      vy  = pvptr[ip-1];
      vyp = pvptr[ip-1] +           2*nx;
      vz  = pvptr[ip-1] + 2*nx*nz;
      vzp = pvptr[ip]   + 2*nx*nz;

      n1 = nx - 1;
      while(n1--)
         {
         vy[0] = vyp[0] + vzp[0] - vz[0] + vzp[nx] - vz[nx];

         vy++; vyp++;
         vz++; vzp++;
         }
*/
/* END NEW 11/20/96 */
      }
   }

if(iflag == TAU)
   {
   pend = 2;        /* do plane 1 for all iy */

   if(bndflag == YN)  /* do planes 1 and 2 at iy=ny-1 boundary */
      pend = 3;

   for(ip=1;ip<pend;ip++)
      {
      n1 = nx - 2;

      tzz  = pvptr[ip] + 5*nx*nz +   nx + 1;
      zero(tzz,n1);

      tzz  = pvptr[ip] + 5*nx*nz +        1;
      tzzp = pvptr[ip] + 5*nx*nz + 2*nx + 1;
      txz  = pvptr[ip] + 7*nx*nz +        1;
      txzp = pvptr[ip] + 7*nx*nz +   nx + 1;
      tyz  = pvptr[ip] + 8*nx*nz +        1;
      tyzp = pvptr[ip] + 8*nx*nz +   nx + 1;

      copy_sgn(tzzp,tzz,n1,-1);
      copy_sgn(txzp,txz,n1,-1);
      copy_sgn(tyzp,tyz,n1,-1);
      }
   }
}

fsurfSHEAR(nx,nz,pvptr,medp,bndflag,iflag)
float **medp, **pvptr;
int nx, nz, bndflag, iflag;
{
float *vx, *vxp, *vy, *vyp, *vz, *vzp;
float *vz1, *vz2, *vz3, *vz4, *vz5, *vz6;
float *tzz, *tzzp, *txz, *txzp, *tyz, *tyzp;
float *lp, *l2mp;
int n1, ip, pstr, pend;
float a0, a1, sum, ff;

float two = 2.0;
float four = 4.0;

if(iflag == VEL)
   {
   pstr = 2;
   pend = 3;

   if(bndflag == YZERO)
      {
      pstr = 1;
      pend = 3;
      }

   if(bndflag == YN)
      {
      pstr = 2;
      pend = 4;
      }

   for(ip=pstr;ip<pend;ip++)
      {
      vx  = pvptr[ip]   +           nx;
      vxp = pvptr[ip]   +         2*nx;
      vz  = pvptr[ip]   + 2*nx*nz + nx;
      vzp = pvptr[ip]   + 2*nx*nz + nx + 1;

      n1 = nx - 1;
      while(n1--)
         {
         vx[0] = vxp[0] + vzp[0] - vz[0];

         vx++; vxp++;
         vz++; vzp++;
         }

      vy  = pvptr[ip-1] +   nx*nz +   nx;
      vyp = pvptr[ip-1] +   nx*nz + 2*nx;
      vz  = pvptr[ip-1] + 2*nx*nz +   nx;
      vzp = pvptr[ip]   + 2*nx*nz +   nx;

      n1 = nx;
      while(n1--)
         {
         vy[0] = vyp[0] + vzp[0] - vz[0];

         vy++; vyp++;
         vz++; vzp++;
         }

/* vz above boundary */

      vx  = pvptr[ip-1] +           2*nx;
      vxp = pvptr[ip-1] +           2*nx + 1;
      vy  = pvptr[ip-2] +   nx*nz + 2*nx + 1;
      vyp = pvptr[ip-1] +   nx*nz + 2*nx + 1;
      vz  = pvptr[ip-1] + 2*nx*nz +      + 1;

      vz1 = pvptr[ip-1] + 2*nx*nz + 2*nx + 1;
      vz2 = pvptr[ip-1] + 2*nx*nz +   nx + 2;
      vz3 = pvptr[ip-1] + 2*nx*nz +   nx;
      vz4 = pvptr[ip]   + 2*nx*nz +   nx;
      vz5 = pvptr[ip-2] + 2*nx*nz +   nx;
      vz6 = pvptr[ip-1] + 2*nx*nz +   nx + 1;

      l2mp = medp[ip] +            nx + 1;
      lp   = medp[ip] +    nx*nz + nx + 1;

      n1 = nx-2;
      while(n1--)
         {
	 ff = lp[0]/l2mp[0];

         vz[0] = vz1[0] - ff*(two*(vxp[0] - vx[0] + vyp[0] - vy[0])
	             - vz2[0] - vz3[0] - vz4[0] - vz5[0] + four*vz6[0]);

         lp++; l2mp++;
         vx++; vxp++;
         vy++; vyp++;
         vz++; vz1++; vz2++; vz3++; vz4++; vz5++; vz6++;
         }
      }
   }

if(iflag == TAU)
   {
   pstr = 2;
   pend = 3;

   if(bndflag == YZERO)
      {
      pstr = 1;
      pend = 3;
      }

   if(bndflag == YN)
      {
      pstr = 2;
      pend = 4;
      }

   for(ip=pstr;ip<pend;ip++)
      {
      txz  = pvptr[ip]   + 7*nx*nz +   nx;
      zero(txz,nx);

      tyz  = pvptr[ip-1] + 8*nx*nz +   nx;
      zero(tyz,nx);

      tzz  = pvptr[ip]   + 5*nx*nz +   nx;
      tzzp = pvptr[ip]   + 5*nx*nz + 2*nx;
      copy_sgn(tzzp,tzz,nx,-1);

      tzz  = pvptr[ip]   + 5*nx*nz;
      tzzp = pvptr[ip]   + 5*nx*nz + 3*nx;
      copy_sgn(tzzp,tzz,nx,-1);

      txz  = pvptr[ip]   + 7*nx*nz;
      txzp = pvptr[ip]   + 7*nx*nz + 2*nx;
      copy_sgn(txzp,txz,nx,-1);

      tyz  = pvptr[ip-1] + 8*nx*nz;
      tyzp = pvptr[ip-1] + 8*nx*nz + 2*nx;
      copy_sgn(tyzp,tyz,nx,-1);
      }
   }
}

/*
	tstepvP3() updates the velocity variables vx, vy, and vz for the new
	time step.  The varibles vx and vz are centered in plane 2 and the
	variable vy is centered in plane 1.
	
	These variables are absorbed
	at the x, y and z boundaries.

	if iflag=YZERO -> 2nd order y derivative of tyy, affects vy near iy=0;
                          also calls function tstepvbndP3() to update velocities
                          at iy=0 boundary.
	if iflag=YN    -> 2nd order y derivative of tyy, txy and txz, affects
                          vx, vy and vz near iy=ny-1;  also calls function
			  tstepvbndP3() to update vx and vz at iy=ny-1 boundary.
*/

void tstepvP3(int nx,int nz,float **pvptr,float **medptr,float *pbnd,struct fdcoefs *fdc,int fs,int ybndflag,struct interface *intfac,struct nodeinfo *ni)
{
float *vx, *vy, *vz, *txx, *tzz, *txz;
float *tyy0, *tyy1, *tyy2, *tyy3;
float *txy0, *txy1, *txy2, *txy3;
float *tyz0, *tyz1, *tyz2, *tyz3;
float *vxi, *vyi, *vzi;
float *taup, *vp, *med, *med1, *med2;
float *tp0, *tp1, *tp2, *tp3;
float *bxp, *byp, *bzp;
int ix, iz, ord, ordx, off1, off2, off3;
int ixend, izstart, izend;

float c0, c1;
 
/* do interior solution */

bxp  = medptr[2] + 5*nx*nz; /* center on plane 2 for vx */
byp  = medptr[1] + 6*nx*nz; /* center on plane 1 for vy */
bzp  = medptr[2] + 7*nx*nz; /* center on plane 2 for vz */

vx  = pvptr[2];
vy  = pvptr[1] +   nx*nz;   /* center on plane 1 */
vz  = pvptr[2] + 2*nx*nz;
txx = pvptr[2] + 3*nx*nz;
tzz = pvptr[2] + 5*nx*nz;
txz = pvptr[2] + 7*nx*nz;

tyy0 = pvptr[0] + 4*nx*nz;
tyy1 = pvptr[1] + 4*nx*nz;
tyy2 = pvptr[2] + 4*nx*nz;
tyy3 = pvptr[3] + 4*nx*nz;

txy0 = pvptr[0] + 6*nx*nz;
txy1 = pvptr[1] + 6*nx*nz;
txy2 = pvptr[2] + 6*nx*nz;
txy3 = pvptr[3] + 6*nx*nz;

tyz0 = pvptr[0] + 8*nx*nz;
tyz1 = pvptr[1] + 8*nx*nz;
tyz2 = pvptr[2] + 8*nx*nz;
tyz3 = pvptr[3] + 8*nx*nz;

/* retain previous time step on inner ring for abs() */

vxi = pbnd + 3*nx*nz;
vyi = pbnd + 3*nx*nz + 2*(nx+nz);
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

store_iring(vx,vxi,vy,vyi,vz,vzi,nx,nz);

/* retain previous time step on plane if near iy=0 or iy=ny-1 */

if(ybndflag == YZERO || ybndflag == YN)
   storev(nx,nz,pvptr,pbnd,ybndflag);

ixend = nx-2;
if(ni->plusId_x == -1)
   ixend = nx-3;

izstart = 2;
if(ni->minusId_z == -1)
   izstart = 1;

izend = nz-2;
if(ni->plusId_z == -1)
   izend = nz-1;

for(iz=izstart;iz<izend;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

/* x derivative */

   taup = txx + iz*nx + 2;
   vp   = vx  + iz*nx + 1;
   med  = bxp + iz*nx + 1;

   if(ybndflag != SKIP_PLANE_Y2)
      diff(ixend,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

   taup = txy1 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;

   if(ybndflag != SKIP_PLANE_Y1)
      diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);
   
   taup = txz + iz*nx + 1;
   vp   = vz  + iz*nx + 1;
   med  = bzp + iz*nx + 1;

   if(iz < nz-2 && ybndflag != SKIP_PLANE_Y2)
      diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

/* y derivative */

   tp0 = txy0 + iz*nx + 1;
   tp1 = txy1 + iz*nx + 1;
   tp2 = txy2 + iz*nx + 1;
   tp3 = txy3 + iz*nx + 1;
   vp  = vx   + iz*nx + 1;
   med = bxp  + iz*nx + 1;

   ordx = ord;
   if(ybndflag == YN) /* true -> near y=ny boundary */
      ordx = 2;

   if(ybndflag != SKIP_PLANE_Y2)
      diff(ixend,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);
   
   tp0  = tyy0 + iz*nx + 1;
   tp1  = tyy1 + iz*nx + 1;
   tp2  = tyy2 + iz*nx + 1;
   tp3  = tyy3 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;
   
   ordx = ord;
   if(ybndflag == YZERO || ybndflag == YN) /* true -> near y=0 or y=ny boundary */
      ordx = 2;

   if(ybndflag != SKIP_PLANE_Y1)
      diff(nx-2,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);

   tp0 = tyz0 + iz*nx + 1;
   tp1 = tyz1 + iz*nx + 1;
   tp2 = tyz2 + iz*nx + 1;
   tp3 = tyz3 + iz*nx + 1;
   vp  = vz   + iz*nx + 1;
   med = bzp  + iz*nx + 1;
   
   ordx = ord;
   if(iz < nz-2 && ybndflag != SKIP_PLANE_Y2)
      {
      if(ybndflag == YN) /* true -> near y=ny boundary */
         ordx = 2;

      diff(nx-2,vp,med,tp2,tp1,tp3,tp0,fdc,ordx,YDERIV);
      }
   
/* z derivative */

/* RWG 2012-07-04

For Txz and Tyz, apply fourth-order Dz derivatives at iz=1.  Old way used 2nd-order, which works
quite well for forward calculations.  But when force is inserted at free-surface (at Vx
or Vy node) this results in a slight amplification relative to forward calculation (about 10%).
This was discovered using SGT tests for CyberShake.  Using fourth-order operators for the
derivatives on the shear stress components reduces the misfit to a couple of percent.

The fourth-order derivatives are computed using anti-symmetry for the shear stress values
above the free-surface:

Txz(k-1) = -Txz(k),     Tyz(k-1) = -Tyz(k)
Txz(k-2) = -Txz(k+1),   Tyz(k-2) = -Tyz(k+1)

*/

   if(iz == 1 || iz == nz-2)
      ordx = 2;
   else
      ordx = ord;

   taup = txz + iz*nx + 1;
   vp   = vx  + iz*nx + 1;
   med  = bxp + iz*nx + 1;

   if(ybndflag != SKIP_PLANE_Y2)
      {
      if(iz == 1)
         {
         c0 = 2.0*fdc->dtoh;
         c1 = 0.0;
         if(ord == 4)
            {
            c0 = 2.0*fdc->c0;
            c1 = 2.0*fdc->c1;
            }

         for(ix=0;ix<ixend;ix++)
            vp[ix] = vp[ix] + med[ix]*(c0*taup[ix] + c1*taup[ix+nx]);
         }
      else
         diff(ixend,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
      }

   taup = tyz1 + iz*nx + 1;
   vp   = vy   + iz*nx + 1;
   med  = byp  + iz*nx + 1;
   
   if(ybndflag != SKIP_PLANE_Y1)
      {
      if(iz == 1)
         {
         c0 = 2.0*fdc->dtoh;
         c1 = 0.0;
         if(ord == 4)
            {
            c0 = 2.0*fdc->c0;
            c1 = 2.0*fdc->c1;
            }

         for(ix=0;ix<nx-2;ix++)
            vp[ix] = vp[ix] + med[ix]*(c0*taup[ix] + c1*taup[ix+nx]);
         }
      else
         diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
      }

   taup = tzz + (iz+1)*nx + 1;
   vp   = vz  +     iz*nx + 1;
   med  = bzp +     iz*nx + 1;

   if(iz < nz-2 && ybndflag != SKIP_PLANE_Y2)
      {
      if(ni->minusId_z == -1 && iz == 1 && fs)
         ordx = ord;
      if(ni->plusId_z == -1 && iz == nz-3)
         ordx = 2;

      diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
      }
   }

/* apply absorbing boundary along x and z edges */
abs_xzbndP3(nx,nz,pvptr,pbnd,&fdc->dtoh,medptr,ybndflag,fs,intfac,ni);

/* apply boundary condition at iy=0,ny-1 */
if(ybndflag == YZERO || ybndflag == YN)
   tstepvbndP3(nx,nz,pvptr,medptr,pbnd,fdc,fs,ybndflag,intfac,ni);

/* do freesurface for vz; vx, vy not needed */
/* 2012-07-31 do this as separate call from main_mpi.c
if(ni->minusId_z == -1 && fs)
   fsurf(nx,nz,pvptr,medptr,ybndflag,VEL);
*/
   /*
   fsurfSHEAR(nx,nz,pvptr,medptr,ybndflag,VEL);
   */
}

/*
	When iflag=YZERO (near iy=0 boundary), tstepvbndP3() updates the velocity
        variables vx and vz for the new time step on plane 1 using 2nd order
        y derivatives.
      
	The values of vx, vy and vz on plane 0 (iflag=YZERO) or plane 3
        (iflag=YN) are updated using an absorbing boundary condition in the
        function abs_ybndP3().
*/

void tstepvbndP3(int nx,int nz,float **pvptr,float **medptr,float *pbnd,struct fdcoefs *fdc,int fs,int ybndflag,struct interface *intfac,struct nodeinfo *ni)
{
float *vx, *vy, *vz, *txx, *tzz, *txz;
float *txy0, *txy1;
float *tyz0, *tyz1;
float *tyy2, *tyy3;
float *txy2;
float *tyz2;
float *vxi, *vyi, *vzi;
float *taup, *vp, *med;
float *tp0, *tp1, *tp2, *tp3;
float *bxp, *byp, *bzp;
int ix, iz, ord, ordx, off1, off2, off3, ip;
int ixend, izstart, izend;
 
/*
         ybndflag=YZERO -> do interior solution for vx and vz; vy already done in tstepv().
         ybndflag=YN    -> go to abs_ybnd(), vx,vy and vz already done in tstepv().
*/

if(ybndflag == YZERO) /* set pointers for vx and vz update */
   {
/* vy, vyi not solved, but needed for store_iring() */
   bxp = medptr[1] + 5*nx*nz;
   byp = medptr[1] + 6*nx*nz;
   bzp = medptr[1] + 7*nx*nz;

   vx  = pvptr[1];
   vy  = pvptr[1] +   nx*nz;
   vz  = pvptr[1] + 2*nx*nz;
   txx = pvptr[1] + 3*nx*nz;
   tzz = pvptr[1] + 5*nx*nz;
   txz = pvptr[1] + 7*nx*nz;

   txy0 = pvptr[0] + 6*nx*nz;
   txy1 = pvptr[1] + 6*nx*nz;

   tyz0 = pvptr[0] + 8*nx*nz;
   tyz1 = pvptr[1] + 8*nx*nz;
 
   /* retain previous time step on inner ring for abs() */

   vxi = pbnd + 3*nx*nz;
   vyi = pbnd + 3*nx*nz + 2*(nx+nz);
   vzi = pbnd + 3*nx*nz + 4*(nx+nz);
   store_iring(vx,vxi,vy,vyi,vz,vzi,nx,nz);

   ixend = nx-2;
   if(ni->plusId_x == -1)
      ixend = nx-3;

   izstart = 2;
   if(ni->minusId_z == -1)
      izstart = 1;

   izend = nz-2;
   if(ni->plusId_z == -1)
      izend = nz-1;

   for(iz=izstart;iz<izend;iz++)
      {
      if(iz < fdc->izord2)
         ord = fdc->order;
      else
         ord = fdc->ordx;

/* x derivative */

      taup = txx + iz*nx + 2;
      vp   = vx  + iz*nx + 1;
      med  = bxp + iz*nx + 1;

      diff(ixend,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);

      taup = txz + iz*nx + 1;
      vp   = vz  + iz*nx + 1;
      med  = bzp + iz*nx + 1;

      if(iz < nz-2)
         diff(nx-2,vp,med,taup,taup-1,taup+1,taup-2,fdc,ord,XDERIV);
   
/* 2nd order y derivative */

      tp0 = txy0 + iz*nx + 1;
      tp1 = txy1 + iz*nx + 1;
      vp  = vx   + iz*nx + 1;
      med = bxp  + iz*nx + 1;

      diff(ixend,vp,med,tp1,tp0,taup,taup,fdc,2,YDERIV);

      tp0 = tyz0 + iz*nx + 1;
      tp1 = tyz1 + iz*nx + 1;
      vp  = vz   + iz*nx + 1;
      med = bzp  + iz*nx + 1;
   
      if(iz < nz-2)
         diff(nx-2,vp,med,tp1,tp0,taup,taup,fdc,2,YDERIV);
   
/* z derivative */

      if(iz == 1 || iz == nz-2)
         ordx = 2;
      else
         ordx = ord;

      taup = txz + iz*nx + 1;
      vp   = vx  + iz*nx + 1;
      med  = bxp + iz*nx + 1;

      diff(ixend,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);

      taup = tzz + (iz+1)*nx + 1;
      vp   = vz  +     iz*nx + 1;
      med  = bzp +     iz*nx + 1;

      if(iz < nz-2)
         {
      /* only have 4th-order information above free-surface for tzz */
         if(ni->minusId_z == -1 && iz == 1 && fs)
            ordx = ord;
         if(ni->plusId_z == -1 && iz == nz-3)
            ordx = 2;

         diff(nx-2,vp,med,taup,taup-nx,taup+nx,taup-2*nx,fdc,ordx,ZDERIV);
	 }
      }

   /* apply absorbing boundary along x and z edges */
   abs_xzbnd1P3(nx,nz,pvptr,pbnd,&fdc->dtoh,medptr,fs,intfac,ni);
   }

/* absorb vx and vz on plane at iy=0,ny-1 */
/* absorb vy on plane at iy=0,ny-2 */
abs_ybndP3(nx,nz,pvptr,medptr,pbnd,fdc,fs,ybndflag,intfac,ni);
}

/*
	tsteppP3() updates the stress variables txx, tyy, tzz, txy, txz and tyz
	for the new time step.  The varibles txx, tyy, tzz and txz are centered
	in plane 2 and the variables txy and tyz are centered in plane 1.
*/

void tsteppP3(int nx,int nz,float **pvptr,float **medptr,struct fdcoefs *fdc,int fs,int ybndflag,int eflag,struct nodeinfo *ni)
{
float *vx0, *vx1, *vx2, *vx3;
float *vy0, *vy1, *vy2, *vy3;
float *vz0, *vz1, *vz2, *vz3;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *taup, *vp, *med;
float *tp1, *tp2, *tp3, *med1, *med2;
float *cp1, *cp2, *cp3, *pkp, *skp;
float *vp0, *vp1, *vp2, *vp3;
float *l2mp, *lp, *mxyp, *mxzp, *myzp;
float *ak, *pk, *sk, *aky, *sky;
int ix, iz, ord, ordx, j;
int izstart, izend;

l2mp  = medptr[2];
lp    = medptr[2] +   nx*nz;
mxyp  = medptr[1] + 2*nx*nz; /* centered on plane 1 for txy */
mxzp  = medptr[2] + 3*nx*nz; /* centered on plane 2 for txz */
myzp  = medptr[1] + 4*nx*nz; /* centered on plane 1 for tyz */

ak    = medptr[2] +  8*nx*nz;
aky   = medptr[1] +  8*nx*nz; /* centered on plane 1 for txy,tyz */
pk    = medptr[2] +  9*nx*nz;
sk    = medptr[2] + 10*nx*nz;
sky   = medptr[1] + 10*nx*nz; /* centered on plane 1 for txy,tyz */
 
vx0 = pvptr[0];
vx1 = pvptr[1];
vx2 = pvptr[2];
vx3 = pvptr[3];

vy0 = pvptr[0] +   nx*nz;
vy1 = pvptr[1] +   nx*nz;
vy2 = pvptr[2] +   nx*nz;
vy3 = pvptr[3] +   nx*nz;

vz0 = pvptr[0] + 2*nx*nz;
vz1 = pvptr[1] + 2*nx*nz;
vz2 = pvptr[2] + 2*nx*nz;
vz3 = pvptr[3] + 2*nx*nz;

txx = pvptr[2] + 3*nx*nz;
tyy = pvptr[2] + 4*nx*nz;
tzz = pvptr[2] + 5*nx*nz;
txy = pvptr[1] + 6*nx*nz; /* centered on plane 1 */
txz = pvptr[2] + 7*nx*nz;
tyz = pvptr[1] + 8*nx*nz; /* centered on plane 1 */

cxx = pvptr[2] +  9*nx*nz;
cyy = pvptr[2] + 10*nx*nz;
czz = pvptr[2] + 11*nx*nz;
cxy = pvptr[1] + 12*nx*nz; /* centered on plane 1 */
cxz = pvptr[2] + 13*nx*nz;
cyz = pvptr[1] + 14*nx*nz; /* centered on plane 1 */

if(eflag == 0)
   {
   if(ybndflag != SKIP_PLANE_Y2)
      {
      mv_mult1(txx,cxx,ak,nx*nz);
      mv_mult1(tyy,cyy,ak,nx*nz);
      mv_mult1(tzz,czz,ak,nx*nz);
      mv_mult1(txz,cxz,ak,nx*nz);
      }
   if(ybndflag != SKIP_PLANE_Y1)
      {
      mv_mult1(txy,cxy,aky,nx*nz);
      mv_mult1(tyz,cyz,aky,nx*nz);
      }
   }

izstart = 2;
if(ni->minusId_z == -1)
   izstart = 1;

izend = nz-1;
if(ni->plusId_z == -1)
   izend = nz;

for(iz=izstart;iz<izend;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

   ordx = ord;

/* x derivative */

   if(ybndflag != SKIP_PLANE_Y2)
      {
      tp1  = txx  + iz*nx + 1;
      tp2  = tyy  + iz*nx + 1;
      tp3  = tzz  + iz*nx + 1;
      vp   = vx2  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cxx  + iz*nx + 1;
      cp2  = cyy  + iz*nx + 1;
      cp3  = czz  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
         {
         if(ni->minusId_z == -1 && iz == 1 && fs)
            ordx = 2;
         else
            ordx = ord;

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-1,vp+1,vp-2,fdc,ordx,XDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-1,vp+1,vp-2,fdc,ordx,XDERIV);
         }

      taup = txz  + (iz-1)*nx;
      vp   = vz2  + (iz-1)*nx + 1;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }

   if(iz < (nz-1) && ybndflag != SKIP_PLANE_Y1) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp   = vy1  + iz*nx + 1;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sky  + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }
  
/* y derivative */

   if(ybndflag != SKIP_PLANE_Y2)
      {
      tp1  = tyy  + iz*nx + 1;
      tp2  = tzz  + iz*nx + 1;
      tp3  = txx  + iz*nx + 1;
      vp0  = vy0  + iz*nx + 1;
      vp1  = vy1  + iz*nx + 1;
      vp2  = vy2  + iz*nx + 1;
      vp3  = vy3  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cyy  + iz*nx + 1;
      cp2  = czz  + iz*nx + 1;
      cp3  = cxx  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
         {
         if((ni->minusId_z == -1 && iz == 1 && fs) || ybndflag == YN)
            ordx = 2;
         else
            ordx = ord;

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp2,vp1,vp3,vp0,fdc,ordx,YDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp2,vp1,vp3,vp0,fdc,ordx,YDERIV);
         }
      }

   if(ybndflag != SKIP_PLANE_Y1)
      {
      taup = tyz  + (iz-1)*nx + 1;
      vp0  = vz0  + (iz-1)*nx + 1;
      vp1  = vz1  + (iz-1)*nx + 1;
      vp2  = vz2  + (iz-1)*nx + 1;
      vp3  = vz3  + (iz-1)*nx + 1;
      med  = myzp + (iz-1)*nx + 1;

      cp1  = cyz  + (iz-1)*nx + 1;
      skp  = sky  + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
      if(eflag == 0)
         diff_mv(nx-2,taup,cp1,med,skp,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
      else
         diff(nx-2,taup,med,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);

      if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
         {
         taup = txy  + iz*nx;
         vp0  = vx0  + iz*nx;
         vp1  = vx1  + iz*nx;
         vp2  = vx2  + iz*nx;
         vp3  = vx3  + iz*nx;
         med  = mxyp + iz*nx;

         cp1  = cxy  + iz*nx;
         skp  = sky  + iz*nx;

         if(eflag == 0)
            diff_mv(nx-1,taup,cp1,med,skp,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
         else
            diff(nx-1,taup,med,vp2,vp1,vp3,vp0,fdc,ord,YDERIV);
         }
      }
  
/* do x-strips: z derivative */

   if(ybndflag != SKIP_PLANE_Y2)
      {
      tp1  = tzz  + iz*nx + 1;
      tp2  = txx  + iz*nx + 1;
      tp3  = tyy  + iz*nx + 1;
      vp   = vz2  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = czz  + iz*nx + 1;
      cp2  = cxx  + iz*nx + 1;
      cp3  = cyy  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1)
         {
         ordx = ord;
         if(iz == 1 || iz == nz-2)
            ordx = 2;

vp2 = vp-2*nx;
/* XXXXYYYY
*/
if(iz == 1)
   {
   ordx = ord;
   vp2 = vz2  + (nz-1)*nx + 1;
   }

         if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-nx,vp+nx,vp2,fdc,ordx,ZDERIV);
         else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-nx,vp+nx,vp2,fdc,ordx,ZDERIV);
         }

      ordx = ord;
      if(iz == 1 || iz == nz-1)
         ordx = 2;

      taup = txz  + (iz-1)*nx;
      vp   = vx2  +     iz*nx;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }

   if(ybndflag != SKIP_PLANE_Y1)
      {
      ordx = ord;
      if(iz == 1 || iz == nz-1)
         ordx = 2;

      taup = tyz  + (iz-1)*nx + 1;
      vp   = vy1  +     iz*nx + 1;
      med  = myzp + (iz-1)*nx + 1;

      cp1  = cyz  + (iz-1)*nx + 1;
      skp  = sky  + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
      if(eflag == 0)
         diff_mv(nx-2,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-2,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }
   }

if(eflag == 0)
   {
   if(ybndflag != SKIP_PLANE_Y2)
      {
      mv_mult2(txx,cxx,nx*nz);
      mv_mult2(tyy,cyy,nx*nz);
      mv_mult2(tzz,czz,nx*nz);
      mv_mult2(txz,cxz,nx*nz);
      }
   if(ybndflag != SKIP_PLANE_Y1)
      {
      mv_mult2(txy,cxy,nx*nz);
      mv_mult2(tyz,cyz,nx*nz);
      }
   }

 /* apply 2nd order boundary condition at iy=0,ny-1 */
if(ybndflag == YZERO || ybndflag == YN)
   tsteppbndP3(nx,nz,pvptr,medptr,fdc,fs,ybndflag,eflag,ni);

if(ni->minusId_z == -1 && fs)    /* do freesurface for tzz, txz, and tyz */
   fsurf(nx,nz,pvptr,medptr,ybndflag,TAU);
   /*
   fsurfSHEAR(nx,nz,pvptr,medptr,ybndflag,TAU);
   */
}

/*
	tsteppbndP3() updates the stress variables txx, tyy, tzz, txy, txz
        and tyz for the new time step when iy=0 or iy=ny-1.  The varibles
        txx, tyy, tzz and txz are centered in plane 1,3 when
        iy=0,ny-1; respectively, and the variables txy and tyz are centered
        in plane 0,2 when iy=0,ny-1; respectively.  The field txz is not needed
        when iy=ny-1.
*/

void tsteppbndP3(int nx,int nz,float **pvptr,float **medptr,struct fdcoefs *fdc,int fs,int ybndflag,int eflag,struct nodeinfo *ni)
{
float *vx0, *vx1;
float *vy0, *vy1;
float *vz0, *vz1;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *taup, *vp, *med;
float *tp1, *tp2, *tp3, *med1, *med2;
float *cp1, *cp2, *cp3, *pkp, *skp;
float *vp0, *vp1;
float *l2mp, *lp, *mxyp, *mxzp, *myzp;
float *ak, *pk, *sk;
int ix, iz, ord, ordx, ip, ipm;
int ixend, izstart, izend;

/*
   ip = 0 when ybndflag = YZERO (iy=0) -> solve all.
   ip = 2 when ybndflag = YN (iy=ny-1) -> solve only txy and tyz.

   ipm = 1 when ybndflag = YZERO
   ipm = 0 when ybndflag = YN
           => this shifts media pointer to be 1 plane inside
	      boundary (consistent with absorbing boundary condition)
*/

if(ybndflag == YZERO)
   {
   ipm = 1;
   ip = 0;
   }
else if(ybndflag == YN)
   {
   ipm = 0;
   ip = 2;
   }

l2mp  = medptr[ip+ipm];
lp    = medptr[ip+ipm] +   nx*nz;
mxyp  = medptr[ip+ipm] + 2*nx*nz;
mxzp  = medptr[ip+ipm] + 3*nx*nz;
myzp  = medptr[ip+ipm] + 4*nx*nz;

ak    = medptr[ip+ipm] +  8*nx*nz;
pk    = medptr[ip+ipm] +  9*nx*nz;
sk    = medptr[ip+ipm] + 10*nx*nz;

vx0 = pvptr[ip];
vx1 = pvptr[ip+1];

vy0 = pvptr[ip]   + nx*nz;
vy1 = pvptr[ip+1] + nx*nz;

vz0 = pvptr[ip]   + 2*nx*nz;
vz1 = pvptr[ip+1] + 2*nx*nz;

txx = pvptr[ip+1] + 3*nx*nz;
tyy = pvptr[ip+1] + 4*nx*nz;
tzz = pvptr[ip+1] + 5*nx*nz;
txy = pvptr[ip]   + 6*nx*nz; /* centered on plane 0,2 */
txz = pvptr[ip+1] + 7*nx*nz;
tyz = pvptr[ip]   + 8*nx*nz; /* centered on plane 0,2 */

cxx = pvptr[ip+1] +  9*nx*nz;
cyy = pvptr[ip+1] + 10*nx*nz;
czz = pvptr[ip+1] + 11*nx*nz;
cxy = pvptr[ip]   + 12*nx*nz; /* centered on plane 0,2 */
cxz = pvptr[ip+1] + 13*nx*nz;
cyz = pvptr[ip]   + 14*nx*nz; /* centered on plane 0,2 */

if(eflag == 0)
   {
   mv_mult1(txx,cxx,ak,nx*nz);
   mv_mult1(tyy,cyy,ak,nx*nz);
   mv_mult1(tzz,czz,ak,nx*nz);
   mv_mult1(txy,cxy,ak,nx*nz);
   mv_mult1(txz,cxz,ak,nx*nz);
   mv_mult1(tyz,cyz,ak,nx*nz);
   }

izstart = 2;
if(ni->minusId_z == -1)
   izstart = 1;

izend = nz-1;
if(ni->plusId_z == -1)
   izend = nz;

for(iz=izstart;iz<izend;iz++)
   {
   if(iz < fdc->izord2)
      ord = fdc->order;
   else
      ord = fdc->ordx;

/* x derivative */

   if(ybndflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = txx  + iz*nx + 1;
      tp2  = tyy  + iz*nx + 1;
      tp3  = tzz  + iz*nx + 1;
      vp   = vx1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cxx  + iz*nx + 1;
      cp2  = cyy  + iz*nx + 1;
      cp3  = czz  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
	 {
	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
	 }
      }

   if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp   = vy0  + iz*nx + 1;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sk   + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }

   if(ybndflag == YZERO) /* txz not needed at iy = ny-1 */
      {
      taup = txz  + (iz-1)*nx;
      vp   = vz1  + (iz-1)*nx + 1;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      else
         diff(nx-1,taup,med,vp,vp-1,vp+1,vp-2,fdc,ord,XDERIV);
      }
  
/* y derivative */

   if(ybndflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = tyy  + iz*nx + 1;
      tp2  = tzz  + iz*nx + 1;
      tp3  = txx  + iz*nx + 1;
      vp0  = vy0  + iz*nx + 1;
      vp1  = vy1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = cyy  + iz*nx + 1;
      cp2  = czz  + iz*nx + 1;
      cp3  = cxx  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1) /* txx,tyy,tzz not needed at iz=nz-1 */
	 {
	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp1,vp0,vp,vp,fdc,2,YDERIV);
	 }
      }

   taup = tyz  + (iz-1)*nx + 1;
   vp0  = vz0  + (iz-1)*nx + 1;
   vp1  = vz1  + (iz-1)*nx + 1;
   med  = myzp + (iz-1)*nx + 1;

   cp1  = cyz  + (iz-1)*nx + 1;
   skp  = sk   + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */

   if(eflag == 0)
      diff_mv(nx-2,taup,cp1,med,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
   else
      diff(nx-2,taup,med,vp1,vp0,vp,vp,fdc,2,YDERIV);

   if(iz < (nz-1)) /* txy not needed at iz = nz-1 */
      {
      taup = txy  + iz*nx;
      vp0  = vx0  + iz*nx;
      vp1  = vx1  + iz*nx;
      med  = mxyp + iz*nx;

      cp1  = cxy  + iz*nx;
      skp  = sk   + iz*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp1,vp0,vp,vp,fdc,2,YDERIV);
      else
         diff(nx-1,taup,med,vp1,vp0,vp,vp,fdc,2,YDERIV);
      }
  
/* z derivative */

   if(ybndflag != YN) /* txx,tyy,tzz not needed at iy=ny-1 */
      {
      tp1  = tzz  + iz*nx + 1;
      tp2  = txx  + iz*nx + 1;
      tp3  = tyy  + iz*nx + 1;
      vp   = vz1  + iz*nx + 1;
      med1 = l2mp + iz*nx + 1;
      med2 = lp   + iz*nx + 1;

      cp1  = czz  + iz*nx + 1;
      cp2  = cxx  + iz*nx + 1;
      cp3  = cyy  + iz*nx + 1;
      pkp  = pk   + iz*nx + 1;
      skp  = sk   + iz*nx + 1;

      if(iz < nz-1)
         {
         if(iz == 1 || iz == nz-2)
            ordx = 2;
         else
            ordx = ord;

	 if(eflag == 0)
            diffx_mv(nx-2,tp1,tp2,tp3,cp1,cp2,cp3,med1,med2,pkp,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
	 else
            diffx(nx-2,tp1,tp2,tp3,med1,med2,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
         }
      }

   if(iz == 1 || iz == nz-1)
      ordx = 2;
   else
      ordx = ord;

   if(ybndflag == YZERO) /* txz not needed at iy = ny-1 */
      {
      taup = txz  + (iz-1)*nx;
      vp   = vx1  +     iz*nx;
      med  = mxzp + (iz-1)*nx;

      cp1  = cxz  + (iz-1)*nx;
      skp  = sk   + (iz-1)*nx;

      if(eflag == 0)
         diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      else
         diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
      }

   taup = tyz  + (iz-1)*nx + 1;
   vp   = vy0  +     iz*nx + 1;
   med  = myzp + (iz-1)*nx + 1;

   cp1  = cyz  + (iz-1)*nx + 1;
   skp  = sk   + (iz-1)*nx + 1;

                 /* tyz not needed at ix = nx-1 */
   if(eflag == 0)
      diff_mv(nx-1,taup,cp1,med,skp,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
   else
      diff(nx-1,taup,med,vp,vp-nx,vp+nx,vp-2*nx,fdc,ordx,ZDERIV);
   }

if(eflag == 0)
   {
   mv_mult2(txx,cxx,nx*nz);
   mv_mult2(tyy,cyy,nx*nz);
   mv_mult2(tzz,czz,nx*nz);
   mv_mult2(txy,cxy,nx*nz);
   mv_mult2(txz,cxz,nx*nz);
   mv_mult2(tyz,cyz,nx*nz);
   }
}

void fsurf_vel(int nx,int nz,float **pvptr,float **medp,struct fdcoefs *fdc,int bndflag)
{
float *vx0, *vx1, *vx2, *vx3, *vx4;
float *vy0, *vy1, *vy2, *vy3, *vy4;
float *vz0, *vz1, *vz2, *vz3;
float *lp, *l2mp, *mup;
int n1, ip, pstr, pend;
float a0, a1, sum, ff;
float c0xy, c0z, c1, c2xy, c2z;

float two = 2.0;

/* XXXXYYYY */

if(fdc->order == 2)
   {
   c0xy = 0.0;
   c0z = 0.0;
   }
else if(fdc->order == 4)
   {
   c0xy = -1.0/27.0;
   c0z = -1.0/27.0;
   }

c1 = 1.0/(1.0 + c0xy);
c2z = c0z*c1;
c2xy = c0xy*c1;

a0 = 0.0;
a1 = 1.0;

sum = 1.0/(a0 + a1);

ip = 2;     /* for interior portion of model */

/*
   1st, set vel components to zero at 'h' above FS, this is iz=0 for Vx and Vy.
   For Vz, we use the bottom row since it is not used for abs. boundary and would
   correspond to iz=-1.
*/

vx0 = pvptr[ip];
zero(vx0,nx);

vy0 = pvptr[ip] +   nx*nz;
zero(vy0,nx);

vz0 = pvptr[ip] + 2*nx*nz + (nz-1)*nx;
zero(vz0,nx);

/*
   Next solve for Vz at h/2 above FS.
*/

l2mp = medp[ip] +            nx + 2;
lp   = medp[ip] +    nx*nz + nx + 2;

vx0 = pvptr[ip]   +           nx;
vx1 = pvptr[ip]   +           nx + 1;
vx2 = pvptr[ip]   +           nx + 2;
vx3 = pvptr[ip]   +           nx + 3;

vy0 = pvptr[ip-2] +   nx*nz + nx + 2;
vy1 = pvptr[ip-1] +   nx*nz + nx + 2;
vy2 = pvptr[ip]   +   nx*nz + nx + 2;
vy3 = pvptr[ip+1] +   nx*nz + nx + 2;

/* vz0 is assumed to be =0.0, we solve for vz1 */
vz1 = pvptr[ip]   + 2*nx*nz +        2;
vz2 = pvptr[ip]   + 2*nx*nz +   nx + 2;
vz3 = pvptr[ip]   + 2*nx*nz + 2*nx + 2;

/* 2nd order at ix=1 boundary */
ff = sum*(a0 + a1*lp[-1]/l2mp[-1]);
vz1[-1] = vz2[-1] + ff*(vx2[-1] - vx1[-1] + vy2[-1] - vy1[-1]);

n1 = nx - 3;
while(n1--)
   {
   ff = sum*(a0 + a1*lp[0]/l2mp[0]);

   vz1[0] = vz2[0] + ff*(vx2[0] - vx1[0] + vy2[0] - vy1[0])
                      + c0z*vz3[0] + c0xy*ff*(vx3[0] - vx0[0] + vy3[0] - vy0[0]);

   vx0++; vx1++; vx2++; vx3++;
   vy0++; vy1++; vy2++; vy3++;
          vz1++; vz2++; vz3++;
   lp++; l2mp++;
   }

/* 2nd order at ix=nx-1 boundary, remember pointers already incremented */
ff = sum*(a0 + a1*lp[0]/l2mp[0]);
vz1[0] = vz2[0] + ff*(vx2[0] - vx1[0] + vy2[0] - vy1[0]);

/* 
   Next solve for Vx at h above FS (in row iz=0).
*/

/* vx0 is assumed to be =0.0, don't need vx2; we solve for vx1 */
vx1 = pvptr[ip]                    + 1;
vx3 = pvptr[ip]   +           2*nx + 1;
vx4 = pvptr[ip]   +           3*nx + 1;
vz0 = pvptr[ip]   + 2*nx*nz;
vz1 = pvptr[ip]   + 2*nx*nz        + 1;
vz2 = pvptr[ip]   + 2*nx*nz        + 2;
vz3 = pvptr[ip]   + 2*nx*nz        + 3;

/* 2nd order at ix=0 boundary */
vx1[-1] = vx3[-1] + c1*(vz2[-1] - vz1[-1] + vz2[nx-1] - vz1[nx-1]);

n1 = nx - 3;
while(n1--)
   {
   vx1[0] = vx3[0] + c1*(vz2[0] - vz1[0] + vz2[nx] - vz1[nx])
                      + c2xy*vx4[0] + c2z*(vz3[0] - vz0[0] + vz3[nx] - vz0[nx]);

          vx1++;        vx3++; vx4++;
   vz0++; vz1++; vz2++; vz3++;
   }

/* 2nd order at ix=nx-1 boundary, remember pointers already incremented */
vx1[0] = vx3[0] + c1*(vz2[0] - vz1[0] + vz2[nx] - vz1[nx]);

/*
   Next solve for Vx at h above FS (in row iz=0).
*/

/* vy0 is assumed to be =0.0, don't need vy2; we solve for vy1 */
vy1 = pvptr[ip-1] +   nx*nz        + 1;
vy3 = pvptr[ip-1] +   nx*nz + 2*nx + 1;
vy4 = pvptr[ip-1] +   nx*nz + 3*nx + 1;
vz0 = pvptr[ip-2] + 2*nx*nz        + 1;
vz1 = pvptr[ip-1] + 2*nx*nz        + 1;
vz2 = pvptr[ip]   + 2*nx*nz        + 1;
vz3 = pvptr[ip+1] + 2*nx*nz        + 1;

n1 = nx - 1;
while(n1--)   {
   vy1[0] = vy3[0] + c1*(vz2[0] - vz1[0] + vz2[nx] - vz1[nx])
                      + c2xy*vy4[0] + c2z*(vz3[0] - vz0[0] + vz3[nx] - vz0[nx]);

          vy1++;        vy3++; vy4++;
   vz0++; vz1++; vz2++; vz3++;
   }

if(bndflag == YZERO)
   {
   ip = 1;
/*
   1st, set vel components to zero at 'h' above FS, this is iz=0 for Vx and Vy.
   For Vz, we use the bottom row since it is not used for abs. boundary and would
   correspond to iz=-1.
*/

   vx0 = pvptr[ip];
   zero(vx0,nx);

   vy0 = pvptr[ip] +   nx*nz;
   zero(vy0,nx);

   vz0 = pvptr[ip] + 2*nx*nz + (nz-1)*nx;
   zero(vz0,nx);

/*
   Next solve for Vz at h/2 above FS, all 2nd order.
*/

   l2mp = medp[ip] +            nx + 1;
   lp   = medp[ip] +    nx*nz + nx + 1;

   vx1 = pvptr[ip]   +           nx;
   vx2 = pvptr[ip]   +           nx + 1;

   vy1 = pvptr[ip-1] +   nx*nz + nx + 1;
   vy2 = pvptr[ip]   +   nx*nz + nx + 1;

/* vz0 is assumed to be =0.0, we solve for vz1 */
   vz1 = pvptr[ip]   + 2*nx*nz +        1;
   vz2 = pvptr[ip]   + 2*nx*nz +   nx + 1;

   n1 = nx - 1;
   while(n1--)
      {
      ff = sum*(a0 + a1*lp[0]/l2mp[0]);

      vz1[0] = vz2[0] + ff*(vx2[0] - vx1[0] + vy2[0] - vy1[0]);

      vx1++; vx2++;
      vy1++; vy2++;
      vz1++; vz2++;
      lp++; l2mp++;
      }
   }
}
