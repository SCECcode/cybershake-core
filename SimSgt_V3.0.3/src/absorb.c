/*
   absorb.c contains the following functions:

      absorb()
      abs_xzbnd()
      abs_xzbnd1()
      abs_ybnd()
      storev()
      store_iring()
      store()
*/

#include "include.h"

/*
	storev() copies the previous time step of vx, vy, and vz into
	the array pbnd.  These fields will subsequently be used in abs_ybnd()
	to absorb the wave field at the boundaries iy=0 and iy=ny-1.
*/

storev(nx,nz,pvptr,pbnd,iflag)
float **pvptr, *pbnd;
int nx, nz, iflag;
{
float *vx, *vy, *vz;
float *vxi, *vyi, *vzi;
int ip, yshft;

/*
   if iflag = YZERO, retain previous time step on plane 1 when iy=0
   if iflag = YN, when iy=ny-1, retain previous time step on plane 1 for vy
		  and on plane 2 for vx and vz.
*/
if(iflag == YZERO)
   {
   yshft = 0;
   ip = 1;
   }
else if(iflag == YN)
   {
   yshft = 1;
   ip = 2;
   }

vx  = pvptr[ip];
vy  = pvptr[ip-yshft] +   nx*nz;
vz  = pvptr[ip]       + 2*nx*nz;
vxi = pbnd;
vyi = pbnd +   nx*nz;
vzi = pbnd + 2*nx*nz;

store(nx*nz,vx,1,vxi,1);
store(nx*nz,vy,1,vyi,1);
store(nx*nz,vz,1,vzi,1);
}

/*
	When iflag=YZERO (iy=0 boundary), abs_ybnd() updates the velocity
        variables vx, vy and vz for the new time step on plane 0 using an
        absorbing boundary condition.
      
	When iflag=YN (iy=ny-1 boundary), abs_ybnd() updates the velocity
        variable vy for the new time step on plane 2 and the velocity
	variables vx and vz for the new time step on plane 3 using an
        absorbing boundary condition.
      
	The values of vx, vy and vz for the previous time step on plane 1
        (iflag=YZERO) or plane 1 or 2 (iflag=YN) are stored in the array
	pbnd in the following manner:

           pbnd         -> vx
           pbnd+nx*nz   -> vy
           pbnd+2*nx*nz -> vz
*/

abs_ybnd(nx,nz,pvptr,medptr,pbnd,fdc,fs,iflag,intfac)
struct interface *intfac;
struct fdcoefs *fdc;
float **pvptr, **medptr, *pbnd;
int nx, nz, fs, iflag;
{
float *vx1, *vy1, *vz1;
float *vx0, *vy0, *vz0;
float *vxi, *vyi, *vzi;
float *vxir, *vyir, *vzir;
float *med;
int iz, ip, iv0, iv1, ipabs, isabs, iabs, kf;
int off1, off2, off3, yshft;

/*
         iflag=YZERO -> update vx, vy and vz on plane 0 using current and
                        previous time step values on plane 1.
         iflag=YN    -> update vy on plane 2 using current and
                        previous time step values on plane 1 and
                        update vx and vz on plane 3 using current and
                        previous time step values on plane 2.

	 v1 (index pointer iv1) is the velocity plane being updated and
	 v0 (index pointer iv1) is the plane one grid level into the model,
	 eg., for iflag=YZERO, v1 is plane 0 and v0 is plane 1.
*/

if(iflag == YZERO)
   {
   yshft = 0;
   ip = 1;
   iv0 = 1;
   iv1 = 0;
   }
if(iflag == YN)
   {
   yshft = 1;
   ip = 2;
   iv0 = 2;
   iv1 = 3;
   }

 /* set pointers for vx, vy and vz update on interior of plane */
 
vx1 = pvptr[iv1];             /* update target planes */
vy1 = pvptr[iv1-yshft] +   nx*nz;
vz1 = pvptr[iv1]       + 2*nx*nz;
vx0 = pvptr[iv0];             /* current time step interior planes */
vy0 = pvptr[iv0-yshft] +   nx*nz;
vz0 = pvptr[iv0]       + 2*nx*nz;
vxi = pbnd;                   /* previous time step interior planes */
vyi = pbnd +   nx*nz;
vzi = pbnd + 2*nx*nz;

/*
   If interfacing,

   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx1,vy1,vz1)
      -along the inner edge of boundary for previous time step (vxi,vyi,vzi)
      -along the inner edge of boundary for current time step (vx0,vy0,vz0)

   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx1,vy1,vz1)
      -along the inner edge of boundary for current time step (vx0,vy0,vz0)
*/

/*
   do x-strips for interior points; stride=1, don't do corner point
*/
 
if(intfac->go)
   {
   if(iflag == YZERO)
      intfac->yshft = 3;
   if(iflag == YN)
      intfac->yshft = 0;

   med = medptr[ip];

   kf = 1; /* do all three components */
   sub_infc_bnd(nx,nz,vx1,vx0,vxi,vy1,vy0,vyi,vz1,vz0,vzi,kf,med,intfac,OPEN_READ);
   }
 
/* retain previous time step on inner ring for x and z boundary edges */

vxir = pbnd + 3*nx*nz;
vyir = pbnd + 3*nx*nz + 2*(nx+nz);
vzir = pbnd + 3*nx*nz + 4*(nx+nz);
store_iring(vxi,vxir,vyi,vyir,vzi,vzir,nx,nz);

iabs = VAVG;

for(iz=1;iz<nz-1;iz++)
   {
   ipabs = PABS;
   isabs = SABS;
   if(iz == 1 || iz == (nz-2))
      {
      ipabs = VAVG;
      isabs = VAVG;
      }

   off1 = iz*nx + 1;
   med = medptr[ip] + off1;

   absorb(nx-1,vx1+off1,vx0+off1,vxi+off1,med,nx*nz,&fdc->dtoh,1,isabs,NOCRN);
   if(iz < nz-2)
      absorb(nx,vz1+off1,vz0+off1,vzi+off1,med,nx*nz,&fdc->dtoh,1,isabs,NOCRN);
   absorb(nx,vy1+off1,vy0+off1,vyi+off1,med,nx*nz,&fdc->dtoh,1,ipabs,NOCRN);
   }

/*
   do x and z boundary edges; use vxir, vyir, vzir for inner ring values
*/

                           /* iz=0 boundary */
if(!fs)
   {
   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[ip] + off2;
   absorb(nx,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   absorb(nx-1,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   absorb(nx,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   }

                          /* ix=nx-1 boundary */

off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
med = medptr[ip] + off2;
absorb(nz,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);

off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
med = medptr[ip] + off2;
absorb(nz-1,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);

off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
med = medptr[ip] + off2;
absorb(nz,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);

                          /* iz=nz-1 boundary */

off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
med = medptr[ip] + off2;
absorb(nx,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);

off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
med = medptr[ip] + off2;
absorb(nx-1,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);

off1 = nz*nx - 2; off2 = (nz-1)*nx - 2; off3 = nx + nz + 1;
med = medptr[ip] + off2;
absorb(nx,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);

                           /* ix=0 boundary */

off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
med = medptr[ip] + off2;
absorb(nz,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);

off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
med = medptr[ip] + off2;
absorb(nz-1,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);

off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
med = medptr[ip] + off2;
absorb(nz,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);

if(intfac->go)
   {
   if(iflag == YZERO)
      intfac->yshft = 3;
   if(iflag == YN)
      intfac->yshft = 0;

   med = medptr[ip];

   kf = 1; /* do all three components */
   add_infc_bnd(nx,nz,vx1,vx0,vy1,vy0,vz1,vz0,kf,med,intfac,DUMY);
   }
}

abs_xzbnd1(nx,nz,pvptr,pbnd,dtoh,medptr,fs,intfac)
struct interface *intfac;
int nx, nz, fs;
float **pvptr, *pbnd, **medptr, *dtoh;
{
float *vx, *vz;
float *vxi, *vzi;
int off1, off2, off3, kf;
float *med;
int iabs;

vx  = pvptr[1];
vz  = pvptr[1] + 2*nx*nz;

vxi = pbnd + 3*nx*nz;
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

iabs = VAVG;

/*
   If interfacing,
   
   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx,vz)
      -along the inner edge of boundary for previous time step (vxi,vzi)
      -along the inner edge of boundary for current time step (vx,vz)
   
   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx,vz)
      -along the inner edge of boundary for current time step (vx,vz)
*/

if(intfac->go)
   {
   intfac->yshft = 2;
   med = medptr[1];

   kf = 2; /* only do vx and vz components (send vx as dummy pointer for vy) */
   sub_infc_bnd(nx,nz,vx,vx,vxi,vx,vx,vx,vz,vz,vzi,kf,med,intfac,OPEN_READ);
   }

/* 1st order absorbing for x and z boundaries */

                           /* iz=0 boundary */

if(!fs)
   {   
   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[1] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,1,iabs,CRN_2D);
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,1,iabs,CRN_2D);
   }
 
                          /* ix=nx-1 boundary */
 
off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
med = medptr[1] + off2;
absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,nx,iabs,CRN_2D);
 
off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
med = medptr[1] + off2;
absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,nx,iabs,CRN_2D);

                          /* iz=nz-1 boundary */
 
off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
med = medptr[1] + off2;
absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-1,iabs,CRN_2D);
 
off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
med = medptr[1] + off2;
absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-1,iabs,CRN_2D);

                           /* ix=0 boundary */
 
off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
med = medptr[1] + off2;
absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-nx,iabs,CRN_2D);
 
off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
med = medptr[1] + off2;
absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-nx,iabs,CRN_2D);

if(intfac->go)
   {
   med = medptr[1];

   kf = 2; /* only do vx and vz components (send vx as dummy pointer for vy) */
   add_infc_bnd(nx,nz,vx,vx,vx,vx,vz,vz,kf,med,intfac,DUMY);
   }
}

abs_xzbnd(nx,nz,pvptr,pbnd,dtoh,medptr,iflag,fs,intfac)
struct interface *intfac;
int nx, nz, iflag, fs;
float **pvptr, *pbnd, **medptr, *dtoh;
{
float *vx, *vy, *vz;
float *vxi, *vyi, *vzi;
int off1, off2, off3;
int ipabs, isabs, vyabs, kf;
float *med;

vx  = pvptr[2];
vy  = pvptr[1] +   nx*nz;   /* center on plane 1 */
vz  = pvptr[2] + 2*nx*nz;

vxi = pbnd + 3*nx*nz;
vyi = pbnd + 3*nx*nz + 2*(nx+nz);
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

/*
   If interfacing,
   
   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx,vy,vz)
      -along the inner edge of boundary for previous time step (vxi,vyi,vzi)
      -along the inner edge of boundary for current time step (vx,vy,vz)
   
   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx,vy,vz)
      -along the inner edge of boundary for current time step (vx,vy,vz)
*/

/* 1st order absorbing for x and z boundaries */

ipabs = PABS; /* absorb P waves */
isabs = SABS; /* absorb S waves */
vyabs = SABS; /* absorb S waves (special for vy) */

if(iflag == YZERO || iflag == YN) /* vy near corner edge -> avg. velocities */
   vyabs = VAVG;
if(iflag == YN) /* vx and vz near corner edge -> average velocities */
   {
   ipabs = VAVG;
   isabs = VAVG;
   }

if(intfac->go)
   {
   intfac->yshft = 1;
   med = medptr[2];

   kf = 1; /* do all three components */
   sub_infc_bnd(nx,nz,vx,vx,vxi,vy,vy,vyi,vz,vz,vzi,kf,med,intfac,OPEN_READ);
   }
 
                           /* iz=0 boundary */
 
if(!fs)
   {

   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[2] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,1,ipabs,CRN_2D);
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,1,isabs,CRN_2D);
 
   med = medptr[1] + off2;
   absorb(nx,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,1,vyabs,CRN_2D);
   }
 
                          /* ix=nx-1 boundary */
 
off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
med = medptr[2] + off2;
absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,nx,ipabs,CRN_2D);
 
off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
med = medptr[2] + off2;
absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,nx,isabs,CRN_2D);
 
off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
med = medptr[1] + off2;
absorb(nz,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,nx,vyabs,CRN_2D);

                          /* iz=nz-1 boundary */
 
off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
med = medptr[2] + off2;
absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-1,ipabs,CRN_2D);
 
off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
med = medptr[2] + off2;
absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-1,isabs,CRN_2D);
 
off1 = nz*nx - 2; off2 = (nz-1)*nx - 2; off3 = nx + nz + 1;
med = medptr[1] + off2;
absorb(nx,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,-1,vyabs,CRN_2D);

                           /* ix=0 boundary */
 
off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
med = medptr[2] + off2;
absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-nx,ipabs,CRN_2D);
 
off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
med = medptr[2] + off2;
absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-nx,isabs,CRN_2D);
 
off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
med = medptr[1] + off2;
absorb(nz,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,-nx,vyabs,CRN_2D);

if(intfac->go)
   {
   intfac->yshft = 1;
   med = medptr[2];

   kf = 1; /* do all three components */
   add_infc_bnd(nx,nz,vx,vx,vy,vy,vz,vz,kf,med,intfac,DUMY);
   }
}

/*
      store_iring() stores the previous time step values of velocity on
      a ring one grid point interior to the grid boundary.  The ring is
      broken into four segments which are stored sequentially in the array
      pbnd.  These values are later used in absorb() to apply the absorbing
      boundary condition.  The order in which these values are stored has
      some redundancy, but it allows for a more simplistic storage format.
      The values are stored sequentially beginning in the upper left corner
      (see diagram below) and proceding in a clockwise fashion, as indicated
      by the '#' symbols.
 
                               iz=0  
      
                             second row starts here
                                        |
                                       \|/
                        o # o o o o o o # o         
first row starts here -># # # # # # # # # #         
                        o # o o o o o o # o         
                        o # o o o o o o # o         
                        o # o o o o o o # o   ix=nx-1
             ix=0       o # o o o o o o # o         
                        o # o o o o o o # o         
                        o # o o o o o o # o         
                        # # # # # # # # # #<- third row starts here
                        o # o o o o o o # o
                         /|\
                          |
                fourth row starts here

			      iz=nz-1

      Thus the solution procedes in the following manner: 1) iz=0 boundary,
      2) ix=nx-1 boundary, 3) iz=nz-1 boundary, 4) ix=0 boundary.
      The first and last values in each row are never used in absorb(), but
      simply occupy spots in the array.  In addition, there is overlap between
      adjacent rows for each of the interior corner points.  Storing the values
      this ways allows the array to be filled with 2 rows of length nx (the
      first and third rows) and 2 rows of length nz (the second and fourth
      rows), which makes for easier bookkeeping.
*/

store_iring(vx,vxi,vy,vyi,vz,vzi,nx,nz)
float *vx, *vxi, *vy, *vyi, *vz, *vzi;
int nx, nz;
{
store(nx,vx+nx,1,vxi,1);
store(nx,vy+nx,1,vyi,1);
store(nx,vz+nx,1,vzi,1);

store(nz,vx+nx-3,nx,vxi+nx,1);
store(nz,vy+nx-2,nx,vyi+nx,1);
store(nz,vz+nx-2,nx,vzi+nx,1);

store(nx,vx+(nz-1)*nx-1,-1,vxi+nx+nz,1);
store(nx,vy+(nz-1)*nx-1,-1,vyi+nx+nz,1);
store(nx,vz+(nz-2)*nx-1,-1,vzi+nx+nz,1);

store(nz,vx+(nz-1)*nx+1,-nx,vxi+2*nx+nz,1);
store(nz,vy+(nz-1)*nx+1,-nx,vyi+2*nx+nz,1);
store(nz,vz+(nz-1)*nx+1,-nx,vzi+2*nx+nz,1);
}

store(n,q1,stride1,q2,stride2)
register float *q1, *q2;
int n, stride1, stride2;
{
while(n--)
   {
   q2[0] = q1[0];
   q1 = q1 + stride1;
   q2 = q2 + stride2;
   }
}

/*
      absorb() applies plane wave absorbing boundary condition to wave field.
      The variables are as follows:

	 vo     -> wave field at previous time step on outer edge of boundary,
	           this will be updated to current time step
	 vi     -> wave field at current time step on inner edge of boundary
	 vim    -> wave field at previous time step on inner edge of boundary

       The boundary points are updated in a sequential fashion which procedes
       around the boundary of the plane in a clockwise manner [see comments
       for function store_iring()].  The pointers vo, vi and vim must point
       to the second element in the array as illustrated below (ie., skip the
       point labeled by 'o').

                      pointers set here
                              |
                             \|/
                       vo ->o * * * * * * * * *
                  vi, vim ->o # # # # # # # # o

		      solution direction -->

       The solution is calculated for the points labeled by '*', using the
       points labeled by '#'.  At the end of the row, the interior corner
       point is used to obtain the update at the boundary corner point,
       using a 45 degree rotation of the boundary condition [h -> h/sqrt(2)].

       Note that the pointers vo and vi are set within the xz plane, while
       the pointer vim is set to the vector pbnd.  As the solution procedes
       aroung the edge of the plane, the pointer to vim is incremented
       sequentially by one (vim++).  However, since the xz planes are stored
       in strips of x values (x is fastest variable), the pointers vo and vi
       cannot always be incremented by one to get to the next point.  These
       pointers are reset using the variable 'stride'.  For each boundary
       segement (iz=0, ix=nx-1, iz=nz-1, ix=0), the appropriate
       values for the variable stirde are illustrated in the diagram below.
 
                           iz=0, stride=1
      
                        o # o o o o o o # o         
           start here -># # # # # # # # # #         
                        o # o o o o o o # o         
                        o # o o o o o o # o         
                        o # o o o o o o # o
      ix=0, stride=-nx  o # o o o o o o # o   ix=nx-1, stride=nx
                        o # o o o o o o # o         
                        o # o o o o o o # o         
                        # # # # # # # # # #
                        o # o o o o o o # o

		         iz=nz-1, stride=-1
        
	Thus, stride increments the pointers vo and vi within the xz plane
	in a clockwise fashion around the border of the model space.

	 wvtype -> wave type being propagated at boundary:

	 if wvtype=PABS,  then P waves are absorbed
	 if wvtype=SABS,  then S waves are absorbed
	 if wvtype=VAVG,  then P and S velocities are averaged (edge points
			  near iy=0,ny-1)

	 iflag -> flag for doing corner point

         if iflag=0,      then corner point is not done (interior of planes
                          at iy=0,ny-1)
         if iflag=CRN_2D, then corner point is done (corner points of
			    of xz planes where iy!=0,ny-1)
         if iflag=CRN_3D, then dtoh -> dtoh/sqrt(2) for edge points of
			    xz planes at iy=0,ny-1 and corner point of these
			    planes are done with dtoh -> dtoh/sqrt(3)

      The variable n sould be set equal to the total number of grid points
      along that strip of the model (eg., n=nx).  The first corner point is
      not done [eg., ix=0, it is the last corner point from the previous call
      to absorb()], and the last corner point (eg., ix=nx-1) along with the
      two points adjacent to the first and last corner points (eg., ix=1 and
      ix=nx-2) are done outside the while loop.  Thus, the loop is done over
      ninc = n-4 points.
*/

absorb2(n,vo,vi,vim,medp,off,dtoh,stride,wvtype,iflag)
register float *vo, *vi, *vim;
float *medp, *dtoh;
int n, stride, wvtype, off, iflag;
{
return;
}

absorb_SYM(n,vo,vi,vim,medp,off,dtoh,stride,wvtype,iflag)
register float *vo, *vi, *vim;
float *medp, *dtoh;
int n, stride, wvtype, off, iflag;
{
int ninc;

ninc = n-3;
while(ninc--)
   {
   vo[0] = vi[0];
 
   vo = vo + stride;
   vi = vi + stride;
   }
vo[0] = vi[0];
 
if(iflag == CRN_2D || iflag == CRN_3D)
   {
   vo = vo + stride;
 
   vo[0] = vi[0];
   }
}

absorb(n,vo,vi,vim,medp,off,dtoh,stride,wvtype,iflag)
register float *vo, *vi, *vim;
float *medp, *dtoh;
int n, stride, wvtype, off, iflag;
{
float *l2mp, *lp, *mp, *bp;
float cosA, c1, c2, dtohc, edgfac, f1;
float half = 0.5;
float one = 1.0;
float sqrt2 = 1.4142136;
float sqrt3 = 1.7320508;
float bcvtol = 4.0;
int ninc;

float vp, vs;
int i;

l2mp = medp;
lp   = medp +   off;
mp   = medp + 2*off; /* use mxyp */
bp   = medp + 5*off; /* use bxp */
 
cosA = one;
cosA = cos(3.14159*45.0/180.0);

edgfac = one;
if(iflag == CRN_3D)
   edgfac = one/sqrt2;

/* RWG test 06/14/00 */
/* don't apply any angle corrections at adjacent corner points */

cosA = one;
edgfac = one;

/* do beginning adjacent corner point */
/* average alpha & beta since partitioning isn't known */

c1 = edgfac*half*(*dtoh)*(sqrt(l2mp[0]*bp[0]) + sqrt(mp[0]*bp[0]));

c2 = (cosA - c1)/(cosA + c1);
vo[0] = vim[0] + c2*(vo[0] - vi[0]);
 
vo = vo + stride;
vi = vi + stride;
vim++;
 
l2mp  = l2mp  + stride;
mp = mp + stride;
bp = bp  + stride;
 
ninc = n-4; /* loop over interior points */

if(wvtype == PABS) /* do P waves */
   {
   while(ninc--)
      {
      c1 = (*dtoh)*sqrt(l2mp[0]*bp[0]);
 
      c2 = (one - c1)/(one + c1);
      vo[0] = vim[0] + c2*(vo[0] - vi[0]);
 
      vo = vo + stride;
      vi = vi + stride;
      vim++;
 
      l2mp = l2mp + stride;
      mp = mp + stride;
      bp = bp + stride;
      }
   }
else if(wvtype == SABS) /* do S waves */
   {
   while(ninc--)
      {
      c1 = (*dtoh)*sqrt(mp[0]*bp[0]);
      c2 = (*dtoh)*sqrt(l2mp[0]*bp[0]);

      if(c2/c1 > bcvtol)
	 {
      /*
	 fprintf(stderr,"bcvtol exceeded: vp=%f vs=%f\n",c2/(*dtoh),c1/(*dtoh));
 */
	 c1 = half*(c1+c2);
	 }
 
      c2 = (one - c1)/(one + c1);
      vo[0] = vim[0] + c2*(vo[0] - vi[0]);
 
      vo = vo + stride;
      vi = vi + stride;
      vim++;
 
      l2mp = l2mp + stride;
      mp = mp + stride;
      bp = bp + stride;
      }
   }
else if(wvtype == VAVG) /* average P and S velocities */
   {
   while(ninc--)
      {
      c1 = edgfac*half*(*dtoh)*(sqrt(l2mp[0]*bp[0]) + sqrt(mp[0]*bp[0]));
 
      c2 = (one - c1)/(one + c1);
      vo[0] = vim[0] + c2*(vo[0] - vi[0]);
 
      vo = vo + stride;
      vi = vi + stride;
      vim++;
 
      l2mp = l2mp + stride;
      mp = mp + stride;
      bp = bp + stride;
      }
   }
else
   {
   fprintf(stderr,"**** Invalid choice, wvtype=%d in absorb(), exiting...\n",wvtype);
   exit(-1);
   }
 
/* do ending adjacent corner point */
/* average alpha & beta since partitioning isn't known */
 
c1 = edgfac*half*(*dtoh)*(sqrt(l2mp[0]*bp[0]) + sqrt(mp[0]*bp[0]));

c2 = (cosA - c1)/(cosA + c1);
vo[0] = vim[0] + c2*(vo[0] - vi[0]);
 
/* do corner point if iflag=CRN_2D or iflag=CRN_3D */
 
if(iflag == CRN_2D)
   {
   vo = vo + stride;
 
   dtohc = (*dtoh)/sqrt2;
   c1 = half*dtohc*(sqrt(l2mp[0]*bp[0]) + sqrt(mp[0]*bp[0]));

   c2 = (one - c1)/(one + c1);
   vo[0] = vim[0] + c2*(vo[0] - vi[0]);
   }
else if(iflag == CRN_3D)
   {
   vo = vo + stride;
 
   dtohc = (*dtoh)/sqrt3;
   c1 = half*dtohc*(sqrt(l2mp[0]*bp[0]) + sqrt(mp[0]*bp[0]));

   c2 = (one - c1)/(one + c1);
   vo[0] = vim[0] + c2*(vo[0] - vi[0]);
   }
}

void abs_xzbndP3(int nx,int nz,float **pvptr,float *pbnd,float *dtoh,float **medptr,int ybndflag,int fs,struct interface *intfac,struct nodeinfo *ni)
{
float *vx, *vy, *vz;
float *vxi, *vyi, *vzi;
int off1, off2, off3;
int ipabs, isabs, vyabs, kf;
float *med;

vx  = pvptr[2];
vy  = pvptr[1] +   nx*nz;   /* center on plane 1 */
vz  = pvptr[2] + 2*nx*nz;

vxi = pbnd + 3*nx*nz;
vyi = pbnd + 3*nx*nz + 2*(nx+nz);
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

/*
   If interfacing,
   
   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx,vy,vz)
      -along the inner edge of boundary for previous time step (vxi,vyi,vzi)
      -along the inner edge of boundary for current time step (vx,vy,vz)
   
   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx,vy,vz)
      -along the inner edge of boundary for current time step (vx,vy,vz)
*/

/* 1st order absorbing for x and z boundaries */

ipabs = PABS; /* absorb P waves */
isabs = SABS; /* absorb S waves */
vyabs = SABS; /* absorb S waves (special for vy) */

if(ybndflag == YZERO || ybndflag == YN) /* vy near corner edge -> avg. velocities */
   vyabs = VAVG;
if(ybndflag == YN) /* vx and vz near corner edge -> average velocities */
   {
   ipabs = VAVG;
   isabs = VAVG;
   }

if(intfac->go)
   {
   intfac->yshft = 1;
   med = medptr[2];

   kf = 1; /* do all three components */
   sub_infc_bnd(nx,nz,vx,vx,vxi,vy,vy,vyi,vz,vz,vzi,kf,med,intfac,OPEN_READ);
   }
 
                           /* iz=0 boundary */
 
if(ni->minusId_z == -1 && !fs)
   {
   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[2] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,1,ipabs,CRN_2D);
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,1,isabs,CRN_2D);
 
   med = medptr[1] + off2;
   absorb(nx,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,1,vyabs,CRN_2D);
   }
 
                          /* ix=nx-1 boundary */
 
if(ni->plusId_x == -1)
   {
   off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
   med = medptr[2] + off2;
   absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,nx,ipabs,CRN_2D);
 
   off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
   med = medptr[2] + off2;
   absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,nx,isabs,CRN_2D);
 
   off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
   med = medptr[1] + off2;
   absorb(nz,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,nx,vyabs,CRN_2D);
   }

                          /* iz=nz-1 boundary */
 
if(ni->plusId_z == -1)
   {
   off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
   med = medptr[2] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-1,ipabs,CRN_2D);
 
   off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
   med = medptr[2] + off2;
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-1,isabs,CRN_2D);
 
   off1 = nz*nx - 2; off2 = (nz-1)*nx - 2; off3 = nx + nz + 1;
   med = medptr[1] + off2;
   absorb(nx,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,-1,vyabs,CRN_2D);
   }

                           /* ix=0 boundary */
 
if(ni->minusId_x == -1)
   {
   off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
   med = medptr[2] + off2;
   absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-nx,ipabs,CRN_2D);
 
   off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
   med = medptr[2] + off2;
   absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-nx,isabs,CRN_2D);
 
   off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
   med = medptr[1] + off2;
   absorb(nz,vy+off1,vy+off2,vyi+off3,med,nx*nz,dtoh,-nx,vyabs,CRN_2D);
   }

if(intfac->go)
   {
   intfac->yshft = 1;
   med = medptr[2];

   kf = 1; /* do all three components */
   add_infc_bnd(nx,nz,vx,vx,vy,vy,vz,vz,kf,med,intfac,DUMY);
   }
}

void abs_xzbnd1P3(int nx,int nz,float **pvptr,float *pbnd,float *dtoh,float **medptr,int fs,struct interface *intfac,struct nodeinfo *ni)
{
float *vx, *vz;
float *vxi, *vzi;
int off1, off2, off3, kf;
float *med;
int iabs;

vx  = pvptr[1];
vz  = pvptr[1] + 2*nx*nz;

vxi = pbnd + 3*nx*nz;
vzi = pbnd + 3*nx*nz + 4*(nx+nz);

iabs = VAVG;

/*
   If interfacing,
   
   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx,vz)
      -along the inner edge of boundary for previous time step (vxi,vzi)
      -along the inner edge of boundary for current time step (vx,vz)
   
   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx,vz)
      -along the inner edge of boundary for current time step (vx,vz)
*/

if(intfac->go)
   {
   intfac->yshft = 2;
   med = medptr[1];

   kf = 2; /* only do vx and vz components (send vx as dummy pointer for vy) */
   sub_infc_bnd(nx,nz,vx,vx,vxi,vx,vx,vx,vz,vz,vzi,kf,med,intfac,OPEN_READ);
   }

/* 1st order absorbing for x and z boundaries */

                           /* iz=0 boundary */

if(ni->minusId_z == -1 && !fs)
   {   
   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[1] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,1,iabs,CRN_2D);
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,1,iabs,CRN_2D);
   }
 
                          /* ix=nx-1 boundary */
 
if(ni->plusId_x == -1)
   {
   off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
   med = medptr[1] + off2;
   absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,nx,iabs,CRN_2D);
 
   off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
   med = medptr[1] + off2;
   absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,nx,iabs,CRN_2D);
   }

                          /* iz=nz-1 boundary */
 
if(ni->plusId_z == -1)
   {
   off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
   med = medptr[1] + off2;
   absorb(nx,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-1,iabs,CRN_2D);
 
   off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
   med = medptr[1] + off2;
   absorb(nx-1,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-1,iabs,CRN_2D);
   }

                           /* ix=0 boundary */
 
if(ni->minusId_x == -1)
   {
   off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
   med = medptr[1] + off2;
   absorb(nz,vx+off1,vx+off2,vxi+off3,med,nx*nz,dtoh,-nx,iabs,CRN_2D);
 
   off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
   med = medptr[1] + off2;
   absorb(nz-1,vz+off1,vz+off2,vzi+off3,med,nx*nz,dtoh,-nx,iabs,CRN_2D);
   }

if(intfac->go)
   {
   med = medptr[1];

   kf = 2; /* only do vx and vz components (send vx as dummy pointer for vy) */
   add_infc_bnd(nx,nz,vx,vx,vx,vx,vz,vz,kf,med,intfac,DUMY);
   }
}

/*
	When ybndflag=YZERO (iy=0 boundary), abs_ybndP3() updates the velocity
        variables vx, vy and vz for the new time step on plane 0 using an
        absorbing boundary condition.
      
	When ybndflag=YN (iy=ny-1 boundary), abs_ybndP3() updates the velocity
        variable vy for the new time step on plane 2 and the velocity
	variables vx and vz for the new time step on plane 3 using an
        absorbing boundary condition.
      
	The values of vx, vy and vz for the previous time step on plane 1
        (ybndflag=YZERO) or plane 1 or 2 (ybndflag=YN) are stored in the array
	pbnd in the following manner:

           pbnd         -> vx
           pbnd+nx*nz   -> vy
           pbnd+2*nx*nz -> vz
*/

void abs_ybndP3(int nx,int nz,float **pvptr,float **medptr,float *pbnd,struct fdcoefs *fdc,int fs,int ybndflag,struct interface *intfac,struct nodeinfo *ni)
{
float *vx1, *vy1, *vz1;
float *vx0, *vy0, *vz0;
float *vxi, *vyi, *vzi;
float *vxir, *vyir, *vzir;
float *med;
int iz, ip, iv0, iv1, ipabs, isabs, iabs, kf;
int off1, off2, off3, yshft;
int izstart, izend;

/*
         ybndflag=YZERO -> update vx, vy and vz on plane 0 using current and
                        previous time step values on plane 1.
         ybndflag=YN    -> update vy on plane 2 using current and
                        previous time step values on plane 1 and
                        update vx and vz on plane 3 using current and
                        previous time step values on plane 2.

	 v1 (index pointer iv1) is the velocity plane being updated and
	 v0 (index pointer iv1) is the plane one grid level into the model,
	 eg., for ybndflag=YZERO, v1 is plane 0 and v0 is plane 1.
*/

if(ybndflag == YZERO)
   {
   yshft = 0;
   ip = 1;
   iv0 = 1;
   iv1 = 0;
   }
if(ybndflag == YN)
   {
   yshft = 1;
   ip = 2;
   iv0 = 2;
   iv1 = 3;
   }

 /* set pointers for vx, vy and vz update on interior of plane */
 
vx1 = pvptr[iv1];             /* update target planes */
vy1 = pvptr[iv1-yshft] +   nx*nz;
vz1 = pvptr[iv1]       + 2*nx*nz;
vx0 = pvptr[iv0];             /* current time step interior planes */
vy0 = pvptr[iv0-yshft] +   nx*nz;
vz0 = pvptr[iv0]       + 2*nx*nz;
vxi = pbnd;                   /* previous time step interior planes */
vyi = pbnd +   nx*nz;
vzi = pbnd + 2*nx*nz;

/*
   If interfacing,

   1) Before call to absorb(), subtract out interface field from points
      -along the outer edge of boundary for previous time step (vx1,vy1,vz1)
      -along the inner edge of boundary for previous time step (vxi,vyi,vzi)
      -along the inner edge of boundary for current time step (vx0,vy0,vz0)

   2) After call to absorb(), add in interface field to points
      -along the outer edge of boundary for current time step (vx1,vy1,vz1)
      -along the inner edge of boundary for current time step (vx0,vy0,vz0)
*/

/*
   do x-strips for interior points; stride=1, don't do corner point
*/
 
if(intfac->go)
   {
   if(ybndflag == YZERO)
      intfac->yshft = 3;
   if(ybndflag == YN)
      intfac->yshft = 0;

   med = medptr[ip];

   kf = 1; /* do all three components */
   sub_infc_bnd(nx,nz,vx1,vx0,vxi,vy1,vy0,vyi,vz1,vz0,vzi,kf,med,intfac,OPEN_READ);
   }
 
/* retain previous time step on inner ring for x and z boundary edges */

vxir = pbnd + 3*nx*nz;
vyir = pbnd + 3*nx*nz + 2*(nx+nz);
vzir = pbnd + 3*nx*nz + 4*(nx+nz);
store_iring(vxi,vxir,vyi,vyir,vzi,vzir,nx,nz);

iabs = VAVG;

izstart = 2;
if(ni->minusId_z == -1)
   izstart = 1;

izend = nz-2;
if(ni->plusId_z == -1)
   izend = nz-1;

for(iz=izstart;iz<izend;iz++)
   {
   ipabs = PABS;
   isabs = SABS;
   if(iz == 1 || iz == (nz-2))
      {
      ipabs = VAVG;
      isabs = VAVG;
      }

   off1 = iz*nx + 1;
   med = medptr[ip] + off1;

   absorb(nx-1,vx1+off1,vx0+off1,vxi+off1,med,nx*nz,&fdc->dtoh,1,isabs,NOCRN);
   if(iz < nz-2)
      absorb(nx,vz1+off1,vz0+off1,vzi+off1,med,nx*nz,&fdc->dtoh,1,isabs,NOCRN);
   absorb(nx,vy1+off1,vy0+off1,vyi+off1,med,nx*nz,&fdc->dtoh,1,ipabs,NOCRN);
   }

/*
   do x and z boundary edges; use vxir, vyir, vzir for inner ring values
*/

                           /* iz=0 boundary */

if(ni->minusId_z == -1 && !fs)
   {
   off1 = 1; off2 = nx + 1; off3 = 1;
   med = medptr[ip] + off2;
   absorb(nx,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   absorb(nx-1,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   absorb(nx,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,1,iabs,CRN_3D);
   }

                          /* ix=nx-1 boundary */

if(ni->plusId_x == -1)
   {
   off1 = 2*nx - 2; off2 = 2*nx - 3; off3 = nx + 1;
   med = medptr[ip] + off2;
   absorb(nz,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);

   off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
   med = medptr[ip] + off2;
   absorb(nz-1,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);

   off1 = 2*nx - 1; off2 = 2*nx - 2; off3 = nx + 1;
   med = medptr[ip] + off2;
   absorb(nz,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,nx,iabs,CRN_3D);
   }

                          /* iz=nz-1 boundary */

if(ni->plusId_z == -1)
   {
   off1 = (nz-1)*nx - 2; off2 = (nz-2)*nx - 2; off3 = nx + nz + 1;
   med = medptr[ip] + off2;
   absorb(nx,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);

   off1 = nz*nx - 3; off2 = (nz-1)*nx - 3; off3 = nx + nz + 2;
   med = medptr[ip] + off2;
   absorb(nx-1,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);

   off1 = nz*nx - 2; off2 = (nz-1)*nx - 2; off3 = nx + nz + 1;
   med = medptr[ip] + off2;
   absorb(nx,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,-1,iabs,CRN_3D);
   }

                           /* ix=0 boundary */

if(ni->minusId_x == -1)
   {
   off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
   med = medptr[ip] + off2;
   absorb(nz,vx1+off1,vx0+off2,vxir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);

   off1 = (nz-3)*nx; off2 = (nz-3)*nx + 1; off3 = 2*nx + nz + 2;
   med = medptr[ip] + off2;
   absorb(nz-1,vz1+off1,vz0+off2,vzir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);

   off1 = (nz-2)*nx; off2 = (nz-2)*nx + 1; off3 = 2*nx + nz + 1;
   med = medptr[ip] + off2;
   absorb(nz,vy1+off1,vy0+off2,vyir+off3,med,nx*nz,&fdc->dtoh,-nx,iabs,CRN_3D);
   }

if(intfac->go)
   {
   if(ybndflag == YZERO)
      intfac->yshft = 3;
   if(ybndflag == YN)
      intfac->yshft = 0;

   med = medptr[ip];

   kf = 1; /* do all three components */
   add_infc_bnd(nx,nz,vx1,vx0,vy1,vy0,vz1,vz0,kf,med,intfac,DUMY);
   }
}
