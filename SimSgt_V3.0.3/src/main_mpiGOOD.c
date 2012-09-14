/*************************************************************************/
/*                                                                       */
/*         3-D Elastic Finite-difference Modeling Code                   */
/*                                                                       */
/*                    copyright (c) 1993                                 */
/*                     Robert W. Graves                                  */
/*                 Woodward-Clyde Consultants                            */
/*                   566 El Dorado Street                                */
/*                    Pasadena, CA 91101                                 */
/*                                                                       */
/*                    tel (818) 449-7650                                 */
/*                    fax (818) 449-3536                                 */
/*                                                                       */
/*    Permission to copy all or part of this work is granted,            */
/*    provided that the copies are not made or distributed               */
/*    for resale, and that the copyright notice and this                 */
/*    notice are retained.                                               */
/*                                                                       */
/*    NOTE:  This work is provided on an "as is" basis.                  */
/*           The author provides no warranty whatsoever,                 */
/*           expressed or implied, regarding the work,                   */
/*           including warranties with respect to its                    */
/*           fitness for any particular purpose.  The code               */
/*           is provided on the condition that it will not               */
/*           be transfered to third parties without the                  */
/*           explicit consent of the author.                             */
/*           In addition, the code has not been                          */
/*           tested extensively, therefore there are no                  */
/*           guarantees about its performance.  I am currently           */
/*           in the process of testing and benchmarking the              */
/*           code and would appreciate any comments and/or               */
/*           feedback from those who use it.                             */
/*                                                                       */
/*                             -RWG 05/05/93                             */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/*			version 3.0: xx/xx/11                            */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/*    A note about units.                                                */
/*                                                                       */
/*      - Global input units are:                                        */
/*                                                                       */
/*        # wavefield velocity is-                                       */
/*           vx,vy,vz       => cm/s                                      */
/*                                                                       */
/*        # seismic velocity and density are-                            */
/*           Vp,Vs     => km/s                                           */
/*           density   => gm/cm^3                                        */
/*                                                                       */
/*        # Lame parameters and stress are-                              */
/*           Tnn,Tij,lambda,mu => km*km*gm/(s*s*cm*cm*cm)                */
/*                             = 10 Nt/(m^2)                             */
/*                             = 10 Pa                                   */
/*                                                                       */
/*        # dimension and time measures-                                 */
/*           h (distance, length) => km                                  */
/*           dt (time)            => s                                   */
/*                                                                       */
/*        # sources are CGS-                                             */
/*           moment    => dyne-cm                                        */
/*           force     => dyne                                           */
/*           slip      => cm                                             */
/*           slip-rate => cm/s                                           */
/*                                                                       */
/*    If the input variables are expressed in these units, then all      */
/*    internal calculations can be done without any additional unit      */
/*    conversions.                                                       */
/*                                                                       */
/*       - For output, units are:                                        */
/*                                                                       */
/*        # wavefield velocities (single station & time slices) are-     */
/*           vx,vy,vz       => cm/s                                      */
/*                                                                       */
/*        # strain Greens tensor (SGT, reciprocal GF) are output in      */
/*          units of strain per unit force-                              */
/*           gxx,gyy,gzz,gxy,gxz,gyz => 1/(dyne)                         */
/*                                                                       */
/*          thus when convolved with a moment-rate function (dyne*cm/s)  */
/*          the result is velocity in (cm/s)-                            */
/*           dM/dt*gij => (dyne*cm/s) * 1/(dyne) => cm/s                 */
/*                                                                       */
/*************************************************************************/
/*****************    GENERAL INFORMATION    *****************************/
/*                                                                       */
/*    Axes orientation is as follows:                                    */
/*                                                                       */
/*         x-axis -> horizontal plane, positive to the east              */
/*         y-axis -> horizontal plane, positive to the south             */
/*         z-axis -> vertical plane, positive downward                   */
/*                                                                       */
/*    In order to run this code on a machine which has                   */
/*    limited core memory availble (ie., less than is needed             */
/*    to store the entire wave field and media parameters),              */
/*    the program uses disk memory to store the majority of              */
/*    the model as the wave field is being updated.                      */
/*    Field and media parameters are read in from disk in                */
/*    slices of xz planes.  This optimizes the number of                 */
/*    time steps that are calculated per IO operation.                   */
/*                                                                       */
/*    The variables which control these operations are:                  */
/*         span    -> number of time steps per IO operation              */
/*         spanlen =  4*span + 2                                         */
/*                 -> number of 2D slices of model needed in             */
/*                    in core at one time.                               */
/*         memlen  =  spanlen*nx*nz*size_float                           */
/*                 -> memory needed for "spanlen" model slices           */
/*                    for one variable (18 variables needed              */
/*                    for elastic case with Q).                          */
/*                                                                       */
/*    Field varibles are stored in pvf array as follows:                 */
/*                                                                       */
/*         pvf          -> vx                                            */
/*         pvf+   nx*nz -> vy                                            */
/*         pvf+ 2*nx*nz -> vz                                            */
/*         pvf+ 3*nx*nz -> txx                                           */
/*         pvf+ 4*nx*nz -> tyy                                           */
/*         pvf+ 5*nx*nz -> tzz                                           */
/*         pvf+ 6*nx*nz -> txy                                           */
/*         pvf+ 7*nx*nz -> txz                                           */
/*         pvf+ 8*nx*nz -> tyz                                           */
/*         pvf+ 9*nx*nz -> cxx                                           */
/*         pvf+10*nx*nz -> cyy                                           */
/*         pvf+11*nx*nz -> czz                                           */
/*         pvf+12*nx*nz -> cxy                                           */
/*         pvf+13*nx*nz -> cxz                                           */
/*         pvf+14*nx*nz -> cyz                                           */
/*                                                                       */
/*    The variables cij above are the memory variables for anelasticity. */
/*                                                                       */
/*    Media varibles are stored in medf array as follows:                */
/*                                                                       */
/*         medf          -> lambda + 2*mu                                */
/*         medf+   nx*nz -> lambda                                       */
/*         medf+ 2*nx*nz -> muxy                                         */
/*         medf+ 3*nx*nz -> muxz                                         */
/*         medf+ 4*nx*nz -> muyz                                         */
/*         medf+ 5*nx*nz -> bouyx                                        */
/*         medf+ 6*nx*nz -> bouyy                                        */
/*         medf+ 7*nx*nz -> bouyz                                        */
/*         medf+ 8*nx*nz -> ak                                           */
/*         medf+ 9*nx*nz -> pk                                           */
/*         medf+10*nx*nz -> sk                                           */
/*         medf+11*nx*nz -> lambda (raw, never altered from input)       */
/*         medf+12*nx*nz -> mu     (raw, never altered from input)       */
/*                                                                       */
/*    The variables ak, pk, and sk above are required to implement the   */
/*    memory variable formulation for anelasticity.                      */
/*                                                                       */
/*    The values for rigidity (mu) and bouyancy (1/density)              */
/*    have been averaged to avoid instability when media                 */
/*    contrasts are large (see Randall et al., Nov., 1991,               */
/*    Geophysics).  The averages are precomputed and stored              */
/*    in the array medf in order to make the calculations more           */
/*    efficient.  Refering to the table above, the media                 */
/*    parameters are averaged in the following manner:                   */
/*                                                                       */
/*       muxy -> avg. over xy plane for txy computation                  */
/*       muxz -> avg. over xz plane for txz computation                  */
/*       muyz -> avg. over yz plane for tyz computation                  */
/*       bouyx -> avg. over x plane for vx computation                   */
/*       bouyy -> avg. over y plane for vy computation                   */
/*       bouyz -> avg. over z plane for vz computation                   */
/*                                                                       */
/*    An arithmetic average is used for bouyancy and a                   */
/*    harmonic average is used for rigidity.  All averages               */
/*    are computed in the forward sense.  That is                        */
/*                                                                       */
/*       bouyx[i] = 2/(rho[i] + rho[i+1])                                */
/*       .                                                               */
/*       .                                                               */
/*       .                                                               */
/*       muxy[i] =                                                       */
/*          4/(1/mu1[i] + 1/mu1[i+1] + 1/mu2[i] + 1/mu2[i+1])            */
/*       .                                                               */
/*       .                                                               */
/*       .                                                               */
/*                                                                       */
/*    No averaging is performed at the maximal edges (ix=nx,             */
/*    iy=ny,iz=nz) of the model grid.                                    */
/*                                                                       */
/*    Absorbing boundary conditions are applied to velocities            */
/*    only, the outer ring of stress values are not updated.             */
/*    To implement the absorbing boundaries, the velocity                */
/*    field at the previous time step is required.  This is              */
/*    stored in the array pbnd.                                          */
/*                                                                       */
/*    For the boundaries at iy=0 and iy=ny-1, the planes of              */
/*    vx, vy and vz at iy=1 and iy=ny-2, respectively, are               */
/*    stored prior to updating the velocity field using the              */
/*    function storev().  The absorbing boundary conditions              */
/*    are applied using the function abs_ybnd() which is                 */
/*    called from the function tstepvbnd().  For these planes,           */
/*    the velocity fields are stored in the array pbnd in the            */
/*    following manner:                                                  */
/*                                                                       */
/*         pbnd         -> vx                                            */
/*         pbnd+nx*nz   -> vy                                            */
/*         pbnd+2*nx*nz -> vz                                            */
/*                                                                       */
/*    For the boundaries of x and z, the boundary condition              */
/*    for vx, vy and vz are applied using the function                   */
/*    absorb() which is called from tstepv().  The inner ring            */
/*    of values are stored using the function store_iring()              */
/*    in the following manner:                                           */
/*                                                                       */
/*         pbnd+3*nx*nz           -> vx                                  */
/*         pbnd+3*nx*nz+2*(nx+nz) -> vy                                  */
/*         pbnd+3*nx*nz+4*(nx+nz) -> vz                                  */
/*                                                                       */
/*    Velocity field is updated first, then stress field.                */
/*                                                                       */
/*************************************************************************/
/************************    NOTES    ************************************/
/*                                                                       */
/*  NOTE: (03/16/93 -RWG)                                                */
/*    The free surface boundary condition has not been tested            */
/*    completely.                                                        */
/*                                                                       */
/*  NOTE: (03/08/97 -RWG)                                                */
/*    The free surface boundary condition has been tested in many ways,  */
/*    and seems to perform very well.                                    */
/*                                                                       */
/*************************************************************************/
/*******************    MAJOR REVISIONS    *******************************/
/*                                                                       */
/*  02/12/93 -RWG                                                        */
/*    Changed calls to ord2(), ord4(), etc. so that stride=1             */
/*    in all cases.  I discovered that when nx and/or nz are             */
/*    large and stride=nx (ie., x derivatives in previous                */
/*    coding), the computation time can be slowed by up to a             */
/*    factor of 2.  I think this is caused by memory swapping            */
/*    when the arrays get large.  This may be a hardware                 */
/*    limitation or it may be set by the system software.  I             */
/*    have only noticed the problem on the Sparc2, it does not           */
/*    seem to occur on the Sun4.                                         */
/*                                                                       */
/*  04/06/93 -RWG                                                        */
/*    Boundaries along ix=nx-1, iy=ny-1 and iz=nz-1 are moved            */
/*    inward 1/2 grid step to stabilize the absorbing                    */
/*    boundaries.  Thus, the interior solution is not                    */
/*    calculated for vx along ix=nx-2, for vy along iy=ny-2              */
/*    or for vz along iz=nz-2.  These values are obtained                */
/*    using the absorbing boundary condition (first-order                */
/*    Clayton-Enquist).  The other variables affected are txx,           */
/*    tyy and tzz.                                                       */
/*                                                                       */
/*  07/29/94 -RWG                                                        */
/*    Set up memory configuration, so that model can either be           */
/*    stored internally in core memory or externally on disk,            */
/*    with user controlled IO.  Information for taking model             */
/*    slices in or out of memory is contained in the following           */
/*    structure:                                                         */
/*                                                                       */
/*            struct modelstorage                                        */
/*               {                                                       */
/*               size_t maxmem;                                          */
/*               int intmem;                                             */
/*               size_t ln_pvslice;                                      */
/*               size_t sz_pvslice;                                      */
/*               int pv_fdr;                                             */
/*               int pv_fdw;                                             */
/*               float *pv_buf;                                          */
/*               float *pv_bufptr;                                       */
/*               size_t ln_medslice;                                     */
/*               size_t sz_medslice;                                     */
/*               int med_fdr;                                            */
/*               int med_fdw;                                            */
/*               float *med_buf;                                         */
/*               float *med_bufptr;                                      */
/*               };                                                      */
/*                                                                       */
/*    The variables 'maxmem' (maximum core memory limit, in              */
/*    bytes) and 'intmem' (internal memory/external memory               */
/*    flag) can be input as as getpar parameters.  Default               */
/*    values are:                                                        */
/*                                                                       */
/*            maxmem = MEM_LIMIT (defined in defs.h)                     */
/*            intmem = 0 (use external storage)                          */
/*                                                                       */
/*    By setting 'intmem=1', the model will be stored entirely           */
/*    in core memory.  This will significantly reduce system             */
/*    CPU time, but does require a large amount of available             */
/*    memory.  The memory limit imposed by 'maxmem' is a soft            */
/*    limit, active only within this code.  I typically set              */
/*    the value of 'maxmem' to be one-half of the total                  */
/*    internal memory (excluding swap space) available on the            */
/*    machine that the code is running.                                  */
/*                                                                       */
/*    BE CAREFUL!!  When the code is run with the parameter              */
/*    'intmem=1', and other processes are accessing memory               */
/*    resources, then there exists a much greater potential              */
/*    for memory swapping to occur.  When this happens, the              */
/*    run-time will increase dramatically.  Unless you can be            */
/*    sure that this will not occur, then it is better to run            */
/*    the code with 'intmem=0'.  Furthermore, system overhead            */
/*    for user controlled IO when 'intmem=0' is usually on the           */
/*    order of 1-2 percent of the total CPU time.  Thus there            */
/*    appears to be very little advantage in using only core             */
/*    memory.                                                            */
/*                                                                       */
/*  03/08/97 -RWG                                                        */
/*                                                                       */
/*   -Changed moment-tensor source insertion to stresses (ala            */
/*    S. Day) instead of velocities (body forces as described            */
/*    in my paper).  The results are virtually identical, but            */
/*    the stress formulation is somewhat more compact.                   */
/*                                                                       */
/*   -Also added the option to specify input weights in the 'faultfile'  */
/*    in terms of relative slip, instead of relative moment.  This lets  */
/*    the scaling to moment to be done using the actual velocity model   */
/*    that the FD code is running.  This is important if the fault       */
/*    intersects a basin or encounters some lateral velocity variation.  */
/*    To envoke this option, set the parameter "relative_slip" in the    */
/*    getpar parameter file (eg., 'e3d.par') to a non-zero value, ie.,   */
/*                                                                       */
/*       relative_slip=1                                                 */
/*                                                                       */
/*    The default is still relative moments (relative_slip=0).           */
/*    Remember, that if you use this option, the weights in the          */
/*    'faultfile' must also be given in terms of relative slip.  There   */
/*    is a corresponding parameter ("relative_slip") that can be passed  */
/*    to the codes 'lisa2mtgrid' and 'lisa2mtgrid_vrup' so that these    */
/*    codes will output the FD source weights in terms of relative slip. */
/*    See me for more details.                                           */
/*                                                                       */
/*  05/09/98 -RWG                                                        */
/*                                                                       */
/*   -To alleviate instability at absorbing boundaries when media        */
/*    has very strong lateral contrast at close proximity (eg., 1 grid   */
/*    point) exactly parallel to boundary, the input velocity model is   */
/*    modified by replicating the plane of grid points adjacent to the   */
/*    boundary.  Thus, the iy=2 plane is copied to iy=1 and iy=0.  A     */
/*    similar type of replication is done for x and z (excluding free    */
/*    surface).  The replication (or padding) width is controlled by the */
/*    parameter NBND_PAD, and associated functions are copy_medslice()   */
/*    and xzpad_medslice().                                              */
/*                                                                       */
/*  07/98 -RWG                                                           */
/*                                                                       */
/*   -Version 1.6 has the added option to output moment-tensor           */
/*    components for use in calculating reciprocal GF's.  See mtheader   */
/*    info.                                                              */
/*                                                                       */
/*  02/99 -RWG                                                           */
/*                                                                       */
/*   -Version 1.7 has modified the output routines for MT components to  */
/*    include information about the model origin (lat,lon), grid         */
/*    spacing, and rigidity at the output location.  See mtheader        */
/*    structure information.                                             */
/*                                                                       */
/*  04/00 -RWG                                                           */
/*                                                                       */
/*   -Version 1.8 has modified IO routines for external model storage    */
/*    routines for large file compatibility.  Changes to the source code */
/*    are trivial (see open() calls in iofunc.c.  Compiler options are   */
/*    more noticable (see LF_FLAGS in makefile) and are definitely       */
/*    system and OS version dependent.  For SUNs the large file option   */
/*    is only available on solaris 2.6 and later OS versions.            */
/*                                                                       */
/*   -Version 1.8 has modified cpu time usage routines based on the      */
/*    system function times().  This was done to circumvent portability  */
/*    problems (still needs to be tested) and problems using procfs.h    */
/*    with large file (64-bit) configuration.                            */
/*                                                                       */
/*  06/00 -RWG                                                           */
/*                                                                       */
/*   -Version 1.9 has option to center media variables either at the     */
/*    normal stress node ("media_boundary=0") or at the NULL unit cell   */
/*    node ("media_boundary=1").  The NULL unit cell node is the only    */
/*    corner of the unit cell which does not have any field variables    */
/*    defined, and it lies diagonally opposite of the normal stress      */
/*    node.  Based on the definition of the exact location of the media  */
/*    variables, different formulas are needed to calculate the          */
/*    effective media parameters that are used for the time updates of   */
/*    the various wave field variables (see Zahradnik et al, 1993).      */
/*    The effective media parameters are calculated in the function      */
/*    avgmedslice() in genmodel.c.                                       */
/*                                                                       */
/*    * media_boundary=0                                                 */
/*    With the media variables centered at the normal stress node, the   */
/*    media boundaries go through the velocity nodes.  This results in a */
/*    0.5*h depth difference in the location of the media boundaries and */
/*    the location of the free surface (see Graves, 1996).  Thus, the    */
/*    thinnest possible layer has depth of '0.5*h', and all other layers */
/*    have a depth of '(k+0.5)*h' [k=integer].                           */
/*                                                                       */
/*    * media_boundary=1                                                 */
/*    With the media variables centered at the NULL unit cell node, the  */
/*    media boundaries cross at the normal stress node.  This results in */
/*    a co-location of the media boundaries and the location of the free */
/*    surface.  In this case, the thinnest possible layer has depth of   */
/*    'h', and all other layers have a depth of 'k*h' [k=integer].       */
/*                                                                       */
/*   -Version 1.9 has cleaned-up the coding used to generate the         */
/*    effective media parameters and write these values to the temporary */
/*    model files.                                                       */
/*                                                                       */
/*   -Version 1.9 has changed the input format specification for a       */
/*    finite fault in the file 'faultfile'.  The new format specifies    */
/*    the point source locations in terms of distances from the model    */
/*    origin.  The distances are floating point coordinates (x,y,z).     */
/*    Earlier versions specified the point source locations as grid      */
/*    indices (ix,iy,iz) relative to the model origin.  For the          */
/*    horizontal directions these are simply related by                  */
/*                                                                       */
/*                x = ix*h    y=iy*h                                     */
/*                                                                       */
/*    However, for the vertical direction the 1 grid downward shift for  */
/*    the free-surface implementation is NOT included in the new         */
/*    floating point faultfile coordinates.  This means that a point     */
/*    located exactly at the free surface in versions 1.8 and earlier    */
/*    is given by iz=1, whereas in version 1.9, this point is located    */
/*    at z=0.0.  The relation between these two specification is given   */
/*                                                                       */
/*                iz = (int)(z/h + 0.5) + 1                              */
/*                                                                       */
/*    The additional "1" in the above equation provides the downward     */
/*    shift to account for the free-surface.                             */
/*                                                                       */
/*    Note that the above format specification only applies to the       */
/*    coordinates specified in the 'faultfile'.  All other coordinate    */
/*    specifications in version 1.9 are given in terms of grid indices   */
/*    (ix,iy,iz), and these must include the 1 grid shift for the        */
/*    free-surface (just the same as earlier versions).  The reason for  */
/*    the format change in the 'faultfile' is to make the 'faultfile'    */
/*    specification compatible between uniform gird versions (v1.9) and  */
/*    non-uniform-grid version (v2.0) of emod3d.                         */
/*                                                                       */
/*  10/01 -RWG                                                           */
/*                                                                       */
/*   -Version 1.10 has implemented the coarse-grained memory variable    */
/*    formulation to model anelaticity.  See Day and Bradley (BSSA,      */
/*    June 2001) and my notes.                                           */
/*                                                                       */
/*  05/05 -RWG                                                           */
/*                                                                       */
/*   -Version 1.11 has implemented the RGF output option with the MPI    */
/*    parallel code formulation.                                         */
/*                                                                       */
/*  09/05 -RWG                                                           */
/*                                                                       */
/*   -Version 1.12:                                                      */
/*                                                                       */
/*     # Cleaned-up indexing for media planes in init_media() and        */
/*       eff_media().  Originally Sidao added two planes to media array  */
/*       (one in begining and one at the end) because he thought they    */
/*       were needed to perform media averaging.  Thus, media array was  */
/*       dimensioned at ny+2 instead of ny.  The indexing was adjusted   */
/*       in eff_media() in Versions 1.11 and earlier.                    */
/*                                                                       */
/*       Now, media array is dimensioned at ny and all indexing is set   */
/*       locally to iy=0,...,ny-1 or gloabally to iy=ny1,...,ny2-1.      */
/*       Functions affected by this are init_media(), eff_media(), and   */
/*       slip2mom().                                                     */
/*                                                                       */
/*     # reformatted 'struct momenttensor'                               */
/*                                                                       */
/*       * 'struct pntsrcs' (point sources) now has a pointer to         */
/*         'struct momenttensor'.  Memory (for nsource) is allocated in  */
/*         init_dc(). Should be backward compatible with previous        */
/*         versions.                                                     */
/*                                                                       */
/*  12/05 -RWG                                                           */
/*                                                                       */
/*   -Version 1.13:                                                      */
/*                                                                       */
/*     # Implemented restart capabilities.                               */
/*                                                                       */
/*       * Ouput options are reduced to 1) point time histories (nseis)  */
/*         with all_in_one=1 the only option, 2) time slice files, and   */
/*         3) SGT output.  Initialization of output is now done in       */
/*         init_outp().                                                  */
/*                                                                       */
/*  05/10 -RWG                                                           */
/*                                                                       */
/*   -Version 1.15:                                                      */
/*                                                                       */
/*     # Implemented buffering of xy-time slice output                   */
/*                                                                       */
/*  xx/11 -RWG                                                           */
/*                                                                       */
/*   -Version 3.1:                                                       */
/*                                                                       */
/*                                                                       */
/*                                                                       */
/*                                                                       */
/*************************************************************************/

#include "include.h"

int main(int ac,char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */

int nx1, nx2, globnx;  /* start and end ix plane, total number of points */
int ny1, ny2, globny;  /* start and end iy plane, total number of points */
int nz1, nz2, globnz;  /* start and end iz plane, total number of points */

struct nodeinfo ninfo;
int nproc, pnlen;
int nodeType = -1;
int segmentId, leftId, rightId;

int izord2 = -1;
int blen, off;
int vmodel_swapb = 0;
float *snd1, *rcv1, *snd2, *rcv2;   /* send/receive pointers */
int sndtag, rcvtag;

/* Rest of variables */
struct fdcoefs fdcoefs;
FILE *fpr, *fpw, *fopfile();
float *pvfield, **pvf, *pvptr[4];                    /* field varibles */
float *medfield, **medf, *medptr[4], *tmpmedfield;   /* media variables */
float *stfunc, *stfp;
float *pbnd;
float *tmp_ptr;
float vtmp, dt, h;
int nx, ny, nz, nt, iy, it, i, ib, bndflag;

int icnt, ystart, yend, itspt, iyspt;
size_t memlen;
int nbpad, print_tol;
int nt_p2;

float tzero;
int bfilt = 0;
float fhi = 0.0;
float flo = 1.0e+15;

float dep_val;
int intsrc = 0;
int nbounce = 1;
int ispan, span1, spanlen;
int span = 1;

char version[16], vid[16];
int vnc;

int size_float = sizeof(float);
struct modelstorage modstor;
struct media_input medinp;

int elas_only = 0;
int eflag = 0;
int eff_media();

float one = 1.0;
float pi = 3.141592654;
struct qvalues qval;
float vfzero = -1.0;
float qfzero = -1.0;
float fmin = -1.0;
float fmax = -1.0;

                /*  damping boundary parameters  */

float *randf, *rfptr;
int rfinc = 50;
int rfcnt, rfspan;

char string[1024];
int freesurf = 1;        /* freesurf = 0 -> absorb */
int order = ORDER;
float vsmin = -1.0;
int report = 1;		/* reporting interval for usage */

struct tms use1, use2;
int rtime0, rtime1;
struct tms tu1, dtu;
int rt1, drt;
float tsize1, tsize2;

                /*  dump output parameters  */

struct dump_output_info dumpinfo;

                /*  restart parameters  */

struct restartinfo rstart;

                /*  rlimit parameters  */

struct rlimit rlp;
int maxfiles;

                /*  model run parameters  */

struct runparams rpars;

                /*  output parameters  */

struct outputfields outp;

struct outputparams *outptr;
struct tsoutputparams *xyoutptr;
struct tsheader tshead;
struct seisoutputparams *soutptr;
struct seisheader *sheadptr;
struct sgtoutputparams *sgtoutptr;

char name[FILE_STR_LEN];
char logfile[FILE_STR_LEN];
char tmpdir[DIR_STR_LEN];
char tmpdirlist[DIR_STR_LEN];
char vmoddir[DIR_STR_LEN];
char logdir[DIR_STR_LEN];
char stfdir[DIR_STR_LEN];

float modelrot = 0.0;  /* rotation of y-axis from south (clockwise positive) */
float modellat = 34.0;    /* latitude of model origin */
float modellon = -118.0;    /* longitude of model origin */

/* source information */

int intface = 0;
struct interface intfac;

int plane = 0;
struct planesrc plsrc;

int isrc;
char slipout[256], xyzfrmt[8], tfrmt[8];
char sourcelocs[256];
char stype[256], stf_file[512];
struct pntsrcs srcs;

float xmom, ymom, zmom;
float tdelay = 0.0;

float invnh;
float normf = 1.0e+05; /* convert km to cm */

/* version number */

sprintf(vid,"3.0-mpi");

vnc = 0;
while(vid[vnc] != '\0')
   vnc++;

/*  MPI-RWG: Start-up MPI */

mpi_init(&ac,&av,&ninfo.nproc,&ninfo.segmentId,ninfo.procname,&pnlen);

ninfo.nproc_x = -1;
ninfo.nproc_y = -1;
ninfo.nproc_z = -1;
ninfo.min_nproc = 4;

/* defaults */

srcs.nsource = 1;
srcs.eqsrc = 0; 
srcs.expl = 0; 
srcs.psrc = 0;
srcs.bforce = 0;
srcs.dblcpl = 0;
srcs.pointmt = 0;
srcs.ffault = 0;
srcs.relative_slip = 1;
srcs.absolute_slip = 0;
srcs.area = 1.0;

slipout[0] = '\0';
sprintf(xyzfrmt,"10.4f");
sprintf(tfrmt,"9.2e");

sourcelocs[0] = '\0';

sprintf(stype,"gaus");
sprintf(name,"generic");
sprintf(outp.seisdir,".");
sprintf(vmoddir,".");
sprintf(tmpdir,"/tmp");
sprintf(logdir,".");
sprintf(stfdir,".");
stf_file[0] = '\0';
tmpdirlist[0] = '\0';

medinp.model_style = 0;
medinp.xorg = 0.0;
medinp.yorg = 0.0;
medinp.media_boundary = 0;
medinp.dampwidth = 0;
medinp.qbndmax = 50.0;        /* Qs is 50 just inside damping region */
medinp.qbndmin = 5.0;        /* Qs is 5 at computational boundary */

qval.qpfrac = -1.0;
qval.qsfrac = -1.0;
qval.qpqs_factor = -1.0;

/* default runparams */

/*
   geoproj=0: RWG spherical projection with local kmplat, kmplon
          =1: RWG great circle projection
          =2: UTM coordinate projection
*/

rpars.geoproj = 1;
rpars.xshift = -1.0e+15;
rpars.yshift = -1.0e+15;

modelrot = 0.0;  /* rotation of y-axis from south (clockwise positive) */
modellat = 34.0;    /* latitude of model origin */
modellon = -118.0;    /* longitude of model origin */

/* memory defaults */

modstor.maxmem = MEM_LIMIT;
modstor.intmem = 1;

setpar(ac,av);
mstpar("version","s",version);		/* version number */

getpar("name","s",name);
getpar("logdir","s",logdir);

makedir(logdir);
sprintf(logfile,"%s/%s-%.4d.rlog",logdir,name,ninfo.segmentId);
freopen(logfile,"w",stderr);

if(strncmp(version,vid,vnc) != 0)
   {
   fprintf(stderr,"***** This is version %s of emod3d-mpi:\n",vid);
   fprintf(stderr,"      Input parameters specify version %s\n",version);
   fprintf(stderr,"      Check input, exiting...\n");
   mpi_exit(-1);
   }

fprintf(stderr,"***** This is version %s of emod3d-mpi *****\n\n",vid);

/* storage parameters */
/*
getpar("span","d",&span);
getpar("intmem","d",&modstor.intmem);
*/
span = 1;
modstor.intmem = 1;

getpar("maxmem","d",&modstor.maxmem);
modstor.maxmem = modstor.maxmem*1000000; /* convert Mbytes -> bytes */

mstpar("nx","d",&nx);             /* model parameters */
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("h","f",&h);
mstpar("nt","d",&nt);
mstpar("dt","f",&dt);

getpar("min_nproc","d",&ninfo.min_nproc);
getpar("nproc_x","d",&ninfo.nproc_x);
getpar("nproc_y","d",&ninfo.nproc_y);
getpar("nproc_z","d",&ninfo.nproc_z);

/* MPI-RWG: compute info for distributed computation */

/*
   Set-up x,y,z indexing parameters:

   globny = global number of planes
   ny = local number of planes (including overlaps)
   ny1 = global index of first plane in this node (including overlap)
   ny2 = {(global index of last plane in this node) + 1} (including overlap)

   Thus, ny = ny2 - ny1
   
   Later, parameters 'iyleft' and 'iyright' will be set.  These are
   defined as:

   iyleft  = global index of first responsible plane in this node
   iyright = global index of last responsible plane in this node

   'Responsible plane' means this node is responsible for computing full
   time update on this plane.

   In general:

      iyleft = ny1 + 2
      iyright = ny2 - 3

   however, for HEAD node: iyleft = 0
   and for TAIL node: iyright = ny2 - 1 (=globny-1)

*/

globnx = nx;
globny = ny;
globnz = nz;

get_nproc(globnx,globny,globnz,&ninfo);		/* determines processor distribution */
get_nodetype_neighbor(&ninfo);			/* determines processor neighbors */

fprintf(stderr,"**** Running on node= %s\n",ninfo.procname);
fprintf(stderr,"          Total nproc= %-5d     nproc_x= %-5d nproc_y= %-5d nproc_z= %-5d\n",ninfo.nproc,ninfo.nproc_x,ninfo.nproc_y,ninfo.nproc_z);
fprintf(stderr,"               nodeId= %-5d\n",ninfo.segmentId);
fprintf(stderr,"            minusId_x= %-5d     plusId_x= %-5d\n",ninfo.minusId_x,ninfo.plusId_x);
fprintf(stderr,"            minusId_y= %-5d     plusId_y= %-5d\n",ninfo.minusId_y,ninfo.plusId_y);
fprintf(stderr,"            minusId_z= %-5d     plusId_z= %-5d\n\n",ninfo.minusId_z,ninfo.plusId_z);
fflush(stderr);

/* find nx1,nx2,nx
        ny1,ny2,ny
	nz1,nz2,nz

   elfag=0 -> balance ny on interior nodes
   elfag=1 -> balance ny across all nodes
*/

eflag = 1;
get_n1n2(eflag,globnx,ninfo.nproc_x,ninfo.procId_x,ninfo.minusId_x,ninfo.plusId_x,&nx1,&nx2,&nx);
get_n1n2(eflag,globny,ninfo.nproc_y,ninfo.procId_y,ninfo.minusId_y,ninfo.plusId_y,&ny1,&ny2,&ny);
get_n1n2(eflag,globnz,ninfo.nproc_z,ninfo.procId_z,ninfo.minusId_z,ninfo.plusId_z,&nz1,&nz2,&nz);

if(ny < 10)
   {
   fprintf(stderr,"**** local ny= %d, must be >= 10\n",ny);
   fprintf(stderr,"**** Need to reduce number of NY processors, exiting...\n");
   fflush(stderr);
   mpi_exit(-1);
   }

fprintf(stderr,"**** Distribution of model subsets:\n");
fprintf(stderr,"             local_nx= %-5d     (nx1= %-5d nx2= %-5d global_nx= %-5d)\n",nx,nx1,nx2,globnx);
fprintf(stderr,"             local_ny= %-5d     (ny1= %-5d ny2= %-5d global_ny= %-5d)\n",ny,ny1,ny2,globny);
fprintf(stderr,"             local_nz= %-5d     (nz1= %-5d nz2= %-5d global_nz= %-5d)\n",nz,nz1,nz2,globnz);
fprintf(stderr,"                    h= %f\n",h);
fprintf(stderr,"                    nt= %d\n",nt);
fprintf(stderr,"                    dt= %f\n\n",dt);
fflush(stderr);

set_model_storage(&modstor,&ninfo,nx,ny,nz,&span,&spanlen,name,tmpdir);
memlen = spanlen*nx*nz*sizeof(float);
span1 = spanlen - 1;

fprintf(stderr,"**** Memory stats for model (approximate):\n");
fprintf(stderr,"                     span= %-5d time steps\n",span);
fprintf(stderr,"              span length= %-5d planes\n",spanlen);
fprintf(stderr,"                     core= %8.1f Mb\n",(float)(N_MED_VARS+N_WAVE_VARS)*memlen/1.0e+06);
fprintf(stderr,"          total for model= %8.1f Mb\n\n",(float)((N_MED_VARS+N_WAVE_VARS))*nx*ny*nz*size_float/1.0e+06);
fflush(stderr);

if(modstor.intmem)
   {
   pvfield = modstor.pv_buf;
   medfield = modstor.med_buf;
   }
else
   {
   pvfield = (float *) check_malloc (N_WAVE_VARS*memlen);
   medfield = (float *) check_malloc (N_MED_VARS*memlen);
   }

pvf = (float **) check_malloc (spanlen*sizeof(float *));
medf = (float **) check_malloc (spanlen*sizeof(float *));
pbnd = (float *) check_malloc (3*(nx*nz + 2*(nx+nz))*size_float);

mpi_final("PROGRAM xxxxXXXX IS FINISHED");
exit(0);

/* MPI-RWG: end of distributed parameter set-up */

getpar("model_style","d",&medinp.model_style);
getpar("vmoddir","s",vmoddir);
getpar("vfzero","f",&vfzero);
getpar("qfzero","f",&qfzero);
getpar("fmin","f",&fmin);
getpar("fmax","f",&fmax);

getpar("elas_only","d",&elas_only);

if(medinp.model_style == 0)
   {
   mstpar("model","s",string);
   sprintf(medinp.modfile,"%s/%s",vmoddir,string);

   getpar("xorg","f",&medinp.xorg);
   getpar("yorg","f",&medinp.yorg);
   }
else if(medinp.model_style == 1)
   {
   mstpar("pmodfile","s",string);
   sprintf(medinp.pmodfile,"%s/%s",vmoddir,string);

   mstpar("smodfile","s",string);
   sprintf(medinp.smodfile,"%s/%s",vmoddir,string);

   mstpar("dmodfile","s",string);
   sprintf(medinp.dmodfile,"%s/%s",vmoddir,string);

   getpar("vmodel_swapb","d",&vmodel_swapb);

   getpar("qpfrac","f",&qval.qpfrac);
   getpar("qsfrac","f",&qval.qsfrac);
   getpar("qpqs_factor","f",&qval.qpqs_factor);

   if(qval.qpfrac <= 0.0 || qval.qsfrac <= 0.0)
      {
      mstpar("qmodfile","s",string);
      sprintf(medinp.qmodfile,"%s/%s",vmoddir,string);
      }
   }

getpar("media_boundary","d",&medinp.media_boundary);

getpar("order","d",&order);         /* spatial differencing order */
if(order != 2 && order != 4)
   order = ORDER;
getpar("izord2","d",&izord2);
if(izord2 > nz)
   izord2 = nz;
getpar("vsmin","f",&vsmin);

getpar("freesurf","d",&freesurf); /* boundary parameters */
if(freesurf)
   freesurf = 1;
getpar("dampwidth","d",&medinp.dampwidth); /* width of damping region */
if(medinp.dampwidth > 0)
   {
   getpar("qbndmax","f",&medinp.qbndmax);
   getpar("qbndmin","f",&medinp.qbndmin);
   getpar("rfinc","d",&rfinc);
   }

getpar("tzero","f",&tzero);       /* source parameters */
getpar("nbounce","d",&nbounce);

getpar("stype","s",stype);
if(strncmp("-1",stype,2) == 0)
   mstpar("stf_file","s",stf_file);

getpar("intsrc","d",&intsrc);
getpar("bfilt","d",&bfilt);
if(bfilt)
   {
   mstpar("flo","f",&flo);
   mstpar("fhi","f",&fhi);
   }

getpar("plane","d",&plane);
getpar("intface","d",&intface);

getpar("stfdir","s",stfdir);
makedir(stfdir);

if(!plane && !intface)
   {
   getpar("eqsrc","d",&srcs.eqsrc);
   getpar("expl","d",&srcs.expl);
   getpar("psrc","d",&srcs.psrc);
   getpar("bforce","d",&srcs.bforce);
   getpar("dblcpl","d",&srcs.dblcpl);
   getpar("pointmt","d",&srcs.pointmt);
   getpar("ffault","d",&srcs.ffault);
   getpar("nsource","d",&srcs.nsource);
   getpar("sourcelocs","s",sourcelocs);

   if(srcs.dblcpl == 0 && srcs.pointmt == 0 && srcs.bforce == 0)
      srcs.nsource = 1;

   if(srcs.dblcpl == 0 && srcs.pointmt == 0 && srcs.ffault == 0)
      {
      srcs.nstf = 1;

      srcs.rtime = (float *) check_malloc (sizeof(float));
      srcs.rtime[0] = tzero;

      srcs.stfindx = (int *) check_malloc (srcs.nsource*sizeof(int));
      for(i=0;i<srcs.nsource;i++)
         srcs.stfindx[i] = 0;

      srcs.ix = (int *) check_malloc ((srcs.nsource)*sizeof(int));
      srcs.iy = (int *) check_malloc ((srcs.nsource)*sizeof(int));
      srcs.iz = (int *) check_malloc ((srcs.nsource)*sizeof(int));

      if(srcs.nsource > 1 && sourcelocs[0] != '\0')
         {
	 fpr = fopfile(sourcelocs,"r");

         for(isrc=0;isrc<srcs.nsource;isrc++)
            fscanf(fpr,"%d %d %d",&srcs.ix[isrc],&srcs.iy[isrc],&srcs.iz[isrc]);

	 fclose(fpr);
         }
      else
         {
         mstpar("xsrc","vd",srcs.ix);
         mstpar("ysrc","vd",srcs.iy);
         mstpar("zsrc","vd",srcs.iz);
         }

      for(isrc=0;isrc<srcs.nsource;isrc++)
         {
	 /*SN: change the global system to local system */
	 srcs.iy[isrc] = srcs.iy[isrc] - ny1;

         if(freesurf && srcs.iz[isrc] < 1)
            srcs.iz[isrc] = 1;
         }

      if(srcs.bforce == 1)
	 {
         srcs.fxsrc = (float *) check_malloc ((srcs.nsource)*sizeof(float));
         srcs.fysrc = (float *) check_malloc ((srcs.nsource)*sizeof(float));
         srcs.fzsrc = (float *) check_malloc ((srcs.nsource)*sizeof(float));

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

	 xmom = ymom = zmom = 0.0;
         getpar("xmom","f",&xmom);
         getpar("ymom","f",&ymom);
         getpar("zmom","f",&zmom);

	 if(xmom != 0.0 || ymom != 0.0 || zmom != 0.0)
	    {
            srcs.fxsrc[0] = xmom;
            srcs.fysrc[0] = ymom;
            srcs.fzsrc[0] = zmom;
	    }
	 else
	    {
            mstpar("fxsrc","vf",srcs.fxsrc);
            mstpar("fysrc","vf",srcs.fysrc);
            mstpar("fzsrc","vf",srcs.fzsrc);
	    }

         if(srcs.nsource > 1 && sourcelocs[0] != '\0')
            {
            for(isrc=1;isrc<srcs.nsource;isrc++)
	       {
               srcs.fxsrc[isrc] = srcs.fxsrc[0];
               srcs.fysrc[isrc] = srcs.fysrc[0];
               srcs.fzsrc[isrc] = srcs.fzsrc[0];
	       }
            }

	 invnh = 1.0/(normf*h);
         for(isrc=0;isrc<srcs.nsource;isrc++)
            {
            srcs.fxsrc[isrc] = srcs.fxsrc[isrc]*invnh*invnh*invnh;
            srcs.fysrc[isrc] = srcs.fysrc[isrc]*invnh*invnh*invnh;
            srcs.fzsrc[isrc] = srcs.fzsrc[isrc]*invnh*invnh*invnh;

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

	    if(freesurf && srcs.iz[isrc] == 1)
	       {
	       srcs.fxsrc[isrc] = 2.0*srcs.fxsrc[isrc];
	       srcs.fysrc[isrc] = 2.0*srcs.fysrc[isrc];
	       }
            }
	 }
      }
   }
else
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.dblcpl = 0;
   srcs.pointmt = 0;
   srcs.ffault = 0;
   srcs.nsource = 0;
   }

if(plane)
   {
   mstpar("incidence","f",&plsrc.incidence);
   mstpar("azimuth","f",&plsrc.azimuth);
   mstpar("p","f",&plsrc.p);
   mstpar("sh","f",&plsrc.sh);
   mstpar("sv","f",&plsrc.sv);

   plsrc.bnd_zero_pad = NBND_ZERO_PLANE;
   getpar("bnd_zero_pad","d",&plsrc.bnd_zero_pad);

   plsrc.ixpzero = NBND_ZERO_PLANE;
   plsrc.iypzero = NBND_ZERO_PLANE;
   plsrc.izpzero = NBND_ZERO_PLANE;
   getpar("ixpzero","d",&plsrc.ixpzero);
   getpar("iypzero","d",&plsrc.iypzero);
   getpar("izpzero","d",&plsrc.izpzero);

   plsrc.intercept = -1;
   getpar("intercept","d",&plsrc.intercept);

   srcs.nsource = 0;
   }

if(intface)
   {
   get_infc_par(&intfac);
   switch(intfac.mode)
      {
      case FAULT:
         get_fault_par(&intfac);
         break;
      }
   }
else if(srcs.dblcpl)
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.pointmt = 0;
   srcs.ffault = 0;
   srcs.relative_slip = srcs.absolute_slip = 0;

   get_dc_par(&srcs,&dt,&tzero,ny1);
   }
else if(srcs.pointmt)
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.dblcpl = 0;
   srcs.ffault = 0;
   srcs.relative_slip = srcs.absolute_slip = 0;

   get_pointmt_par(&srcs,&dt,&tzero,ny1);
   }
else if(srcs.ffault == 1)
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.pointmt = 0;

   getpar("slipout","s",slipout);
   getpar("xyzfrmt","s",xyzfrmt);
   getpar("tfrmt","s",tfrmt);

   get_ff_par(&srcs,&dt,&h,freesurf,&tzero);

/*
   RWG 12/19/05:
     All the sub-fault specifications are in global system
     thus, it is necessary to transform them into local systems.

     If relative moments specified, do global to local conversion here,
     otherwise it will be done later in slip2mom().
*/

   if(srcs.absolute_slip == 0 && srcs.absolute_slip == 0)
      globalSource2Local(&srcs,ny1);
      
/*  Finite-fault is now represented as a distribution of double-couples */

   srcs.ffault = 0;
   srcs.dblcpl = 1;
   }
else if(srcs.ffault == 2) /* read SRF format -> moment-tensors */
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.dblcpl = 0;

   getpar("slipout","s",slipout);
   getpar("xyzfrmt","s",xyzfrmt);
   getpar("tfrmt","s",tfrmt);

   mstpar("faultfile","s",srcs.faultfile);
   }

/* other run parameters */

getpar("modelrot","f",&modelrot);
getpar("modellat","f",&modellat);
getpar("modellon","f",&modellon);
getpar("xshift","f",&rpars.xshift);
getpar("yshift","f",&rpars.yshift);
getpar("geoproj","d",&rpars.geoproj);

/* output parameters */
getpar("report","d",&report);
getpar("seisdir","s",outp.seisdir);

outp.sgtout = 0;
getpar("sgtout","d",&outp.sgtout);  /* moment tensor at individual points */
if(outp.sgtout == 1)
   {
   mstpar("sgtcords","s",outp.sgtcords);

   sgtoutptr = &(outp.sgtpnts);
   sgtoutptr->tinc = 1;
   getpar("sgt_tinc","d",&(sgtoutptr->tinc));
   }

outp.nseis = 0;
getpar("nseis","d",&outp.nseis);         /*  individual points */
if(outp.nseis)
   {
   mstpar("seiscords","s",outp.seiscords);

   soutptr = &(outp.spnts);
   soutptr->tinc = 1;
   getpar("seistinc","d",&(soutptr->tinc));
   }

outp.ts_xy = 0;
outp.ts_xz = 0;
outp.ts_yz = 0;
getpar("ts_xy","d",&outp.ts_xy);         /*  time slice parameters  */
getpar("ts_xz","d",&outp.ts_xz);
getpar("ts_yz","d",&outp.ts_yz);

outp.dtts = 1;
outp.dxts = 1;
outp.dyts = 1;
outp.dzts = 1;
getpar("dtts","d",&outp.dtts);
getpar("dxts","d",&outp.dxts);
getpar("dyts","d",&outp.dyts);
getpar("dzts","d",&outp.dzts);

outp.ix_ts = 1;
outp.iy_ts = 1;
outp.iz_ts = 1;
getpar("ix_ts","d",&outp.ix_ts);
getpar("iy_ts","d",&outp.iy_ts);
getpar("iz_ts","d",&outp.iz_ts);

/* restart parameters */

rstart.enable_flag = 0;
rstart.read_flag = 0;
rstart.itinc = -999999;
rstart.it = -999999;
sprintf(rstart.dir,".");

getpar("enable_restart","d",&rstart.enable_flag);
getpar("read_restart","d",&rstart.read_flag);

if(rstart.enable_flag == 1 || rstart.read_flag == 1)
   {
   getpar("restartdir","s",rstart.dir);
   makedir(rstart.dir);
   }

if(rstart.enable_flag == 1)
   getpar("restart_itinc","d",&rstart.itinc);

if(rstart.read_flag == 1)
   {
   getpar("restartname","s",rstart.name);

   sprintf(rstart.readfile,"%s/%s_rst-%.4d.e3d",rstart.dir,rstart.name,segmentId);
   rstart.fdr = opfile(rstart.readfile);
   }

/* dump output parameters */

dumpinfo.enable_flag = 0;
dumpinfo.itinc = -999999;

getpar("enable_output_dump","d",&dumpinfo.enable_flag);

strcpy(dumpinfo.main_dir,outp.seisdir); /* intialize with something just in case */
strcpy(dumpinfo.name,name);

if(dumpinfo.enable_flag == 1)
   {
   if(rstart.enable_flag == 1)
      strcpy(dumpinfo.local_restartdir,rstart.dir);
   else
      sprintf(dumpinfo.local_restartdir,"RESTART_NOT_USED");

/*
   Initialize here, but check in init_outp() to see if files need to
   be copied.
*/
   
   sprintf(dumpinfo.local_outputdir,"NO_SEISMOS_OUPUT");

   mstpar("dump_itinc","d",&dumpinfo.itinc);
   mstpar("main_dump_dir","s",dumpinfo.main_dir);
   makedir(dumpinfo.main_dir);

   sprintf(dumpinfo.command,"DEFAULT");
   getpar("dump_command","s",dumpinfo.command);
   if(strcmp(dumpinfo.command,"DEFAULT") == 0)
      sprintf(dumpinfo.command,"\\cp -rp");
   }

endpar();

fflush(stderr);

getrlimit(RLIMIT_FSIZE,&rlp);
getrlimit(RLIMIT_NOFILE,&rlp);
maxfiles = rlp.rlim_max;

if(report)
   {
   rtime1 = times(&use1);
   rtime0 = rtime1;
   }

if(medinp.media_boundary == 1)
   {
   fprintf(stderr,"**** Media values centered on NULL node\n");
   fprintf(stderr,"     -topmost (thinnest possible) layer has depth of 'h'\n");
   fprintf(stderr,"     -all other layers have a depth of 'k*h' [k=integer]\n\n");
   medinp.media_boundary = NULL_NODE_CENTERED;
   }
else if(medinp.media_boundary == 2)
   {
   fprintf(stderr,"**** Media values centered on VX node\n");
   fprintf(stderr,"     -topmost (thinnest possible) layer has depth of 'h'\n");
   fprintf(stderr,"     -all other layers have a depth of 'k*h' [k=integer]\n\n");
   fprintf(stderr,"     -media averaging is done with arithmetic operators (KBO)\n\n");
   medinp.media_boundary = KBO_AT_VXNODE;
   }
else if(medinp.media_boundary == 3)
   {
   fprintf(stderr,"**** Media values centered on VZ node\n");
   fprintf(stderr,"     -topmost (thinnest possible) layer has depth of 'h'\n");
   fprintf(stderr,"     -all other layers have a depth of 'k*h' [k=integer]\n\n");
   fprintf(stderr,"     -media averaging is done with arithmetic operators (KBO)\n\n");
   medinp.media_boundary = KBO_AT_VZNODE;
   }
else if(medinp.media_boundary == 4)
   {
   fprintf(stderr,"**** Media values centered on VZ node\n");
   fprintf(stderr,"     -topmost (thinnest possible) layer has depth of 'h'\n");
   fprintf(stderr,"     -all other layers have a depth of 'k*h' [k=integer]\n\n");
   fprintf(stderr,"     -media averaging is done with harmonic operators\n\n");
   medinp.media_boundary = HARMONIC_AT_VZNODE;
   }
else if(medinp.media_boundary == -999)
   {
   fprintf(stderr,"**** CAUTION!  No media averaging will be done.\n");
   fprintf(stderr,"     This option should only be used as a test.\n\n");
   }
else 
   {
   fprintf(stderr,"**** Media values centered on normal stress node,\n");
   fprintf(stderr,"     -topmost (thinnest possible) layer has depth of '0.5*h'\n");
   fprintf(stderr,"     -all other layers have a depth of '(k+0.5)*h' [k=integer]\n\n");
   medinp.media_boundary = NORMAL_STRESS_CENTERED;
   }

fprintf(stderr,"**** Initializing run parameters:\n\n");
fflush(stderr);

if(bfilt)
   {
   it = 1.0/(dt*flo);
   tdelay = it*dt;
   }

rpars.segmentId = segmentId;
rpars.nodeType = nodeType;
rpars.nx = nx;
rpars.globny = globny;
rpars.nz = nz;
rpars.nt = nt;
rpars.ny1 = ny1;
rpars.ny2 = ny2;
rpars.span = span;
rpars.freesurf = freesurf;
rpars.h = h;
rpars.dt = dt;
rpars.modelrot = modelrot;
rpars.modellon = modellon;
rpars.modellat = modellat;
rpars.xmom = xmom;
rpars.ymom = ymom;
rpars.zmom = zmom;
rpars.tdelay = tdelay;

init_modelc(&rpars);

fprintf(stderr,"**** Initializing media values:\n");
fflush(stderr);

/* use paired MPI_Barrier() calls so nodes don't all hit the disk at same time
*/
for(ib=0;ib<segmentId;ib++)
   MPI_Barrier(MPI_COMM_WORLD);

init_media(&medinp,medfield,&rpars,&qval,&modstor,vmodel_swapb);

/* compliment to above MPI_Barrier() call
*/
for(ib=segmentId;ib<nproc;ib++)
   MPI_Barrier(MPI_COMM_WORLD);

fprintf(stderr,"**** Initializing wave field values:\n");
fflush(stderr);

if(plane != 1)
   {
   if(modstor.intmem)
      init_field_val(1.0e-20,pvfield,ny*N_WAVE_VARS*nx*nz);
   else
      {
      init_field_val(1.0e-20,pvfield,N_WAVE_VARS*nx*nz);

      for(iy=0;iy<ny;iy++)
         rite_model(&modstor,pvfield,WAVE_FIELD);
      }
   }
else
   {
   init_model_seek(&modstor,MEDIA_FIELD);
   medf[0] = medfield;
   medf[0] = reed_model(&modstor,medf[0],MEDIA_FIELD);
   init_model_seek(&modstor,MEDIA_FIELD);

   plsrc.srcl2m = medf[0][nx*nz - 1];
   plsrc.srclam = medf[0][2*nx*nz - 1];
   plsrc.srcmu = 1.0/medf[0][5*nx*nz - 1];

   plsrc.svel = sqrt(plsrc.srcmu/medf[0][8*nx*nz-1]);
   plsrc.pvel = sqrt(plsrc.srcl2m/medf[0][8*nx*nz-1]);

   pvptr[0] = pvfield;

   for(iy=0;iy<ny;iy++)
      {
      if(modstor.intmem)
	 pvptr[0] = pvfield + iy*N_WAVE_VARS*nx*nz;

      init_field(pvptr[0],N_WAVE_VARS*nx*nz);

/* 04/08/05 RWG: Use globny and glob_iy=ny1+iy for plane wave specification */
/* STILL NEED TO TEST */
      init_plane_src(pvptr[0],nx,globny,nz,ny1+iy,&h,&tzero,&dt,&plsrc,nbounce);

      rite_model(&modstor,pvptr[0],WAVE_FIELD);
      }
   }

if(intface)
   {
   fprintf(stderr,"**** Initializing interfacing wave field:\n");
   intfac.go = 1;

   switch(intfac.mode)
      {
      case FAULT:
         sprintf(intfac.gfname,"%s/%s_gf",tmpdir,name);
	 init_fault(nx,ny,nz,nt,nt_p2,stfunc,&h,&dt,&intfac,medinp.medsten,medinp.medprofs,medinp.nprof,medfield,freesurf);
         break;
      }

   intfac.ny = ny;
   intfac.goplane = (int *) check_malloc (ny*sizeof(int));
   intfac.newgo = (int *) check_malloc (ny*sizeof(int));
   for(iy=0;iy<ny;iy++)
      {
      intfac.goplane[iy] = 1;
      intfac.newgo[iy] = 0;
      }
   } /*SN: interface!=0 */
else
   {
   intfac.go = 0;

   if(srcs.dblcpl)
      {
      fprintf(stderr,"**** Initializing double-couple sources:\n");

      srcs.doublec.h = h;

      if(srcs.relative_slip || srcs.absolute_slip == 1)
         {
	 if(slipout[0] != '\0')
	    sprintf(slipout,"%s-%.4d",slipout,segmentId);

         slip2mom(&srcs,&h,&modstor,medfield,nx,ny1,ny2,ny,nz,slipout,&dt,xyzfrmt,tfrmt,nodeType);
         }

      init_dc(&srcs,&dt);
      intsrc = 0;

      fprintf(stderr,"\n");
      }

   if(srcs.pointmt)
      {
      fprintf(stderr,"**** Initializing point moment-tensor sources:\n");

      init_pointmt(&srcs,&dt,&h);
      intsrc = 0;

      fprintf(stderr,"\n");
      }

   if(srcs.ffault == 2)
      {
      fprintf(stderr,"**** Initializing moment-tensor sources:\n");
      fflush(stderr);

      if(slipout[0] != '\0')
         sprintf(slipout,"%s-%.4d",slipout,segmentId);

      get_srf_par(&srcs,&rpars,bfilt,&flo,&fhi,&tdelay);
      init_mt(&srcs,&h,&modstor,medfield,nx,ny1,ny2,ny,nz,slipout,&dt,xyzfrmt,tfrmt,nodeType);

      if(bfilt)
         {
         fprintf(stderr,"\n");
         fprintf(stderr,"     source function has been filtered:\n");
         fprintf(stderr,"       -STF bandpass: %6.2f Hz < f < %6.2f Hz\n",fhi,flo);
         fprintf(stderr,"       -a time delay of %10.5e sec has been applied to preserve causality\n",tdelay);
         fflush(stderr);
	 }
      }
   else
      {
      fprintf(stderr,"**** Initializing source time function:\n");

      stfunc = (float *) check_malloc (srcs.nstf*nt*sizeof(float));
      init_field_val(1.0e-20,stfunc,srcs.nstf*nt);
      for(i=0;i<srcs.nstf;i++)
         {
         stfp = stfunc + i*nt;
         getsource(stfp,&dt,nt,&srcs.rtime[i],nbounce,intsrc,stype,name,bfilt,&flo,&fhi,i,segmentId,&tdelay,stfdir,stf_file);
         }
      }

   fprintf(stderr,"\n");
   }

if(rstart.enable_flag == 1 || rstart.read_flag == 1)
   fprintf(stderr,"     Restart mode invoked:\n\n");

if(rstart.read_flag == 1)
   {
   reed(rstart.fdr,&rstart.it,sizeof(int));

   fprintf(stderr,"     -reading wavefield at time step= %d from file:\n",rstart.it);
   fprintf(stderr,"      %s\n\n",rstart.readfile);
   fflush(stderr);

   reed(rstart.fdr,&rstart.rpars,sizeof(struct runparams));
   check_rpars(rstart.readfile,&rstart.rpars,&rpars);
   }

if(rstart.enable_flag == 1)
   {
   sprintf(rstart.dumpfile,"NOT_STARTED_YET");

   if(rstart.itinc < 1)
      rstart.itinc = 100;

   i = rstart.itinc%span;
   rstart.itinc = rstart.itinc + i;

   fprintf(stderr,"     -restart wavefield output at every %d time steps\n\n",rstart.itinc);
   fflush(stderr);
   }

if(dumpinfo.enable_flag == 1)
   {
   if(dumpinfo.itinc < 1)
      dumpinfo.itinc = 1000;

   i = dumpinfo.itinc%span;
   dumpinfo.itinc = dumpinfo.itinc + i;
   }

fprintf(stderr,"**** Computing effective media parameters:\n");
fflush(stderr);

eflag = eff_media(&medinp,medfield,nx,globny,ny1,ny2,nz,&modstor,nodeType);

if(elas_only == 1)
   fprintf(stderr,"\n**** 'Elastic Only' mode invoked, no attenuation operators will be used.\n\n");
else
   {
   fprintf(stderr,"**** Computing memory variable coefficients:\n");
   mvar_coefs(&modstor,medfield,nx,ny,nz,&dt,nt,freesurf,&vfzero,&qfzero,&fmin,&fmax,eflag,nodeType,ny1);
   /*
   smv_coefs(&modstor,medfield,nx,ny,nz,&dt,nt,freesurf,&vfzero,&qfzero,&fmin,&fmax,eflag);
   */
   }

fprintf(stderr,"**** Computing FD coefficients (checking vmin & vmax):\n");
set_vminvmax_reedmod(&modstor,medfield,&fdcoefs,nx,ny,nz,&vsmin);

/* MPI-RWG: find global max(izord2), vmin & vmax */

vtmp = fdcoefs.vmax;
mpi_global_val(&vtmp,&fdcoefs.vmax,"Vmax",1,MPI_FLOAT,MPI_MAX);

vtmp = fdcoefs.vmin;
mpi_global_val(&vtmp,&fdcoefs.vmin,"Vmin",1,MPI_FLOAT,MPI_MIN);

set_izord2(&modstor,medfield,&fdcoefs,nx,ny,nz,&vsmin);

if(izord2 < 0)
   izord2 = fdcoefs.izord2;

mpi_global_val(&izord2,&fdcoefs.izord2,"izord2",1,MPI_INT,MPI_MAX);

setcoefs(order,&fdcoefs,&dt,&h,&vsmin);
fflush(stderr);

if(elas_only == 0 || medinp.dampwidth > 0)
   {
   rfcnt = 0;
   rfspan = 6;

   i = rfinc%span;
   rfinc = rfinc + i;
   if(rfinc <= 0)
      rfinc = 1;

   randf = (float *) check_malloc (rfspan*nx*nz*sizeof(float));
   init_field_val(1.0e-20,randf,rfspan*nx*nz);
   }

              /* initialize output files */

fprintf(stderr,"**** Initializing output files:\n\n");
fflush(stderr);

init_outp(&outp,&rpars,name,&srcs,&rstart,&dumpinfo);

         /* set field pointers */

for(ispan=0;ispan<spanlen;ispan++)
   {
   pvf[ispan] = pvfield + N_WAVE_VARS*ispan*nx*nz;
   medf[ispan] = medfield + N_MED_VARS*ispan*nx*nz;
   }

if(rstart.read_flag == 1)
   {
   for(ib=0;ib<segmentId;ib++)    /* use paired MPI_Barrier() calls so nodes don't all hit the disk at same time */
      MPI_Barrier(MPI_COMM_WORLD);

   get_restart(&rstart,&modstor,pvfield,nx,ny,nz);

   for(ib=segmentId;ib<nproc;ib++)    /* compliment to above MPI_Barrier() call */
      MPI_Barrier(MPI_COMM_WORLD);

/* if randfield is used, roll through rfpointers to ensure exact replication */
   if((elas_only == 0 || medinp.dampwidth > 0))
      {
      for(it=0;it<=rstart.it;it=it+span)
         {
         if((it+span)%(rfinc) == 0)
            {
	    for(iy=0;iy<ny;iy++)
	       {
               rfcnt++;
               if(rfcnt == rfspan)
	          rfcnt = 0;
	       }
            }
         }
      }
   }

             /*  Initialize usage status report  */

if(report)
   init_report(rtime0,&rtime1,&use2,&use1,&dtu,&drt,&tsize2,&tsize1,report);

for(it=0;it<nt;it=it+span)
   {
   init_model_seek(&modstor,WAVE_FIELD);
   init_model_seek(&modstor,MEDIA_FIELD);

   if(it > rstart.it)
      {
      for(iy=0;iy<spanlen;iy++)   /* roll in: get field and media values */
         {
         pvf[iy]  = reed_model(&modstor,pvf[iy],WAVE_FIELD);
         medf[iy] = reed_model(&modstor,medf[iy],MEDIA_FIELD);
         }

      for(ispan=0;ispan<span;ispan++) /* roll in at iy=0 boundary */
         {
         itspt = it + ispan;                 /* time step */
         yend = spanlen - 4*(ispan + 1) + 1; /* y limit for 'ispan' update */
         for(iy=0;iy<yend;iy++)    /* do velocity update, iy -> plane 0 */
            {
            for(icnt=0;icnt<4;icnt++)  /* set pointers */
	       {
	       pvptr[icnt] = pvf[icnt+iy];
	       medptr[icnt] = medf[icnt+iy];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iy+2) == srcs.iy[isrc])          /* add in body force */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,FSRC);
	       }

	    if(intfac.go)     /* set interfacing parameters */
	       {
	       intfac.iy = iy + 3;
	       intfac.it = itspt;
	       }

	    bndflag = 0;
  	    if(nodeType == PARL_HEAD && iy == 0)   
               bndflag = YZERO;

            tstepv(nx,nz,pvptr,medptr,pbnd,&fdcoefs,freesurf,bndflag,&intfac);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy+2,itspt,0);
*/

	    if(plane == 2)
               add_plane(stfunc,pvptr,nx,nz,itspt,&plsrc);

            wrout(&outp,pvptr[1],medptr[1],iy+1,ny1+iy+1,outp.nseis,nx,nz,itspt);
	    if(iy == 0)
               wrout(&outp,pvptr[0],medptr[0],iy,ny1+iy,outp.nseis,nx,nz,itspt);

	    } /* end of for(iy-0;iy<yend;iy++) */

/* update the boundary information by communication with the next machine*/
/* notice the difference of out/in planes 
	  between this code and the head code 
	  head: 0,2;  body 2,0
	  head: 1,3;  body 3,1
*/
	
         if(nodeType != PARL_HEAD)
            {
            rt1 = times(&tu1);

	    snd1 = pvf[2];
	    rcv1 = pvf[0];
	    snd2 = pvf[3];
	    rcv2 = pvf[1];
	    sndtag = itspt + SEND_VEL_LEFT;
	    rcvtag = itspt + RECV_VEL_LEFT;
            blen = nx*nz*3*sizeof(float);

            mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,leftId,sndtag,rcvtag);

            dtimes(tu1,rt1,&dtu,&drt,&tsize2,4*blen);
            }

         yend = yend - 2;
         for(iy=0;iy<yend;iy++)    /* do stress update */
            {
            for(icnt=0;icnt<4;icnt++)  /* set pointers */
	       {
	       pvptr[icnt] = pvf[icnt+iy];
	       medptr[icnt] = medf[icnt+iy];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iy+2) == srcs.iy[isrc])          /* add in pressure source */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,PSRC);
	       }

	    bndflag = 0;
	    if(nodeType == PARL_HEAD && iy == 0)   
               bndflag = YZERO;

            tstepp(nx,nz,pvptr,medptr,&fdcoefs,freesurf,bndflag,elas_only);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy+2,itspt,1);
*/
	    } /* end of stress update */

/* communication information */	
/* notice the difference of out/in planes 
	  between this code and the head code 
	  head: 0,2;  body 2,0
	  head: 1,3;  body 3,1
*/

         if(nodeType != PARL_HEAD)
            {
            rt1 = times(&tu1);

            off = nx*nz*3;
	    snd1 = pvf[2] + off;
	    rcv1 = pvf[0] + off;
	    snd2 = pvf[3] + off;
	    rcv2 = pvf[1] + off;
	    sndtag = itspt + SEND_TAU_LEFT;
	    rcvtag = itspt + RECV_TAU_LEFT;
            blen = nx*nz*(N_WAVE_VARS-3)*sizeof(float);

            mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,leftId,sndtag,rcvtag);

            dtimes(tu1,rt1,&dtu,&drt,&tsize2,4*blen);
            }
         } /* end of roll in */

   /* if damping boundaries are used, add random field to prevent underflow */

      if((elas_only == 0 || medinp.dampwidth > 0) && (it+span)%(rfinc) == 0)
         {
         rfptr = randf + rfcnt*nx*nz;
         add_random_field(pvf[0],rfptr,nx,nz);

         rfcnt++;
         if(rfcnt == rfspan)
	    rfcnt = 0;
         }

      rite_model(&modstor,pvf[0],WAVE_FIELD);      /* store updated plane */

		   /* reset field and media pointers */

      tmp_ptr = pvf[0];
      for(ispan=0;ispan<span1;ispan++)
         pvf[ispan] = pvf[ispan+1];
      pvf[span1] = tmp_ptr;

      tmp_ptr = medf[0];
      for(ispan=0;ispan<span1;ispan++)
         medf[ispan] = medf[ispan+1];
      medf[span1] = tmp_ptr;

		  /* do interior planes */

      for(iy=spanlen;iy<ny-1;iy++)
         {
         /* get next plane of field variables and media parameters */

         pvf[span1]  = reed_model(&modstor,pvf[span1],WAVE_FIELD);
         medf[span1] = reed_model(&modstor,medf[span1],MEDIA_FIELD);

         for(ispan=0;ispan<span;ispan++) /* do 'span' time updates */
	    {
            itspt = it + ispan;       /* time step */
            iyspt = iy - 4*ispan;     /* y index of plane 3 */

            for(icnt=0;icnt<4;icnt++)  /* set pointers */
	       {
	       pvptr[icnt] = pvf[span1 - (3-icnt) - 4*ispan];
	       medptr[icnt] = medf[span1 - (3-icnt) - 4*ispan];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iyspt-1) == srcs.iy[isrc])          /* add in body force */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,FSRC);
	       }

	    if(intfac.go)     /* set interfacing parameters */
	       {
	       intfac.iy = iyspt;
	       intfac.it = itspt;
	       }

            /* do velocity update */
            tstepv(nx,nz,pvptr,medptr,pbnd,&fdcoefs,freesurf,0,&intfac);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy+2,itspt,2);
*/

	    if(plane == 2)
               add_plane(stfunc,pvptr,nx,nz,itspt,&plsrc);

            wrout(&outp,pvptr[1],medptr[1],iyspt-2,ny1+iyspt-2,outp.nseis,nx,nz,itspt); /* write output */

            for(icnt=0;icnt<4;icnt++)  /* reset pointers 2 planes back */
	       {
	       pvptr[icnt] = pvf[span1 - (3-icnt) - 4*ispan - 2];
	       medptr[icnt] = medf[span1 - (3-icnt) - 4*ispan - 2];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iyspt-3) == srcs.iy[isrc])      /* add in pressure source */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,PSRC);
	       }

            /* do pressure update -> bndflag=0 */
            tstepp(nx,nz,pvptr,medptr,&fdcoefs,freesurf,0,elas_only);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy,itspt,3);
*/
	    }

   /* if damping boundaries are used, add random field to prevent underflow */

         if((elas_only == 0 || medinp.dampwidth > 0) && (it+span)%(rfinc) == 0)
            {
            rfptr = randf + rfcnt*nx*nz;
            add_random_field(pvf[0],rfptr,nx,nz);

            rfcnt++;
            if(rfcnt == rfspan)
	       rfcnt = 0;
            }

         rite_model(&modstor,pvf[0],WAVE_FIELD);         /* store updated planes */

		   /* reset field and media pointers */

         tmp_ptr = pvf[0];
         for(ispan=0;ispan<span1;ispan++)
            pvf[ispan] = pvf[ispan+1];
         pvf[span1] = tmp_ptr;

         tmp_ptr = medf[0];
         for(ispan=0;ispan<span1;ispan++)
            medf[ispan] = medf[ispan+1];
         medf[span1] = tmp_ptr;
         }

   /* get last set of field and media planes at iy=ny-1 */

      pvf[span1] = reed_model(&modstor,pvf[span1],WAVE_FIELD);
      medf[span1] = reed_model(&modstor,medf[span1],MEDIA_FIELD);

      for(ispan=0;ispan<span;ispan++) /* roll out at iy=ny-1 boundary */
         {
         itspt = it + ispan;          /* time step */
         ystart = (ny-1) - 4*ispan;   /* y start for 'ispan' update */
         for(iy=ystart;iy<ny;iy++)    /* do velocity update, iy -> plane 3 */
            {
            for(icnt=0;icnt<4;icnt++)  /* set pointers */
	       {
	       pvptr[icnt] = pvf[span1 - (3-icnt) - (ny-1-iy)];
	       medptr[icnt] = medf[span1 - (3-icnt) - (ny-1-iy)];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iy-1) == srcs.iy[isrc])          /* add in body force */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,FSRC);
	       }

	    if(intfac.go)     /* set interfacing parameters */
	       {
	       intfac.iy = iy;
	       intfac.it = itspt;
	       }

/* boundary condition flag at iy=ny-1 */

	    bndflag = 0;
	    if(nodeType == PARL_TAIL && iy == ny-1)
	       bndflag = YN;

            tstepv(nx,nz,pvptr,medptr,pbnd,&fdcoefs,freesurf,bndflag,&intfac);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy+2,itspt,4);
*/

	    if(plane == 2)
               add_plane(stfunc,pvptr,nx,nz,itspt,&plsrc);

            wrout(&outp,pvptr[1],medptr[1],iy-2,ny1+iy-2,outp.nseis,nx,nz,itspt);
	    if(iy == ny-1)
	       {
               wrout(&outp,pvptr[2],medptr[2],iy-1,ny1+iy-1,outp.nseis,nx,nz,itspt);
               wrout(&outp,pvptr[3],medptr[3],iy,ny1+iy,outp.nseis,nx,nz,itspt);
	       }

            }

/** now pvptr[2] and pvptr[3] need to be updated
  while pvptr[0] and pvptr[1] need to be sent to the next machine
**/

         if(nodeType != PARL_TAIL)
            {
            rt1 = times(&tu1);

	    snd1 = pvptr[0];
	    rcv1 = pvptr[2];
	    snd2 = pvptr[1];
	    rcv2 = pvptr[3];
	    sndtag = itspt + SEND_VEL_RIGHT;
	    rcvtag = itspt + RECV_VEL_RIGHT;
            blen = nx*nz*3*sizeof(float);

            mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,rightId,sndtag,rcvtag);

            dtimes(tu1,rt1,&dtu,&drt,&tsize2,4*blen);
            }

         ystart = ystart - 2;
         for(iy=ystart;iy<ny;iy++)    /* do pressure update */
            {
            for(icnt=0;icnt<4;icnt++)  /* set pointers */
	       {
	       pvptr[icnt] = pvf[span1 - (3-icnt) - (ny-1-iy)];
	       medptr[icnt] = medf[span1 - (3-icnt) - (ny-1-iy)];
	       }

	    for(isrc=0;isrc<srcs.nsource;isrc++)
	       {
               if((iy-1) == srcs.iy[isrc])          /* add in pressure source */
                  add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,itspt,&srcs,PSRC);
	       }

	/* boundary condition flag at iy=ny-1 */
	/* SN: conditional compilation */

	    bndflag = 0;
	    if(nodeType==PARL_TAIL && iy == ny-1) 
	       bndflag = YN;

            tstepp(nx,nz,pvptr,medptr,&fdcoefs,freesurf,bndflag,elas_only);

/*
check_nan(pvptr[2],medptr[2],nx,nz,ny1-iy+2,itspt,5);
*/
	    } /* end of stress update */

/* exchange stress with machine on right */

         if(nodeType != PARL_TAIL)
            {
            rt1 = times(&tu1);

            off = nx*nz*3;
	    snd1 = pvptr[0] + off;
	    rcv1 = pvptr[2] + off;
	    snd2 = pvptr[1] + off;
	    rcv2 = pvptr[3] + off;
	    sndtag = itspt + SEND_TAU_RIGHT;
	    rcvtag = itspt + RECV_TAU_RIGHT;
            blen = nx*nz*(N_WAVE_VARS-3)*sizeof(float);

            mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,rightId,sndtag,rcvtag);

            dtimes(tu1,rt1,&dtu,&drt,&tsize2,4*blen);
            }
         } /* end of  roll out */

/* XXXX  */
/* Do all of MPI communication here. Use array modstor.pv_buf */

   /* if seismogram output is used, write temporary buffer to file */

      if(outp.nseis)
         {
         soutptr = &(outp.spnts);

         soutptr->flushbuf = 1;
         wrout(&outp,pvptr[3],medptr[3],-1,-1,-1,-1,-1,-1);
         soutptr->flushbuf = 0;
         }

   /* if xy time slice output is used, write temporary buffer to file */

      if(outp.ts_xy)
         {
         xyoutptr = &(outp.xyslice);

         xyoutptr->flushbuf = 1;
         wrout(&outp,pvptr[3],medptr[3],-1,-1,-1,-1,-1,-1);
         xyoutptr->flushbuf = 0;
         }

   /* if moment tensor output is used, write temporary buffer to output file */

      if(outp.sgtout)
         {
         sgtoutptr = &(outp.sgtpnts);

         sgtoutptr->flushbuf = 1;
         wrout(&outp,pvptr[3],medptr[3],-1,-1,-1,-1,-1,-1);
         sgtoutptr->flushbuf = 0;
         }

   /* if damping boundaries are used, add random field to prevent underflow */

      if((elas_only == 0 || medinp.dampwidth > 0) && (it+span)%(rfinc) == 0)
         {
         for(ispan=0;ispan<spanlen;ispan++)
            {
            rfptr = randf + rfcnt*nx*nz;
            add_random_field(pvf[ispan],rfptr,nx,nz);

            rfcnt++;
            if(rfcnt == rfspan)
	       rfcnt = 0;
            }
         }

      for(ispan=0;ispan<spanlen;ispan++)    /* store updated wave field planes */
         rite_model(&modstor,pvf[ispan],WAVE_FIELD);

      if(intfac.go)     /* check to continue interfacing */
         check_infc_goflag(&intfac);

      if(rstart.enable_flag == 1 && ((it+span)%(rstart.itinc)==0))
         dump_restart(&rstart,name,segmentId,&rpars,&modstor,pvfield,it,nx,ny,nz);

      if(dumpinfo.enable_flag == 1 && (((it+span)%(dumpinfo.itinc)==0) || (it+span) >= nt))
         dump_files(&outp,segmentId,&dumpinfo);

      }  /* end of restart check */

   if(report && ((it+span)%(report)==0))
      print_report(rtime0,&rtime1,&use2,&use1,&dtu,&drt,&tsize2,&tsize1,it+span);

   } /* end of time step */

if(modstor.intmem == 0)
   {
   sprintf(string,"%s/%s.%d.pv_e3d",tmpdir,name,ny1);
   unlink(string);
   sprintf(string,"%s/%s.%d.med_e3d",tmpdir,name,ny1);
   unlink(string);
   }

if(intface)     /* remove temporary interfacing files */
   rm_infc_files(&intfac);

if(dumpinfo.enable_flag == 1)   /* remove temporary output files */
   rm_dump_files(&outp,segmentId,&dumpinfo);

total_report(rtime0,&use2,tsize2);

mpi_final("PROGRAM emod3d-mpi IS FINISHED");
exit(0);
}
