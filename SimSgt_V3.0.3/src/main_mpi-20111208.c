/*******************************************************************************/
/*                                                                             */
/*         3-D Elastic Finite-difference Modeling Code                         */
/*                                                                             */
/*                    copyright (c) 1993                                       */
/*                     Robert W. Graves                                        */
/*                 Woodward-Clyde Consultants                                  */
/*                   566 El Dorado Street                                      */
/*                    Pasadena, CA 91101                                       */
/*                                                                             */
/*                    tel (818) 449-7650                                       */
/*                    fax (818) 449-3536                                       */
/*                                                                             */
/*    Permission to copy all or part of this work is granted,                  */
/*    provided that the copies are not made or distributed                     */
/*    for resale, and that the copyright notice and this                       */
/*    notice are retained.                                                     */
/*                                                                             */
/*    NOTE:  This work is provided on an "as is" basis.                        */
/*           The author provides no warranty whatsoever,                       */
/*           expressed or implied, regarding the work,                         */
/*           including warranties with respect to its                          */
/*           fitness for any particular purpose.  The code                     */
/*           is provided on the condition that it will not                     */
/*           be transfered to third parties without the                        */
/*           explicit consent of the author.                                   */
/*           In addition, the code has not been                                */
/*           tested extensively, therefore there are no                        */
/*           guarantees about its performance.  I am currently                 */
/*           in the process of testing and benchmarking the                    */
/*           code and would appreciate any comments and/or                     */
/*           feedback from those who use it.                                   */
/*                                                                             */
/*                             -RWG 05/05/93                                   */
/*                                                                             */
/*******************************************************************************/
/*******************************************************************************/
/*                                                                             */
/*			version 3.0: 2011-11-22                                */
/*                                                                             */
/*******************************************************************************/
/*******************************************************************************/
/*                                                                             */
/*    A note about units.                                                      */
/*                                                                             */
/*      - Global input units are:                                              */
/*                                                                             */
/*        # wavefield velocity is-                                             */
/*           vx,vy,vz       => cm/s                                            */
/*                                                                             */
/*        # seismic velocity and density are-                                  */
/*           Vp,Vs     => km/s                                                 */
/*           density   => gm/cm^3                                              */
/*                                                                             */
/*        # Lame parameters and stress are-                                    */
/*           Tnn,Tij,lambda,mu => km*km*gm/(s*s*cm*cm*cm)                      */
/*                             = 1e+09 Nt/(m^2)                                */
/*                             = 1e+09 Pa                                      */
/*                             = 1 GPa                                         */
/*                                                                             */
/*        # dimension and time measures-                                       */
/*           h (distance, length) => km                                        */
/*           dt (time)            => s                                         */
/*                                                                             */
/*        # sources are CGS-                                                   */
/*           moment    => dyne-cm                                              */
/*           force     => dyne                                                 */
/*           slip      => cm                                                   */
/*           slip-rate => cm/s                                                 */
/*                                                                             */
/*    If the input variables are expressed in these units, then all            */
/*    internal calculations can be done without any additional unit            */
/*    conversions.                                                             */
/*                                                                             */
/*       - For output, units are:                                              */
/*                                                                             */
/*        # wavefield velocities (single station & time slices) are-           */
/*           vx,vy,vz       => cm/s                                            */
/*                                                                             */
/*        # strain Greens tensor (SGT, reciprocal GF) are output in            */
/*          units of strain per unit force-                                    */
/*           gxx,gyy,gzz,gxy,gxz,gyz => 1/(dyne)                               */
/*                                                                             */
/*          thus when convolved with a moment-rate function (dyne*cm/s)        */
/*          the result is velocity in (cm/s)-                                  */
/*           dM/dt*gij => (dyne*cm/s) * 1/(dyne) => cm/s                       */
/*                                                                             */
/*******************************************************************************/
/*****************    GENERAL INFORMATION    ***********************************/
/*                                                                             */
/*    Axes orientation is as follows:                                          */
/*                                                                             */
/*         x-axis -> horizontal plane, positive to the east                    */
/*         y-axis -> horizontal plane, positive to the south                   */
/*         z-axis -> vertical plane, positive downward                         */
/*                                                                             */
/*    In order to run this code on a machine which has                         */
/*    limited core memory availble (ie., less than is needed                   */
/*    to store the entire wave field and media parameters),                    */
/*    the program uses disk memory to store the majority of                    */
/*    the model as the wave field is being updated.                            */
/*    Field and media parameters are read in from disk in                      */
/*    slices of xz planes.  This optimizes the number of                       */
/*    time steps that are calculated per IO operation.                         */
/*                                                                             */
/*    The variables which control these operations are:                        */
/*         span    -> number of time steps per IO operation                    */
/*         spanlen =  4*span + 2                                               */
/*                 -> number of 2D slices of model needed in                   */
/*                    in core at one time.                                     */
/*         memlen  =  spanlen*nx*nz*size_float                                 */
/*                 -> memory needed for "spanlen" model slices                 */
/*                    for one variable (18 variables needed                    */
/*                    for elastic case with Q).                                */
/*                                                                             */
/*    Field varibles are stored in pvf array as follows:                       */
/*                                                                             */
/*         pvf          -> vx                                                  */
/*         pvf+   nx*nz -> vy                                                  */
/*         pvf+ 2*nx*nz -> vz                                                  */
/*         pvf+ 3*nx*nz -> txx                                                 */
/*         pvf+ 4*nx*nz -> tyy                                                 */
/*         pvf+ 5*nx*nz -> tzz                                                 */
/*         pvf+ 6*nx*nz -> txy                                                 */
/*         pvf+ 7*nx*nz -> txz                                                 */
/*         pvf+ 8*nx*nz -> tyz                                                 */
/*         pvf+ 9*nx*nz -> cxx                                                 */
/*         pvf+10*nx*nz -> cyy                                                 */
/*         pvf+11*nx*nz -> czz                                                 */
/*         pvf+12*nx*nz -> cxy                                                 */
/*         pvf+13*nx*nz -> cxz                                                 */
/*         pvf+14*nx*nz -> cyz                                                 */
/*                                                                             */
/*    The variables cij above are the memory variables for anelasticity.       */
/*                                                                             */
/*    Media varibles are stored in medf array as follows:                      */
/*                                                                             */
/*         medf          -> lambda + 2*mu                                      */
/*         medf+   nx*nz -> lambda                                             */
/*         medf+ 2*nx*nz -> muxy                                               */
/*         medf+ 3*nx*nz -> muxz                                               */
/*         medf+ 4*nx*nz -> muyz                                               */
/*         medf+ 5*nx*nz -> bouyx                                              */
/*         medf+ 6*nx*nz -> bouyy                                              */
/*         medf+ 7*nx*nz -> bouyz                                              */
/*         medf+ 8*nx*nz -> ak                                                 */
/*         medf+ 9*nx*nz -> pk                                                 */
/*         medf+10*nx*nz -> sk                                                 */
/*         medf+11*nx*nz -> lambda (raw, never altered from input)             */
/*         medf+12*nx*nz -> mu     (raw, never altered from input)             */
/*                                                                             */
/*    The variables ak, pk, and sk above are required to implement the         */
/*    memory variable formulation for anelasticity.                            */
/*                                                                             */
/*    The values for rigidity (mu) and bouyancy (1/density)                    */
/*    have been averaged to avoid instability when media                       */
/*    contrasts are large (see Randall et al., Nov., 1991,                     */
/*    Geophysics).  The averages are precomputed and stored                    */
/*    in the array medf in order to make the calculations more                 */
/*    efficient.  Refering to the table above, the media                       */
/*    parameters are averaged in the following manner:                         */
/*                                                                             */
/*       muxy -> avg. over xy plane for txy computation                        */
/*       muxz -> avg. over xz plane for txz computation                        */
/*       muyz -> avg. over yz plane for tyz computation                        */
/*       bouyx -> avg. over x plane for vx computation                         */
/*       bouyy -> avg. over y plane for vy computation                         */
/*       bouyz -> avg. over z plane for vz computation                         */
/*                                                                             */
/*    An arithmetic average is used for bouyancy and a                         */
/*    harmonic average is used for rigidity.  All averages                     */
/*    are computed in the forward sense.  That is                              */
/*                                                                             */
/*       bouyx[i] = 2/(rho[i] + rho[i+1])                                      */
/*       .                                                                     */
/*       .                                                                     */
/*       .                                                                     */
/*       muxy[i] =                                                             */
/*          4/(1/mu1[i] + 1/mu1[i+1] + 1/mu2[i] + 1/mu2[i+1])                  */
/*       .                                                                     */
/*       .                                                                     */
/*       .                                                                     */
/*                                                                             */
/*    No averaging is performed at the maximal edges (ix=nx,                   */
/*    iy=ny,iz=nz) of the model grid.                                          */
/*                                                                             */
/*    Absorbing boundary conditions are applied to velocities                  */
/*    only, the outer ring of stress values are not updated.                   */
/*    To implement the absorbing boundaries, the velocity                      */
/*    field at the previous time step is required.  This is                    */
/*    stored in the array pvbuf.                                               */
/*                                                                             */
/*    For the boundaries at iy=0 and iy=ny-1, the planes of                    */
/*    vx, vy and vz at iy=1 and iy=ny-2, respectively, are                     */
/*    stored prior to updating the velocity field using the                    */
/*    function storev().  The absorbing boundary conditions                    */
/*    are applied using the function abs_ybnd() which is                       */
/*    called from the function tstepvbnd().  For these planes,                 */
/*    the velocity fields are stored in the array pvbuf in the                 */
/*    following manner:                                                        */
/*                                                                             */
/*         pvbuf         -> vx                                                 */
/*         pvbuf+nx*nz   -> vy                                                 */
/*         pvbuf+2*nx*nz -> vz                                                 */
/*                                                                             */
/*    For the boundaries of x and z, the boundary condition                    */
/*    for vx, vy and vz are applied using the function                         */
/*    absorb() which is called from tstepv().  The inner ring                  */
/*    of values are stored using the function store_iring()                    */
/*    in the following manner:                                                 */
/*                                                                             */
/*         pvbuf+3*nx*nz           -> vx                                       */
/*         pvbuf+3*nx*nz+2*(nx+nz) -> vy                                       */
/*         pvbuf+3*nx*nz+4*(nx+nz) -> vz                                       */
/*                                                                             */
/*    Velocity field is updated first, then stress field.                      */
/*                                                                             */
/*******************************************************************************/
/************************    NOTES    ******************************************/
/*                                                                             */
/*  NOTE: (03/16/93 -RWG)                                                      */
/*    The free surface boundary condition has not been tested                  */
/*    completely.                                                              */
/*                                                                             */
/*  NOTE: (03/08/97 -RWG)                                                      */
/*    The free surface boundary condition has been tested in many ways,        */
/*    and seems to perform very well.                                          */
/*                                                                             */
/*******************************************************************************/
/*******************    MAJOR REVISIONS    *************************************/
/*                                                                             */
/*  02/12/93 -RWG                                                              */
/*    Changed calls to ord2(), ord4(), etc. so that stride=1                   */
/*    in all cases.  I discovered that when nx and/or nz are                   */
/*    large and stride=nx (ie., x derivatives in previous                      */
/*    coding), the computation time can be slowed by up to a                   */
/*    factor of 2.  I think this is caused by memory swapping                  */
/*    when the arrays get large.  This may be a hardware                       */
/*    limitation or it may be set by the system software.  I                   */
/*    have only noticed the problem on the Sparc2, it does not                 */
/*    seem to occur on the Sun4.                                               */
/*                                                                             */
/*  04/06/93 -RWG                                                              */
/*    Boundaries along ix=nx-1, iy=ny-1 and iz=nz-1 are moved                  */
/*    inward 1/2 grid step to stabilize the absorbing                          */
/*    boundaries.  Thus, the interior solution is not                          */
/*    calculated for vx along ix=nx-2, for vy along iy=ny-2                    */
/*    or for vz along iz=nz-2.  These values are obtained                      */
/*    using the absorbing boundary condition (first-order                      */
/*    Clayton-Enquist).  The other variables affected are txx,                 */
/*    tyy and tzz.                                                             */
/*                                                                             */
/*  07/29/94 -RWG                                                              */
/*    Set up memory configuration, so that model can either be                 */
/*    stored internally in core memory or externally on disk,                  */
/*    with user controlled IO.  Information for taking model                   */
/*    slices in or out of memory is contained in the following                 */
/*    structure:                                                               */
/*                                                                             */
/*            struct modelstorage                                              */
/*               {                                                             */
/*               size_t maxmem;                                                */
/*               int intmem;                                                   */
/*               size_t ln_pvslice;                                            */
/*               size_t sz_pvslice;                                            */
/*               int pv_fdr;                                                   */
/*               int pv_fdw;                                                   */
/*               float *pv_buf;                                                */
/*               float *pv_bufptr;                                             */
/*               size_t ln_medslice;                                           */
/*               size_t sz_medslice;                                           */
/*               int med_fdr;                                                  */
/*               int med_fdw;                                                  */
/*               float *med_buf;                                               */
/*               float *med_bufptr;                                            */
/*               };                                                            */
/*                                                                             */
/*    The variables 'maxmem' (maximum core memory limit, in                    */
/*    bytes) and 'intmem' (internal memory/external memory                     */
/*    flag) can be input as as getpar parameters.  Default                     */
/*    values are:                                                              */
/*                                                                             */
/*            maxmem = MEM_LIMIT (defined in defs.h)                           */
/*            intmem = 0 (use external storage)                                */
/*                                                                             */
/*    By setting 'intmem=1', the model will be stored entirely                 */
/*    in core memory.  This will significantly reduce system                   */
/*    CPU time, but does require a large amount of available                   */
/*    memory.  The memory limit imposed by 'maxmem' is a soft                  */
/*    limit, active only within this code.  I typically set                    */
/*    the value of 'maxmem' to be one-half of the total                        */
/*    internal memory (excluding swap space) available on the                  */
/*    machine that the code is running.                                        */
/*                                                                             */
/*    BE CAREFUL!!  When the code is run with the parameter                    */
/*    'intmem=1', and other processes are accessing memory                     */
/*    resources, then there exists a much greater potential                    */
/*    for memory swapping to occur.  When this happens, the                    */
/*    run-time will increase dramatically.  Unless you can be                  */
/*    sure that this will not occur, then it is better to run                  */
/*    the code with 'intmem=0'.  Furthermore, system overhead                  */
/*    for user controlled IO when 'intmem=0' is usually on the                 */
/*    order of 1-2 percent of the total CPU time.  Thus there                  */
/*    appears to be very little advantage in using only core                   */
/*    memory.                                                                  */
/*                                                                             */
/*  03/08/97 -RWG                                                              */
/*                                                                             */
/*   -Changed moment-tensor source insertion to stresses (ala                  */
/*    S. Day) instead of velocities (body forces as described                  */
/*    in my paper).  The results are virtually identical, but                  */
/*    the stress formulation is somewhat more compact.                         */
/*                                                                             */
/*   -Also added the option to specify input weights in the 'faultfile'        */
/*    in terms of relative slip, instead of relative moment.  This lets        */
/*    the scaling to moment to be done using the actual velocity model         */
/*    that the FD code is running.  This is important if the fault             */
/*    intersects a basin or encounters some lateral velocity variation.        */
/*    To envoke this option, set the parameter "relative_slip" in the          */
/*    getpar parameter file (eg., 'e3d.par') to a non-zero value, ie.,         */
/*                                                                             */
/*       relative_slip=1                                                       */
/*                                                                             */
/*    The default is still relative moments (relative_slip=0).                 */
/*    Remember, that if you use this option, the weights in the                */
/*    'faultfile' must also be given in terms of relative slip.  There         */
/*    is a corresponding parameter ("relative_slip") that can be passed        */
/*    to the codes 'lisa2mtgrid' and 'lisa2mtgrid_vrup' so that these          */
/*    codes will output the FD source weights in terms of relative slip.       */
/*    See me for more details.                                                 */
/*                                                                             */
/*  05/09/98 -RWG                                                              */
/*                                                                             */
/*   -To alleviate instability at absorbing boundaries when media              */
/*    has very strong lateral contrast at close proximity (eg., 1 grid         */
/*    point) exactly parallel to boundary, the input velocity model is         */
/*    modified by replicating the plane of grid points adjacent to the         */
/*    boundary.  Thus, the iy=2 plane is copied to iy=1 and iy=0.  A           */
/*    similar type of replication is done for x and z (excluding free          */
/*    surface).  The replication (or padding) width is controlled by the       */
/*    parameter NBND_PAD, and associated functions are copy_medslice()         */
/*    and xzpad_medslice().                                                    */
/*                                                                             */
/*  07/98 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.6 has the added option to output moment-tensor                 */
/*    components for use in calculating reciprocal GF's.  See mtheader         */
/*    info.                                                                    */
/*                                                                             */
/*  02/99 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.7 has modified the output routines for MT components to        */
/*    include information about the model origin (lat,lon), grid               */
/*    spacing, and rigidity at the output location.  See mtheader              */
/*    structure information.                                                   */
/*                                                                             */
/*  04/00 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.8 has modified IO routines for external model storage          */
/*    routines for large file compatibility.  Changes to the source code       */
/*    are trivial (see open() calls in iofunc.c.  Compiler options are         */
/*    more noticable (see LF_FLAGS in makefile) and are definitely             */
/*    system and OS version dependent.  For SUNs the large file option         */
/*    is only available on solaris 2.6 and later OS versions.                  */
/*                                                                             */
/*   -Version 1.8 has modified cpu time usage routines based on the            */
/*    system function times().  This was done to circumvent portability        */
/*    problems (still needs to be tested) and problems using procfs.h          */
/*    with large file (64-bit) configuration.                                  */
/*                                                                             */
/*  06/00 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.9 has option to center media variables either at the           */
/*    normal stress node ("media_boundary=0") or at the NULL unit cell         */
/*    node ("media_boundary=1").  The NULL unit cell node is the only          */
/*    corner of the unit cell which does not have any field variables          */
/*    defined, and it lies diagonally opposite of the normal stress            */
/*    node.  Based on the definition of the exact location of the media        */
/*    variables, different formulas are needed to calculate the                */
/*    effective media parameters that are used for the time updates of         */
/*    the various wave field variables (see Zahradnik et al, 1993).            */
/*    The effective media parameters are calculated in the function            */
/*    avgmedslice() in genmodel.c.                                             */
/*                                                                             */
/*    * media_boundary=0                                                       */
/*    With the media variables centered at the normal stress node, the         */
/*    media boundaries go through the velocity nodes.  This results in a       */
/*    0.5*h depth difference in the location of the media boundaries and       */
/*    the location of the free surface (see Graves, 1996).  Thus, the          */
/*    thinnest possible layer has depth of '0.5*h', and all other layers       */
/*    have a depth of '(k+0.5)*h' [k=integer].                                 */
/*                                                                             */
/*    * media_boundary=1                                                       */
/*    With the media variables centered at the NULL unit cell node, the        */
/*    media boundaries cross at the normal stress node.  This results in       */
/*    a co-location of the media boundaries and the location of the free       */
/*    surface.  In this case, the thinnest possible layer has depth of         */
/*    'h', and all other layers have a depth of 'k*h' [k=integer].             */
/*                                                                             */
/*   -Version 1.9 has cleaned-up the coding used to generate the               */
/*    effective media parameters and write these values to the temporary       */
/*    model files.                                                             */
/*                                                                             */
/*   -Version 1.9 has changed the input format specification for a             */
/*    finite fault in the file 'faultfile'.  The new format specifies          */
/*    the point source locations in terms of distances from the model          */
/*    origin.  The distances are floating point coordinates (x,y,z).           */
/*    Earlier versions specified the point source locations as grid            */
/*    indices (ix,iy,iz) relative to the model origin.  For the                */
/*    horizontal directions these are simply related by                        */
/*                                                                             */
/*                x = ix*h    y=iy*h                                           */
/*                                                                             */
/*    However, for the vertical direction the 1 grid downward shift for        */
/*    the free-surface implementation is NOT included in the new               */
/*    floating point faultfile coordinates.  This means that a point           */
/*    located exactly at the free surface in versions 1.8 and earlier          */
/*    is given by iz=1, whereas in version 1.9, this point is located          */
/*    at z=0.0.  The relation between these two specification is given         */
/*                                                                             */
/*                iz = (int)(z/h + 0.5) + 1                                    */
/*                                                                             */
/*    The additional "1" in the above equation provides the downward           */
/*    shift to account for the free-surface.                                   */
/*                                                                             */
/*    Note that the above format specification only applies to the             */
/*    coordinates specified in the 'faultfile'.  All other coordinate          */
/*    specifications in version 1.9 are given in terms of grid indices         */
/*    (ix,iy,iz), and these must include the 1 grid shift for the              */
/*    free-surface (just the same as earlier versions).  The reason for        */
/*    the format change in the 'faultfile' is to make the 'faultfile'          */
/*    specification compatible between uniform gird versions (v1.9) and        */
/*    non-uniform-grid version (v2.0) of emod3d.                               */
/*                                                                             */
/*  10/01 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.10 has implemented the coarse-grained memory variable          */
/*    formulation to model anelaticity.  See Day and Bradley (BSSA,            */
/*    June 2001) and my notes.                                                 */
/*                                                                             */
/*  05/05 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.11 has implemented the RGF output option with the MPI          */
/*    parallel code formulation.                                               */
/*                                                                             */
/*  09/05 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.12:                                                            */
/*                                                                             */
/*     # Cleaned-up indexing for media planes in init_media() and              */
/*       eff_media().  Originally Sidao added two planes to media array        */
/*       (one in begining and one at the end) because he thought they          */
/*       were needed to perform media averaging.  Thus, media array was        */
/*       dimensioned at ny+2 instead of ny.  The indexing was adjusted         */
/*       in eff_media() in Versions 1.11 and earlier.                          */
/*                                                                             */
/*       Now, media array is dimensioned at ny and all indexing is set         */
/*       locally to iy=0,...,ny-1 or gloabally to iy=ny1,...,ny2-1.            */
/*       Functions affected by this are init_media(), eff_media(), and         */
/*       slip2mom().                                                           */
/*                                                                             */
/*     # reformatted 'struct momenttensor'                                     */
/*                                                                             */
/*       * 'struct pntsrcs' (point sources) now has a pointer to               */
/*         'struct momenttensor'.  Memory (for nsource) is allocated in        */
/*         init_dc(). Should be backward compatible with previous              */
/*         versions.                                                           */
/*                                                                             */
/*  12/05 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.13:                                                            */
/*                                                                             */
/*     # Implemented restart capabilities.                                     */
/*                                                                             */
/*       * Ouput options are reduced to 1) point time histories (nseis)        */
/*         with all_in_one=1 the only option, 2) time slice files, and         */
/*         3) SGT output.  Initialization of output is now done in             */
/*         init_outp().                                                        */
/*                                                                             */
/*  05/10 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.15:                                                            */
/*                                                                             */
/*     # Implemented buffering of xy-time slice output                         */
/*                                                                             */
/*  03/11 -RWG                                                                 */
/*                                                                             */
/*   -Version 1.16:                                                            */
/*                                                                             */
/*     # Hardwired span=1, intmem=1                                            */
/*                                                                             */
/*     # Removed roll-in, interior, roll-out sequencing                        */
/*                                                                             */
/*                                                                             */
/*  11/2011 -RWG                                                               */
/*                                                                             */
/*   -Version 3.0:                                                             */
/*                                                                             */
/*    Version 3.0 adds the flexibility to divide the model along the           */
/*    X and Z directions, in addition to the Y direction.  This results        */
/*    in full 3D parallelization, and provides for more efficient use          */
/*    of memory and much greater scalability to large numbers of               */
/*    processors.                                                              */
/*                                                                             */
/*    Information on model distribution and local coordinates/indices          */
/*    is contained in the structure nodeinfo:                                  */
/*                                                                             */
/*      struct nodeinfo                                                        */
/*         {                                                                   */
/*         int nproc;                                                          */
/*         int nproc_x;                                                        */
/*         int nproc_y;                                                        */
/*         int nproc_z;                                                        */
/*         int min_nproc;                                                      */
/*         int segmentId;                                                      */
/*         int minusId_x;                                                      */
/*         int plusId_x;                                                       */
/*         int minusId_y;                                                      */
/*         int plusId_y;                                                       */
/*         int minusId_z;                                                      */
/*         int plusId_z;                                                       */
/*         int procId_x;                                                       */
/*         int procId_y;                                                       */
/*         int procId_z;                                                       */
/*         int globnx;                                                         */
/*         int nx1;                                                            */
/*         int nx2;                                                            */
/*         int ixminus;                                                        */
/*         int ixplus;                                                         */
/*         int loc_nx;                                                         */
/*         int globny;                                                         */
/*         int ny1;                                                            */
/*         int ny2;                                                            */
/*         int iyminus;                                                        */
/*         int iyplus;                                                         */
/*         int loc_ny;                                                         */
/*         int globnz;                                                         */
/*         int nz1;                                                            */
/*         int nz2;                                                            */
/*         int izminus;                                                        */
/*         int izplus;                                                         */
/*         int loc_nz;                                                         */
/*         char procname[512];                                                 */
/*         };                                                                  */
/*                                                                             */
/*    which is cast as 'struct nodeinfo ninfo' in main().  This generalizes    */
/*    previously used parameters like nodeType, iyleft, iyright, ny1, etc.     */
/*                                                                             */
/*    The parameters in ninfo are set using the functions:                     */
/*                                                                             */
/*      get_nproc(&ninfo);                                                     */
/*      get_nodetype_neighbor(&ninfo);                                         */
/*      get_n1n2(bflag,&ninfo);                                                */
/*                                                                             */
/*    New versions of many functions have been added to accomodate the code    */
/*    changes.  These are denoted with a 'P3' appended to the function name,   */
/*    e.g., get_srf_parP3() is the updated version of get_srf_par(). The       */
/*    basic operations in these functions is the same as in previous code      */
/*    versions, e.g., tstepvP3() updates the velocity components just as done  */
/*    previously by tstepv(). Also, in most cases, the modifications to        */
/*    the functions are rather minor.                                          */
/*                                                                             */
/*    Data exchange among nodes is now done using exchange_pvP3(), which       */
/*    cycles over all boundaries of each node (x0,xN,y0,yN,z0,zN) and          */
/*    sends/receives data as appropriate.  This is done using sequential       */
/*    calls of MPI_Sendrecv(), just as in previous versions.                   */
/*                                                                             */
/*    The outward changes to the code are not major.  There are some           */
/*    new input parameters (and some that have been removed), and all          */
/*    of the input and most of the output file formats are the same            */
/*    as in earlier (V1.16 and before) versions.                               */
/*    See below for more details.                                              */
/*                                                                             */
/*     # New input parameters:                                                 */
/*                                                                             */
/*       * Version 3.0 can be run with exactly the same input                  */
/*         parameters, files and options as previous versions. However,        */
/*         in this default parameterization, the resulting distribution        */
/*         of the model over the desired number of processors may not          */
/*         be very optimal or efficient.  In this mode, the code               */
/*         first tries to divide the model into equal sized cubes given        */
/*         the requested number of processors.  Cubic model subsets are        */
/*         desirable because they maximize the volume to surface area,         */
/*         which in turn minimizes the amount of inter-node                    */
/*         communications for a given calculation. If that doesn't work        */
/*         with the requested number of processors, the code then tries        */
/*         to set the number of processors along the X, Y and Z                */
/*         directions to be the same with the constraint that the total        */
/*         number be equal to the specified target.  This may or may           */
/*         not work in a very optimal way either.                              */
/*                                                                             */
/*         Alternatively, one can explicitly specify the number of             */
/*         processors along the X, Y and Z directions with the                 */
/*         parameters:                                                         */
/*                                                                             */
/*           nproc_x=<int>                                                     */
/*           nproc_y=<int>                                                     */
/*           nproc_z=<int>                                                     */
/*                                                                             */
/*         with the constraint np=nproc_x*nproc_y*nproc_z (np=requested        */
/*         number of processors).                                              */
/*                                                                             */
/*         I recommend using the code check_n1n2 to estimate/iterate to        */
/*         get the optimal parameterization.                                   */
/*                                                                             */
/*     # No longer needed/supported input parameters:                          */
/*                                                                             */
/*       * The following parameters are no longer needed                       */
/*                                                                             */
/*           span=<int>      (removed, implicitly span=1)                      */
/*           intmem=<int>    (hardwired to 1, all internal memory)             */
/*                                                                             */
/*       * The following parameters/options are no longer supported            */
/*                                                                             */
/*           psrc=1   (pressure source; now use Mxx, Myy and Mzz)              */
/*           expl=1   (explosion source; now use dblcpl=1, strike=-9999)       */
/*           eqsrc=1  (old style EQ source; now use dblcpl/pointmt)            */
/*           ffault=1 (old style finite-fault; now use ffault=2 -> SRF)        */
/*                                                                             */
/*     # Output file format changes:                                           */
/*                                                                             */
/*       Output file names are in the same form as previous versions-          */
/*                                                                             */
/*         <name>_<mode>-<procId>.e3d                                          */
/*                                                                             */
/*       However, now <procId> is a left-zero-padded integer of 5 character    */
/*       places (only 4 places were used in earlier versions).  This           */
/*       accomodates up to 99999 individual processors.                        */
/*                                                                             */
/*       * Point seismograms (nseis=1)                                         */
/*                                                                             */
/*         Point seismograms now output all 9 wave field components, 3         */
/*         velocities (just as previously) and 6 stresses.  The order of the   */
/*         components is Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz.                     */
/*         The code 'fdbin2wcc' has been modified to accept the input          */
/*         parameter                                                           */
/*                                                                             */
/*           nseis_comps=<int>   (default is 3 for backward compatibility)     */
/*                                                                             */
/*         When reading V3.0 output with fdbin2wcc, set nseis_comps=9.         */
/*                                                                             */
/*       * Time slices (ts_xy=1, ts_xz=1, or ts_yz=1)                          */
/*                                                                             */
/*         STILL NEED TO WORK ON THIS                                                */
/*                                                                             */
/*       * SGT output (sgtout=1)                                               */
/*                                                                             */
/*         STILL NEED TO WORK ON THIS                                          */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*                                                                             */
/*     # Testing status:                                                       */
/*                                                                             */
/*       Media specification-                                                  */
/*                                                                             */
/*       * model_style=0 (1D layered w/Q):      Nov 2011 (FdRuns_Test1)        */
/*       * model_style=1 (3D binary files w/Q): Nov 2011 (FdRuns_Test2)        */
/*       * elas_only=1:                                                        */
/*                                                                             */
/*       Sources-                                                              */
/*                                                                             */
/*       * pointmt: Nov 2011 (FdRuns_Test1)                                    */
/*       * bforce:  Nov 2011 (FdRuns_Test5)                                    */
/*       * dblcpl:                                                             */
/*       * SRF:     Nov 2011 (FdRuns_Test4a)                                   */
/*                                                                             */
/*       Outputs-                                                              */
/*                                                                             */
/*       * nseis=1 (point seismograms):       Nov 2011 (FdRuns_Test1, etc.)    */
/*       * sgtout=1 (strain greens tensors):                                   */
/*       * ts_xy=1 (xy time slice):                                            */
/*       * ts_xz=1 (xz time slice):                                            */
/*       * ts_yz=1 (yz time slice):                                            */
/*                                                                             */
/*******************************************************************************/

#include "include.h"

int main(int ac,char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */

int nx1, nx2, globnx;  /* start and end ix plane, total number of points */
int ny1, ny2, globny;  /* start and end iy plane, total number of points */
int nz1, nz2, globnz;  /* start and end iz plane, total number of points */

struct nodeinfo ninfo;
int nproc, pnlen;
int bflag = 1;

int izord2 = -1;
size_t blen, off;
int vmodel_swapb = 0;
float *snd1, *rcv1, *snd2, *rcv2;   /* send/receive pointers */
int sndtag, rcvtag;

/* Rest of variables */
struct fdcoefs fdcoefs;
FILE *fpr, *fpw, *fopfile();
float *pvfield, *pvptr[4];                    /* field varibles */
float *medfield, *medptr[4], *tmpmedfield;   /* media variables */
float *stfunc, *stfp;

/*
   pvbuf used as temporary buffer for absorbing boundaries
   and for exchange of wave field among nodes
*/
float *pvbuf;

float *tmp_ptr;
float vtmp, dt, h;
int nx, ny, nz, nt, nbuf, iy, it, i, ib, ybndflag;

int ix, iz, iv, ip;

int icnt;
int nbpad, print_tol;
int nt_p2;

float tzero;
int bfilt = 0;
float fhi = 0.0;
float flo = 1.0e+15;

float dep_val;
int intsrc = 0;
int nbounce = 1;

char version[16], vid[16];
int vnc;

struct modelstorage modstor;
struct media_input medinp;

int elas_only = 0;
int eflag = 0;

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
struct restartinfoP3 rstart_p3;

                /*  rlimit parameters  */

struct rlimit rlp;
int maxfiles;

                /*  model run parameters  */

struct runparams rpars;
struct runparamsP3 rpars_p3;

                /*  output parameters  */

struct outputfields outp;

struct tsoutputparams *tsoutptr;
struct seisoutputparams *soutptr;
struct sgtoutputparams *sgtoutptr;

char name[FILE_STR_LEN];
char logfile[FILE_STR_LEN];
char tmpdir[DIR_STR_LEN];
char tmpdirlist[DIR_STR_LEN];
char vmoddir[DIR_STR_LEN];
char logdir[DIR_STR_LEN];
char stfdir[DIR_STR_LEN];

/* source information */

int intface = 0;
struct interface intfac;

int plane = 0;
struct planesrc plsrc;

int isrc;
char slipout[256], xyzfrmt[8], tfrmt[8];
char stype[256], stf_file[512];
struct pntsrcs srcs;

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
ninfo.min_nproc = 1;

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

rpars_p3.geoproj = 1;
rpars_p3.xshift = -1.0e+15;
rpars_p3.yshift = -1.0e+15;

rpars_p3.modelrot = 0.0;  /* rotation of y-axis from south (clockwise positive) */
rpars_p3.modellat = 34.0;    /* latitude of model origin */
rpars_p3.modellon = -118.0;    /* longitude of model origin */

rpars_p3.xmom = 0.0;
rpars_p3.ymom = 0.0;
rpars_p3.zmom = 0.0;

/* memory defaults */

modstor.maxmem = MEM_LIMIT;
modstor.intmem = 1;

setpar(ac,av);
mstpar("version","s",version);		/* version number */

getpar("name","s",name);
getpar("logdir","s",logdir);

makedir(logdir);
sprintf(logfile,"%s/%s-%.5d.rlog",logdir,name,ninfo.segmentId);
freopen(logfile,"w",stderr);

if(strncmp(version,vid,vnc) != 0)
   {
   fprintf(stderr,"***** This is version %s of emod3d-mpi:\n",vid);
   fprintf(stderr,"      Input parameters specify version %s\n",version);
   fprintf(stderr,"      Check input, exiting...\n");
   mpi_exit(-1);
   }

fprintf(stderr,"***** This is version %s of emod3d-mpi *****\n\n",vid);

getpar("maxmem","d",&modstor.maxmem);
modstor.maxmem = modstor.maxmem*1000000; /* convert Mbytes -> bytes */

mstpar("nx","d",&nx);             /* model parameters */
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("h","f",&h);
mstpar("nt","d",&nt);
mstpar("dt","f",&dt);

getpar("nproc_x","d",&ninfo.nproc_x);
getpar("nproc_y","d",&ninfo.nproc_y);
getpar("nproc_z","d",&ninfo.nproc_z);
getpar("min_nproc","d",&ninfo.min_nproc);

/* MPI-RWG: compute info for distributed computation */

/*
   Set-up x,y,z indexing parameters:

   ninfo.globnx = global number of x grid points for entire model
   ninfo.globny = global number of y grid points for entire model
   ninfo.globnz = global number of z grid points for entire model

   nx = ninfo.loc_nx = local number of x grid points for this node (including overlaps)
   ny = ninfo.loc_ny = local number of y grid points for this node (including overlaps)
   nz = ninfo.loc_nz = local number of z grid points for this node (including overlaps)

   ninfo.nx1 = global index of first x grid location in this node (including overlap)
   ninfo.nx2 = global index of last x grid location in this node + 1 (including overlap)
   ninfo.ny1 = global index of first y grid location in this node (including overlap)
   ninfo.ny2 = global index of last y grid location in this node + 1 (including overlap)
   ninfo.nz1 = global index of first z grid location in this node (including overlap)
   ninfo.nz2 = global index of last z grid location in this node + 1 (including overlap)

   Thus, for example, nx = ninfo.nx2 - ninfo.nx1
   
   ninfo.ixminus  = global x index of first responsible point in this node
   ninfo.ixplus = global x index of last responsible point in this node
   ninfo.iyminus  = global y index of first responsible point in this node
   ninfo.iyplus = global y index of last responsible point in this node
   ninfo.izminus  = global z index of first responsible point in this node
   ninfo.izplus = global z index of last responsible point in this node

   'Responsible point' means this node is responsible for computing full
   time update at this location.

   In general:

      ninfo.ixminus = ninfo.nx1 + 2
      ninfo.ixplus = ninfo.nx2 - 3
      ninfo.iyminus = ninfo.ny1 + 2
      ninfo.iyplus = ninfo.ny2 - 3
      ninfo.izminus = ninfo.nz1 + 2
      ninfo.izplus = ninfo.nz2 - 3

   However, if a node is along a computational boundary of the model:
   
   if ninfo.minusId_x = -1, then ninfo.ixminus = 0
      ninfo.minusId_y = -1, then ninfo.iyminus = 0
      ninfo.minusId_z = -1, then ninfo.izminus = 0
   
   if ninfo.plusId_x = -1, then ninfo.ixplus = ninfo.globnx - 1
      ninfo.plusId_y = -1, then ninfo.iyplus = ninfo.globny - 1
      ninfo.plusId_z = -1, then ninfo.izplus = ninfo.globnz - 1

   blfag=0 -> balance only on interior nodes
   blfag=1 -> balance across all nodes

   generally bflag=1 is most efficient

*/

ninfo.globnx = nx;
ninfo.globny = ny;
ninfo.globnz = nz;

bflag = 1;

get_nproc(&ninfo);				/* determines processor distribution */
get_nodetype_neighbor(&ninfo);			/* determines processor neighbors */
get_n1n2(bflag,&ninfo);				/* find n1, n2, loc_nn, minus, plus indices */

nx = ninfo.loc_nx;
ny = ninfo.loc_ny;
nz = ninfo.loc_nz;

if(nx < 6 || ny < 6 || nz < 6)
   {
   fprintf(stderr,"**** local nx= %d ny= %d nz= %d, all must be >= 6\n",nx,ny,nz);
   fprintf(stderr,"**** Need to reduce number of processors, exiting...\n");
   fflush(stderr);
   mpi_exit(-1);
   }

fprintf(stderr,"**** Running on node= %s\n",ninfo.procname);
fprintf(stderr,"          Total nproc= %-5d     nproc_x= %-5d nproc_y= %-5d nproc_z= %-5d\n",ninfo.nproc,ninfo.nproc_x,ninfo.nproc_y,ninfo.nproc_z);
fprintf(stderr,"               nodeId= %-5d\n",ninfo.segmentId);
fprintf(stderr,"            minusId_x= %-5d     plusId_x= %-5d\n",ninfo.minusId_x,ninfo.plusId_x);
fprintf(stderr,"            minusId_y= %-5d     plusId_y= %-5d\n",ninfo.minusId_y,ninfo.plusId_y);
fprintf(stderr,"            minusId_z= %-5d     plusId_z= %-5d\n\n",ninfo.minusId_z,ninfo.plusId_z);
fflush(stderr);

fprintf(stderr,"**** Distribution of model subsets:\n");
fprintf(stderr,"             local_nx= %-5d     (nx1= %-5d nx2= %-5d global_nx= %-5d)\n",ninfo.loc_nx,ninfo.nx1,ninfo.nx2,ninfo.globnx);
fprintf(stderr,"             local_ny= %-5d     (ny1= %-5d ny2= %-5d global_ny= %-5d)\n",ninfo.loc_ny,ninfo.ny1,ninfo.ny2,ninfo.globny);
fprintf(stderr,"             local_nz= %-5d     (nz1= %-5d nz2= %-5d global_nz= %-5d)\n",ninfo.loc_nz,ninfo.nz1,ninfo.nz2,ninfo.globnz);
fprintf(stderr,"                    h= %f\n",h);
fprintf(stderr,"                    nt= %d\n",nt);
fprintf(stderr,"                    dt= %f\n\n",dt);
fflush(stderr);

blen = set_model_storageP3(&modstor,nx,ny,nz);

pvfield = modstor.pv_buf;
pvbuf = modstor.pv_buf + N_WAVE_VARS*nx*ny*nz;
medfield = modstor.med_buf;

fprintf(stderr,"**** Memory stats for model (approximate):\n");
fprintf(stderr,"          total for model= %8.1f Mb\n\n",(float)(blen/1.0e+06));
fflush(stderr);

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

/*
   Source types 'psrc', 'expl', 'eqsrc' and 'ffault=1' (old finite-fault format) are
   not supported starting with V3.0 (20111116).  Explosion source can specified
   with dblcpl=1 and strike=-9999.
*/

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

if(plane || intface)
   {
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

   get_dc_parP3(&srcs,&tzero,&rpars_p3);
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

   get_pointmt_parP3(&srcs,&tzero,&rpars_p3);
   }
else if(srcs.bforce)
   {
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.pointmt = 0;
   srcs.dblcpl = 0;
   srcs.ffault = 0;
   srcs.relative_slip = srcs.absolute_slip = 0;

   get_bforce_parP3(&srcs,&tzero,&rpars_p3);
   }
else if(srcs.ffault == 1)
   {
   fprintf(stderr,"**** ffault=1 no longer supported, exiting...\n");
   fflush(stderr);
   mpi_exit(-1);
/*
   srcs.eqsrc = 0;
   srcs.expl = 0;
   srcs.psrc = 0;
   srcs.bforce = 0;
   srcs.pointmt = 0;
   srcs.ffault = 0;
   srcs.dblcpl = 1;

   getpar("slipout","s",slipout);
   getpar("xyzfrmt","s",xyzfrmt);
   getpar("tfrmt","s",tfrmt);

   get_ff_par(&srcs,&dt,&h,freesurf,&tzero);
*/
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

getpar("modelrot","f",&rpars_p3.modelrot);
getpar("modellat","f",&rpars_p3.modellat);
getpar("modellon","f",&rpars_p3.modellon);
getpar("xshift","f",&rpars_p3.xshift);
getpar("yshift","f",&rpars_p3.yshift);
getpar("geoproj","d",&rpars_p3.geoproj);

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

rstart_p3.enable_flag = 0;
rstart_p3.read_flag = 0;
rstart_p3.itinc = -999999;
rstart_p3.it = -999999;
sprintf(rstart_p3.dir,".");

getpar("enable_restart","d",&rstart_p3.enable_flag);
getpar("read_restart","d",&rstart_p3.read_flag);

if(rstart_p3.enable_flag == 1 || rstart_p3.read_flag == 1)
   {
   getpar("restartdir","s",rstart_p3.dir);
   makedir(rstart_p3.dir);
   }

if(rstart_p3.enable_flag == 1)
   getpar("restart_itinc","d",&rstart_p3.itinc);

if(rstart_p3.read_flag == 1)
   {
   getpar("restartname","s",rstart_p3.name);

   sprintf(rstart_p3.readfile,"%s/%s_rst-%.5d.e3d",rstart_p3.dir,rstart_p3.name,ninfo.segmentId);
   rstart_p3.fdr = opfile(rstart_p3.readfile);
   }

/* dump output parameters */

dumpinfo.enable_flag = 0;
dumpinfo.itinc = -999999;

getpar("enable_output_dump","d",&dumpinfo.enable_flag);

strcpy(dumpinfo.main_dir,outp.seisdir); /* intialize with something just in case */
strcpy(dumpinfo.name,name);

if(dumpinfo.enable_flag == 1)
   {
   if(rstart_p3.enable_flag == 1)
      strcpy(dumpinfo.local_restartdir,rstart_p3.dir);
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

copy_nodeinfo(&(rpars_p3.ni),&ninfo);
rpars_p3.nx = nx;
rpars_p3.ny = ny;
rpars_p3.nz = nz;
rpars_p3.nt = nt;
rpars_p3.freesurf = freesurf;
rpars_p3.h = h;
rpars_p3.dt = dt;
rpars_p3.tdelay = tdelay;

init_modelcP3(&rpars_p3);

fprintf(stderr,"**** Initializing media values:\n");
fflush(stderr);

/* use paired MPI_Barrier() calls so nodes don't all hit the disk at same time */
/*
for(ib=0;ib<ninfo.segmentId;ib++)
   MPI_Barrier(MPI_COMM_WORLD);
*/

init_mediaP3(&medinp,medfield,&rpars_p3,&qval,vmodel_swapb);

/* compliment to above MPI_Barrier() call */
/*
for(ib=ninfo.segmentId;ib<ninfo.nproc;ib++)
   MPI_Barrier(MPI_COMM_WORLD);
*/

fprintf(stderr,"**** Initializing wave field values:\n");
fflush(stderr);

init_field_val(1.0e-20,pvfield,ny*N_WAVE_VARS*nx*nz);

if(plane || intface)
   {
   if(plane == 1) /* different from plane=2, see add_plane() */
      {
      plsrc.srcl2m = medfield[nx*nz - 1];
      plsrc.srclam = medfield[2*nx*nz - 1];
      plsrc.srcmu = 1.0/medfield[5*nx*nz - 1];

      plsrc.svel = sqrt(plsrc.srcmu/medfield[8*nx*nz-1]);
      plsrc.pvel = sqrt(plsrc.srcl2m/medfield[8*nx*nz-1]);

      for(iy=0;iy<ny;iy++)
         {
         pvptr[0] = pvfield + iy*N_WAVE_VARS*nx*nz;
         init_field(pvptr[0],N_WAVE_VARS*nx*nz);

/* 04/08/05 RWG: Use globny and glob_iy=ny1+iy for plane wave specification */
/* 04/01/2011 Need to fix init_plane_srcP3() */
/* XXXXXX STILL NEED TO TEST */
         init_plane_src(pvptr[0],nx,globny,nz,ny1+iy,&h,&tzero,&dt,&plsrc,nbounce);
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
      }
   }
else
   {
   intfac.go = 0;

   if(srcs.dblcpl)
      {
      fprintf(stderr,"**** Initializing double-couple sources:\n");

      init_dcP3(&srcs,&rpars_p3);
      intsrc = 0;

      fprintf(stderr,"\n");
      }

   if(srcs.pointmt)
      {
      fprintf(stderr,"**** Initializing point moment-tensor sources:\n");

      init_pointmtP3(&srcs,&rpars_p3);
      intsrc = 0;

      fprintf(stderr,"\n");
      }

   if(srcs.bforce)
      {
      fprintf(stderr,"**** Initializing body-force sources:\n");

      init_bforceP3(&srcs,&rpars_p3);
      intsrc = 0;

      fprintf(stderr,"\n");
      }

   if(srcs.ffault == 2)
      {
      fprintf(stderr,"**** Initializing finite-fault moment-tensor sources:\n");
      fflush(stderr);

      if(slipout[0] != '\0')
         sprintf(slipout,"%s-%.5d",slipout,ninfo.segmentId);

      get_srf_parP3(&srcs,&rpars_p3,bfilt,&flo,&fhi,&tdelay);
      init_srfP3(&srcs,medfield,slipout,xyzfrmt,tfrmt,&rpars_p3);

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
         getsource(stfp,&dt,nt,&srcs.rtime[i],nbounce,intsrc,stype,name,bfilt,&flo,&fhi,i,ninfo.segmentId,&tdelay,stfdir,stf_file);
         }
      }

   fprintf(stderr,"\n");
   }

if(rstart_p3.enable_flag == 1 || rstart_p3.read_flag == 1)
   fprintf(stderr,"     Restart mode invoked:\n\n");

if(rstart_p3.read_flag == 1)
   {
   reed(rstart_p3.fdr,&rstart_p3.it,sizeof(int));

   fprintf(stderr,"     -reading wavefield at time step= %d from file:\n",rstart_p3.it);
   fprintf(stderr,"      %s\n\n",rstart_p3.readfile);
   fflush(stderr);

   reed(rstart_p3.fdr,&rstart_p3.rpars_p3,sizeof(struct runparamsP3));
   check_rparsP3(rstart_p3.readfile,&rstart_p3.rpars_p3,&rpars_p3);
   }

if(rstart_p3.enable_flag == 1)
   {
   sprintf(rstart_p3.dumpfile,"NOT_STARTED_YET");

   if(rstart_p3.itinc < 1)
      rstart_p3.itinc = 100;

   fprintf(stderr,"     -restart wavefield output at every %d time steps\n\n",rstart_p3.itinc);
   fflush(stderr);
   }

if(dumpinfo.enable_flag == 1)
   {
   if(dumpinfo.itinc < 1)
      dumpinfo.itinc = 1000;
   }

fprintf(stderr,"**** Computing effective media parameters:\n");
fflush(stderr);

eflag = eff_mediaP3(&medinp,medfield,&rpars_p3);

if(elas_only == 1)
   fprintf(stderr,"\n**** 'Elastic Only' mode invoked, no attenuation operators will be used.\n\n");
else
   {
   fprintf(stderr,"**** Computing memory variable coefficients:\n");
   mvar_coefsP3(medfield,&rpars_p3,&vfzero,&qfzero,&fmin,&fmax,eflag);
   /*
   smv_coefs(&modstor,medfield,nx,ny,nz,&dt,nt,freesurf,&vfzero,&qfzero,&fmin,&fmax,eflag);
   */
   }

fprintf(stderr,"**** Computing FD coefficients (checking vmin & vmax):\n");
set_vminvmax_reedmodP3(medfield,&fdcoefs,nx,ny,nz,&vsmin);

/* MPI-RWG: find global max(izord2), vmin & vmax */

vtmp = fdcoefs.vmax;
mpi_global_val(&vtmp,&fdcoefs.vmax,"Vmax",1,MPI_FLOAT,MPI_MAX);

vtmp = fdcoefs.vmin;
mpi_global_val(&vtmp,&fdcoefs.vmin,"Vmin",1,MPI_FLOAT,MPI_MIN);

set_izord2P3(medfield,&fdcoefs,nx,ny,nz,&vsmin);

/* adjust to global coordinates for checking across all nodes */
fdcoefs.izord2 = fdcoefs.izord2 + ninfo.nz1;

if(izord2 < 0)
   izord2 = fdcoefs.izord2;

mpi_global_val(&izord2,&fdcoefs.izord2,"izord2",1,MPI_INT,MPI_MAX);

setcoefs(order,&fdcoefs,&dt,&h,&vsmin);
fflush(stderr);

/* now adjust back to local coordinates */
fdcoefs.izord2 = fdcoefs.izord2 - ninfo.nz1;

if(elas_only == 0 || medinp.dampwidth > 0)
   {
   rfcnt = 0;
   rfspan = 6;

   if(rfinc <= 0)
      rfinc = 1;

   randf = (float *) check_malloc (rfspan*nx*nz*sizeof(float));
   init_field_val(1.0e-20,randf,rfspan*nx*nz);
   }

              /* initialize output files */

fprintf(stderr,"**** Initializing output files:\n\n");
fflush(stderr);

/*
eventually need to remove old rpars from init_outpP3
*/
init_outpP3(&outp,&rpars_p3,name,&srcs,&rstart_p3,&dumpinfo,&rpars);

/* make sure everyone's at the same place */
MPI_Barrier(MPI_COMM_WORLD);

         /* set field pointers */

if(rstart_p3.read_flag == 1)
   {
   /* use paired MPI_Barrier() calls so nodes don't all hit the disk at same time */
/*
   for(ib=0;ib<ninfo.segmentId;ib++)
      MPI_Barrier(MPI_COMM_WORLD);
*/

   get_restartP3(&rstart_p3,pvfield,nx,ny,nz);

   /* compliment to above MPI_Barrier() call */
/*
   for(ib=ninfo.segmentId;ib<ninfo.nproc;ib++)
      MPI_Barrier(MPI_COMM_WORLD);
*/

/* if randfield is used, roll through rfpointers to ensure exact replication */
   if((elas_only == 0 || medinp.dampwidth > 0))
      {
      for(it=0;it<=rstart_p3.it;it++)
         {
         if((it+1)%(rfinc) == 0)
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

for(it=0;it<nt;it++)
   {
   if(it > rstart_p3.it)
      {

/*
for(iy=0;iy<ny;iy++)
   {
   for(iv=0;iv<9;iv++)
      {
      for(iz=0;iz<nz;iz++)
         {
	 for(ix=0;ix<nx;ix++)
	    {
	    ip = ix + iz*nx + iv*nx*nz + iy*9*nx*nz;
            if(isnan(pvfield[ip]) != 0)
	       {
	       fprintf(stderr,"BEFORE_V it= %d ix= %d iy= %d iz= %d iv= %d pv= %e\n",it,ix,iy,iz,iv,pvfield[ip]);
	       fflush(stderr);
	       }
	    }
         }
      }
   }
*/

      for(iy=2;iy<ny-1;iy++)     /* do velocity update */
         {
         for(icnt=0;icnt<4;icnt++)  /* set pointers */
            {
            pvptr[icnt] = pvfield + ((iy-2) + icnt)*N_WAVE_VARS*nx*nz;
            medptr[icnt] = medfield + ((iy-2) + icnt)*N_MED_VARS*nx*nz;
            }

         for(isrc=0;isrc<srcs.nsource;isrc++)
            {
            if(iy == srcs.iy[isrc])          /* add in body force */
               add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,it,&srcs,FSRC);
            }

         if(intfac.go)     /* set interfacing parameters */
            {
            intfac.iy = iy + 1;
            intfac.it = it;
            }

         ybndflag = 0;
         if(iy == 2)
            {
            ybndflag = SKIP_PLANE_Y1;
            if(ninfo.minusId_y == -1)
               ybndflag = YZERO;
            }
         if(iy == ny-2)
            {
            ybndflag = SKIP_PLANE_Y2;
            if(ninfo.plusId_y == -1)
               ybndflag = YN;
            }

         tstepvP3(nx,nz,pvptr,medptr,pvbuf,&fdcoefs,freesurf,ybndflag,&intfac,&ninfo);

         if(plane == 2)
            add_plane(stfunc,pvptr,nx,nz,it,&plsrc);

/*
         wroutP3(&outp,pvptr[1],medptr[1],iy-1,it,&rpars_p3);
         if(iy == 2 && ninfo.minusId_y == -1)
            wroutP3(&outp,pvptr[0],medptr[0],iy-2,it,&rpars_p3);
         if(iy == ny-2 && ninfo.plusId_y == -1)
            {
            wroutP3(&outp,pvptr[2],medptr[2],iy,it,&rpars_p3);
            wroutP3(&outp,pvptr[3],medptr[3],iy+1,it,&rpars_p3);
            }
*/
         }

/*
for(iy=0;iy<ny;iy++)
   {
   for(iv=0;iv<9;iv++)
      {
      for(iz=0;iz<nz;iz++)
         {
	 for(ix=0;ix<nx;ix++)
	    {
	    ip = ix + iz*nx + iv*nx*nz + iy*9*nx*nz;
            if(isnan(pvfield[ip]) != 0)
	       {
	       fprintf(stderr,"AFTER_V it= %d ix= %d iy= %d iz= %d iv= %d pv= %e\n",it,ix,iy,iz,iv,pvfield[ip]);
	       fflush(stderr);
	       }
	    }
         }
      }
   }
*/

   /* exchange velocities */

      rt1 = times(&tu1);
      blen = exchange_pvP3(pvfield,pvbuf,&ninfo,it,VEL);
      dtimes(tu1,rt1,&dtu,&drt,&tsize2,blen);

/*
for(iy=0;iy<ny;iy++)
   {
   for(iv=0;iv<9;iv++)
      {
      for(iz=0;iz<nz;iz++)
         {
	 for(ix=0;ix<nx;ix++)
	    {
	    ip = ix + iz*nx + iv*nx*nz + iy*9*nx*nz;
            if(isnan(pvfield[ip]) != 0)
	       {
	       fprintf(stderr,"BEFORE_T it= %d ix= %d iy= %d iz= %d iv= %d pv= %e\n",it,ix,iy,iz,iv,pvfield[ip]);
	       fflush(stderr);
	       }
	    }
         }
      }
   }
*/

      for(iy=2;iy<ny-1;iy++)     /* do stress update */
         {
         for(icnt=0;icnt<4;icnt++)  /* set pointers */
            {
            pvptr[icnt] = pvfield + ((iy-2) + icnt)*N_WAVE_VARS*nx*nz;
            medptr[icnt] = medfield + ((iy-2) + icnt)*N_MED_VARS*nx*nz;
            }

         for(isrc=0;isrc<srcs.nsource;isrc++)
            {
            if(iy == srcs.iy[isrc])          /* add in pressure source */
               add_src(stfunc,nt,&dt,pvptr,medptr,isrc,nx,nz,it,&srcs,PSRC);
            }

         ybndflag = 0;
         if(iy == 2)
            {
            ybndflag = SKIP_PLANE_Y1;
            if(ninfo.minusId_y == -1)
               ybndflag = YZERO;
            }
         if(iy == ny-2)
            {
            ybndflag = SKIP_PLANE_Y2;
            if(ninfo.plusId_y == -1)
               ybndflag = YN;
            }

         tsteppP3(nx,nz,pvptr,medptr,&fdcoefs,freesurf,ybndflag,elas_only,&ninfo);

         wroutP3(&outp,pvptr[1],medptr[1],iy-1,it,&rpars_p3);
         if(iy == 2 && ninfo.minusId_y == -1)
            wroutP3(&outp,pvptr[0],medptr[0],iy-2,it,&rpars_p3);
         if(iy == ny-2 && ninfo.plusId_y == -1)
            {
            wroutP3(&outp,pvptr[2],medptr[2],iy,it,&rpars_p3);
            wroutP3(&outp,pvptr[3],medptr[3],iy+1,it,&rpars_p3);
            }
         }

/*
for(iy=0;iy<ny;iy++)
   {
   for(iv=0;iv<9;iv++)
      {
      for(iz=0;iz<nz;iz++)
         {
	 for(ix=0;ix<nx;ix++)
	    {
	    ip = ix + iz*nx + iv*nx*nz + iy*9*nx*nz;
            if(isnan(pvfield[ip]) != 0)
	       {
	       fprintf(stderr,"AFTER_T it= %d ix= %d iy= %d iz= %d iv= %d pv= %e\n",it,ix,iy,iz,iv,pvfield[ip]);
	       fflush(stderr);
	       }
	    }
         }
      }
   }
*/

   /* exchange stresses */

      rt1 = times(&tu1);
      blen = exchange_pvP3(pvfield,pvbuf,&ninfo,it,TAU);
      dtimes(tu1,rt1,&dtu,&drt,&tsize2,blen);

/*
for(iy=0;iy<ny;iy++)
   {
   for(iv=0;iv<9;iv++)
      {
      for(iz=0;iz<nz;iz++)
         {
	 for(ix=0;ix<nx;ix++)
	    {
	    ip = ix + iz*nx + iv*nx*nz + iy*9*nx*nz;
            if(isnan(pvfield[ip]) != 0)
	       {
	       fprintf(stderr,"AFTER_ALL it= %d ix= %d iy= %d iz= %d iv= %d pv= %e\n",it,ix,iy,iz,iv,pvfield[ip]);
	       fflush(stderr);
	       }
	    }
         }
      }
   }

   if(it==1)
   {
      MPI_Barrier(MPI_COMM_WORLD);
   mpi_exit(-1);
   }
*/

   /* if seismogram output is used, write temporary buffer to file */

      soutptr = &(outp.spnts);
      if(outp.nseis && soutptr->iflag == 1)
         {
         soutptr->flushbuf = 1;
         wroutP3(&outp,pvptr[3],medptr[3],-1,-1,&rpars_p3);
         soutptr->flushbuf = 0;
         }

   /* if xy time slice output is used, write temporary buffer to file */

      tsoutptr = &(outp.xyslice);
      if(outp.ts_xy && tsoutptr->iflag == 1)
         {
         tsoutptr->flushbuf = 1;
         wroutP3(&outp,pvptr[3],medptr[3],-1,-1,&rpars_p3);
         tsoutptr->flushbuf = 0;
         }

   /* if xz time slice output is used, write temporary buffer to file */

      tsoutptr = &(outp.xzslice);
      if(outp.ts_xz && tsoutptr->iflag == 1)
         {
         tsoutptr->flushbuf = 1;
         wroutP3(&outp,pvptr[3],medptr[3],-1,-1,&rpars_p3);
         tsoutptr->flushbuf = 0;
         }

   /* if yz time slice output is used, write temporary buffer to file */

      tsoutptr = &(outp.yzslice);
      if(outp.ts_yz && tsoutptr->iflag == 1)
         {
         tsoutptr->flushbuf = 1;
         wroutP3(&outp,pvptr[3],medptr[3],-1,-1,&rpars_p3);
         tsoutptr->flushbuf = 0;
         }

   /* if moment tensor output is used, write temporary buffer to output file */

      sgtoutptr = &(outp.sgtpnts);
      if(outp.sgtout && sgtoutptr->iflag == 1)
         {
         sgtoutptr->flushbuf = 1;
         wroutP3(&outp,pvptr[3],medptr[3],-1,-1,&rpars_p3);
         sgtoutptr->flushbuf = 0;
         }

   /* if damping boundaries are used, add random field to prevent underflow */

      if((elas_only == 0 || medinp.dampwidth > 0) && (it+1)%(rfinc) == 0)
         {
         for(iy=0;iy<ny;iy++)
            {
            rfptr = randf + rfcnt*nx*nz;
            add_random_field(pvfield+iy*N_WAVE_VARS*nx*nz,rfptr,nx,nz);

            rfcnt++;
            if(rfcnt == rfspan)
               rfcnt = 0;
            }
         }

      if(intfac.go)     /* check to continue interfacing */
         check_infc_goflag(&intfac);

      if(rstart_p3.enable_flag == 1 && ((it+1)%(rstart_p3.itinc)==0))
         dump_restartP3(&rstart_p3,name,ninfo.segmentId,&rpars_p3,&modstor,pvfield,it);

      if(dumpinfo.enable_flag == 1 && (((it+1)%(dumpinfo.itinc)==0) || (it+1) >= nt))
         dump_files(&outp,ninfo.segmentId,&dumpinfo);

      }  /* end of restart check */

   if(report && ((it+1)%(report)==0))
      print_report(rtime0,&rtime1,&use2,&use1,&dtu,&drt,&tsize2,&tsize1,it+1);

   } /* end of time step */

if(intface)     /* remove temporary interfacing files */
   rm_infc_files(&intfac);

if(dumpinfo.enable_flag == 1)   /* remove temporary output files */
   rm_dump_files(&outp,ninfo.segmentId,&dumpinfo);

total_report(rtime0,&use2,tsize2);

mpi_final("PROGRAM emod3d-mpi IS FINISHED");
exit(0);
}
