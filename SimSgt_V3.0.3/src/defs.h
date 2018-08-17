#define MEM_LIMIT     2048 /* internal memory limit in Mega-bytes */
#define ORDER            4
#define ORDERX           2 /* spatial differencing order below iz=izord2 */
#define N_WAVE_VARS     15 /* number of wave field variables (9 + 6 mem vars)*/
#define N_MED_VARS      13 /* number of media field variables stored in memory */
#define N_MEDPROF_VARS   6 /* number of media variables in medprof structure */
#define NPROFS          20 /* max. number of media profiles */
#define BASINPROF        1 /* != 0 -> scale vel. profiles with basin depth */
#define MAXSTAT        512 /* max. number of individual output stations */
#define WAVE_FIELD       0 /* flag for wave field storage (case pv) */
#define MEDIA_FIELD      1 /* flag for media field storage (case med) */

/* number of bytes used in station name for single file seismo output */
#define NBYTE_STATNAME   8

/* number of characters for file and directory pathnames */
#define DIR_STR_LEN   1024
#define FILE_STR_LEN   2048
 
/* dummy flag that is not used as a flag for anything else */
#define DUMY            -999
 
/* flags for P or S wave absorbing boundaries */
#define PABS            1 /* do P waves */
#define SABS            2 /* do S waves */
#define VAVG            3 /* average P and S velocities (near edges and
                             corner points) */
 
/* flags for doing corner points in abs() */
#define NOCRN           0 /* don't do corner point (interior of planes
                             at iy=0,ny-1) */
#define CRN_2D          1 /* xz plane corner point */
#define CRN_3D          2 /* iy=0,ny-1 corner edge, corners of 3D model space */ 
/* flags for free-surface */
#define VEL             1 /* solve for vz above free-surface */
#define TAU             2 /* set tzz, txz, tyz at and above free-surface */
 
/* flags for model boundaries */
#define XZERO           1 /* near ix=0 */
#define XN              2 /* near ix=nx-1 */
#define YZERO           3 /* near iy=0 */
#define YN              4 /* near iy=ny-1 */
#define ZZERO           5 /* near iz=0 */
#define ZN              6 /* near iz=nz-1 */
#define NBND_PAD        2 /* replicate outer 3 grid points of velocity model
			     for absorbing boundary stability (2 means 3 
			     because of C-type indexing) */
#define NBND_ZERO_PLANE 5 /* taper initial plane wave within 5 grid
		       	     points of x, y, z boundaries */
#define SKIP_PLANE_Y1	-1 /* skip time update on pvptr[1] plane (in overlap region) */
#define SKIP_PLANE_Y2	-2 /* skip time update on pvptr[2] plane (in overlap region) */
 
/* flags for x, y, or z derivative, functions diff() and diffx() */
#define XDERIV          1
#define YDERIV          2
#define ZDERIV          3
 
/* flags for different sources */
#define FSRC            1 /* body force -> add to velocities */
#define PSRC            2 /* pressure source -> add to normal stresses */
#define NSRCMAX         500 /* maximum number of double couple sources */

#define PI 3.141592653589740
#define AMP 3.0
#define RPERD 0.01745329
#define TSHIFT 1.0

#define         FLAT_CONST      298.256
#define         ERAD            6378.139

#if _FILE_OFFSET_BITS == 64

#define RDONLY_FLAGS    O_RDONLY | O_LARGEFILE
#define RDWR_FLAGS      O_RDWR | O_LARGEFILE
#define CROP_FLAGS      O_CREAT | O_RDWR | O_LARGEFILE
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR | O_LARGEFILE

#define print_lfmode(A,B)  fprintf(stderr,"***** Large file IO enabled, O_LARGEFILE=%d\n      File size limits:  soft=%.3g Gb  hard=%.3g Gb\n\n",O_LARGEFILE,(float)(A)/1.e+09,(float)(B)/1.e+09)

#else

#define RDONLY_FLAGS    O_RDONLY
#define RDWR_FLAGS      O_RDWR
#define CROP_FLAGS      O_CREAT | O_RDWR
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR

#define O_LARGEFILE    0

#define print_lfmode(A,B)  fprintf(stderr,"***** Large file IO disabled, O_LARGEFILE=%d\n      File size limits:  soft=%.3g Gb  hard=%.3g Gb\n\n",O_LARGEFILE,(float)(A)/1.e+09,(float)(B)/1.e+09)

#endif

#define NORMAL_STRESS_CENTERED    0
#define NULL_NODE_CENTERED    1
#define KBO_AT_VXNODE    2
#define KBO_AT_VZNODE    3
#define HARMONIC_AT_VZNODE    4

#define ENDIAN_BUF_LEN 32
