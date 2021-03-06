#ifndef RUPGEN_DEFS_H
#define RUPGEN_DEFS_H

/* Work stages */
typedef enum { RUPGEN_STAGE_CALC=0,
	       RUPGEN_STAGE_WRITE,
	       RUPGEN_STAGE_END } rg_stage_t;

/* Max filename */
#define MAX_FILENAME 128

/* Max id string length */
#define MAX_ID_LEN 128


#define DEFAULT_DHYPO_FRAC       0.75  /* hypo at 0.75 down-dip width */
#define DEFAULT_SHYPO_STEP       20.0  /* hypo spacing at 20 km along strike */
#define DEFAULT_SHYPO_MIN_OFF    1.0   /* hypos start at 1.0 km along strike */
#define DEFAULT_SLIPS_TO_HYPOS   2     /* #slip models = 2 times #hypos */

#define         DEFAULT_DT         0.1
#define         NTMAX              10000
#define         SOMERVILLE_FLAG    1
#define         MAI_FLAG           2
#define         INPUT_CORNERS_FLAG           -1
#define         MINSLIP            1.0e-02

#define DEFAULT_KMODEL    2
#define DEFAULT_FLIP_AT_SURFACE    1
#define DEFAULT_STRETCH_KCORNER    0
#define DEFAULT_CIRCULAR_AVERAGE    0
#define DEFAULT_MODIFIED_CORNERS    0
#define DEFAULT_TRUNCATE_ZERO_SLIP    1
#define DEFAULT_RAND_RAKE_DEGS    60.0
#define DEFAULT_SLIP_SIGMA    0.85

#define DEFAULT_TSFAC  -0.5
#define DEFAULT_TSFAC_COEF  1.8
#define DEFAULT_TSFAC_FACTOR    1

#define DEFAULT_RT_SCALEFAC    1

#define DEFAULT_VR_TO_VS_FRAC    0.8   /* vrup = 0.8 times local Vs */
#define DEFAULT_SHAL_VRUP_FRAC   0.7   /* shallow_vrup = 0.8 times back_vrup */

#define DEFAULT_DEPTH_SCALING_LEVEL  6.5
#define DEFAULT_DEPTH_SCALING_RANGE  1.5

#define DEFAULT_TOP_TAP         0.0
#define DEFAULT_SIDE_TAP         0.05   /* taper slip at 0.1*flen on sides */
#define DEFAULT_BOT_TAP          0.1   /* taper slip at 0.1*fwid on top/bot */

#define RDONLY_FLAGS    O_RDONLY
#define RDWR_FLAGS      O_RDWR
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR

#define NUM_GENSLIP_ARGS 7
#define MAX_STR_ARGV 512

#define GET_STATS 1
#define RUN_GENSLIP 2

#if _FILE_OFFSET_BITS == 64

#undef RDONLY_FLAGS
#undef RDWR_FLAGS
#undef CROPTR_FLAGS

#define RDONLY_FLAGS    O_RDONLY | O_LARGEFILE
#define RDWR_FLAGS      O_RDWR | O_LARGEFILE
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR | O_LARGEFILE

#endif

#endif

