#define MAX_PARAMS 100 /* maximum number of parameters to pass */
#define MAX_PARAM_LENGTH 64 /* maximum length of an individual parameter */

#define         NTMAX              10000
#define         SOMERVILLE_FLAG    1
#define         MAI_FLAG           2
#define         MINSLIP            1.0e-02

#define   DHYPO_FRAC       0.75     /* hypo at 0.75 down-dip width */
#define   SHYPO_STEP       20.0     /* hypo spacing at 20 km along strike */
#define   SHYPO_MIN_OFF    1.0      /* hypos start at 1.0 km along strike */
#define   SLIPS_TO_HYPOS   2    /* no. slip models = 2 times no. of hypos */

#define DEFAULT_VR_TO_VS_FRAC    0.8   /* vrup = 0.8 times local Vs */
#define DEFAULT_SHAL_VRUP_FRAC   0.8   /* shallow_vrup = 0.8 times back_vrup */
#define DEFAULT_TSFAC            -0.5  /* tinit is shifted -0.5 at max_slip */

#define RDONLY_FLAGS    O_RDONLY
#define RDWR_FLAGS      O_RDWR
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR

#define PSA_WRITE 0
#define PSA_FILE_FORWARD 1
#define PSA_NO_WRITE 2

#define MAXPERIODS 305
#define NUM_SCEC_PERIODS 44

#if _FILE_OFFSET_BITS == 64

#undef RDONLY_FLAGS
#undef RDWR_FLAGS
#undef CROPTR_FLAGS

#define RDONLY_FLAGS    O_RDONLY | O_LARGEFILE
#define RDWR_FLAGS      O_RDWR | O_LARGEFILE
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR | O_LARGEFILE

#endif
