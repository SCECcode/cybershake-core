#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "genslip.h"
#include "misc.h"
#include "getpar.h"

int compareFloat(const void* f1, const void* f2) {
        return (*(float*)f1 - *(float*)f2);
}

/*
************************************************************************************

   03/06/2000: version 2.1

   * Added option to read input fault description from Generic Slip Format (GSF)
     file.

   * Added option to specify "Brune" type STF.

************************************************************************************

   07/08/2009: version 3.0

   * Added flag to make SRF output optional.

   * Added tsfac depth and slip scaling for initiation time.

************************************************************************************

   2010/04/09: version 3.1

   * Added option to adjust rupture speed at segment boundaries, see seg_delay

************************************************************************************

   2010/07/21: version 3.2

   * Added option to adjust rise time with sqrt(slip) [previously done using generic_slip2srf]
     This is done using the parameter stfparams.rt_scalefac=1 (=0 gives no slip scaling)

   * Set the following as defaults:

	 kmodel=2
         flip_at_surface=1
	 stretch_kcorner=0
	 circular_average=0
	 modified_corners=0
         truncate_zero_slip=1
         rand_rake_degs=60.0
	 slip_sigma=0.85

	 tsfac_coef=1.8
	 tsfac_factor=1
	 tsfac=-ts_coef*1.0e-09*exp(log(moment)/3.0)

	 rvfrac=0.8
	 shal_vrup=0.7

         shal_vrup_dep=6.5
         shal_vrup_deprange=1.5

         tsfac_dep=6.5
         tsfac_deprange=1.5

         stfparams.trise=1.6*1.0e-09*exp(log(moment)/3.0)
         stfparams.rt_scalefac=1
         stfparams.risetimedep=6.5
         stfparams.risetimedep_range=1.5
         stfparams.risetimefac=2.0

         side_taper=0.05
         bot_taper=0.1
         top_taper=0.0

************************************************************************************

   2013/10/01: version 3.2.1

   * Changed default/recommended usage to directly produce SRF file (instead of GSF
     and then using genericslip2srf to get SRF).  All options existed in previous
     versions. To directly produce SRF, use these steps (see README also):

          1) Modify output format by setting:

             write_srf=1
             write_gsf=0

          so that the code will directly produce the SRF file without having to
          make subsequent call to 'generic_slip2srf'

          2) Additional parameters needed for SRF output are:

             stype=$STYPE
             dt=$DT
             risetime=$RISETIME
             plane_header=1
             risetimefac=$RTFAC
             risetimedep=$RTDEP
             risetimedep_range=$RTDEP_RANGE
             rt_scalefac=$RT_SCALEFAC

          These should already be defined in the script.

          3) The default setting with "write_srf=1" will write to stdout.  So the
          call to the code should be something like:

          genslip-v3.2.1 read_erf=0 write_srf=1 ... rt_scalefac=$RT_SCALEFAC > $SRFDIR/$SRFFILE

          And then remove the subsequent call to 'generic_slip2srf'.

   * Additional modifications based on results of SCEC BBP testing to reduce
     the strength/coherency of longer period radiation.

   * Changed the default tapering parameters to:

         side_taper=0.02
         bot_taper=0.0
         top_taper=0.0

   * Added following parameters to adjust rise time along bottom portion of fault
     (similar to shallow adjustments).  Defaults are:

         stfparams.deep_risetimedep=16.5      # or hypo_depth if greater
         stfparams.deep_risetimedep_range=1.5
         stfparams.deep_risetimefac=1.5

   * Also added parameter to set a minimum background slip level given as a percentage
     of the average slip amount (basically fills-in very low/zero slip patches with
     long rise time low slip. Testing indicates that this has little impact on overall
     results, so for now the default setting is to not use this option:

         slip_water_level=-1

************************************************************************************

   2013/11/19: version 3.3

   * Changed default/recommended slip rate function to 'Mliu' which is modified version
     of UCSB function.  See gen_Mliu_stf() in srf_subs.c for details.

   * Added parameter 'risetime_coef' which is used to calculate average risetime
     from moment using:

         stfparams.trise = risetime_coef*1.0e-09*exp(log(mom)/3.0);

     This option only works if risetime=-1 which is now the default.

   * Changed default/recommended coefficient for average rise time relation from '1.6'
     to '1.45' =>

         risetime_coef = 1.45

     This can be overridden by specifying as a getper parameter, e.g., for use in
     CEUS, use risetime_coef=3.75

   * Added 'alphaT' parameter for scaling of rise time for dip & rake.  Update of alphaT
     parameter defined by Graves and Pitarka (2010) =>

          alphaT = 1.0/(1.0 + fD*fR*Calpha)

          where Calpha = 0.1   (max value when fD=fR=1.0)
         
	        fD = 1.0 - (dip - 45.0)/45.0      for 45.0 < dip <= 90.0
	        fD = 1.0                          for dip <= 45.0
	        fD = 0.0                          otherwise
         
	        fR = 1.0 - (|rake - 90.0|)/90.0   for 0.0 <= rake <= 180.0
	        fR = 0.0                          otherwise

          Note: should have 0 <= dip <= 90 and -180 <= rake <= 180

     Default for Calpha is 0.1, but can be passed to code as getpar() variable.

   * Added random perturbations to tsfac so that it is not 1:1 correlated with slip.
     Perturbations are log normal with ln(sigma)=tsfac_rand*tsfac.  Default for
     tsfac_rand=0.2.

   * Added random perturbations to risetime so that it is not 1:1 correlated with
     sqrt(slip).  Perturbations are log normal with ln(sigma)=rt_rand*trise.  Default
     for rt_rand=0.5.

************************************************************************************
*/
#ifdef _USE_MEMCACHED
int mc_genslip(int ac, char** av, rg_stats_t *stats, struct standrupformat* srf, int state, char* memcached_server)
#else
int genslip(int ac,char **av, rg_stats_t *stats, struct standrupformat* srf, int state)
#endif
{
FILE *fpr, *fpw;
struct complex *cslip, *crake;

float flen, fwid, dx, dy, sval;
float dkx, dky, kx, ky, xx, yy, zz;
float cosA, sinA, cosD, sinD, dd, sn, se;
float sum, sigma, neg_sum;
int nx, ny;
int i, j, k, ip, nt6, it;
char infile[1024], str[1024];
char init_slip_file[1024], outfile[1024];

float bigM;
float xl, yl;

int flip_at_surface = DEFAULT_FLIP_AT_SURFACE;
int stretch_kcorner = DEFAULT_STRETCH_KCORNER;
int ny_in, ny_str, ny_end;

float pi = 3.14159265;
double rperd = 0.017453293;

float mag;
float side_taper = DEFAULT_SIDE_TAP;
float bot_taper = DEFAULT_BOT_TAP;
float top_taper = DEFAULT_TOP_TAP;

int truncate_zero_slip = DEFAULT_TRUNCATE_ZERO_SLIP;

int generate_seed = 1;
long seed = 0;
long starting_seed;
int kmodel = MAI_FLAG;   /* default is mai */
int circular_average = DEFAULT_CIRCULAR_AVERAGE;
int modified_corners = DEFAULT_MODIFIED_CORNERS;

float kx_corner, ky_corner, xmag_exponent, ymag_exponent;

float mag_area_Acoef = -1.0;
float mag_area_Bcoef = -1.0;
int use_median_mag = 0;

int uniformgrid_hypo = 1;
int random_hypo = 1;
int uniform_prob4hypo = 0;
struct hypo_distr_params hpar_as, hpar_dd;
int calc_shypo = 1;

float shypo = -1.0e+15;
float dhypo = -1.0e+15;

int nrup_min = 10;
float target_hypo_spacing = 4.5;
float hypo_s0, hypo_d0, hypo_ds, hypo_dd;
int ih_scnt, ih_dcnt, nhypo_s, nhypo_d;

float avgstk, rake, rt, tsmin;
float xhypo;
struct velmodel rvmod;
double rayp, rupt_rad;

int seg_delay = 0;
int nseg_bounds, ig;
float delh, hx, hg, gwid2, *gwid, shypo_mseg, *xseg, *rvfac_seg;

float rup_delay = 0.0;

float rmin, rmax, ravg, rfac, rmed;
float dkx_rk, dky_rk;
float *sort_rake, *psrc_rake;

float set_rake = -999.0;
float rand_rake_degs = DEFAULT_RAND_RAKE_DEGS;
int test;

float shypo_step = DEFAULT_SHYPO_STEP;
float shypo_min_off = DEFAULT_SHYPO_MIN_OFF;
float dhypo_frac = DEFAULT_DHYPO_FRAC;
int slips_to_hypos = DEFAULT_SLIPS_TO_HYPOS;

float sh0;
int nh = -1;
int ns = -1;
int ih, js;

float rvfrac = DEFAULT_VR_TO_VS_FRAC;
float shal_vrup = DEFAULT_SHAL_VRUP_FRAC;
float htol = 0.1;

float smax, sf, sfac0, sabs, sden, snum;
float savg;
float target_savg = -1.0;

float depmin1, depmax1, depmin2, depmax2;
float deep_vrup_dep = 99.5;
float shal_vrup_dep = DEFAULT_DEPTH_SCALING_LEVEL;
float shal_vrup_deprange = DEFAULT_DEPTH_SCALING_RANGE;

float tsfac = -1.0e+15;;
float tsfac_coef = DEFAULT_TSFAC_COEF;
float tsfac_factor = DEFAULT_TSFAC_FACTOR;
float tsfac_dep = DEFAULT_DEPTH_SCALING_LEVEL;
float tsfac_deprange = DEFAULT_DEPTH_SCALING_RANGE;
float tsfac_rand = 0.2;

/* 20131119 RWG:

   If risetime=-1 as input (default) code will calculate average risetime
   from moment using:

      risetime = risetime_coef*1.0e-09*exp(log(mom)/3.0);

   Default rise time relation coefficient modified from somerville, 1.6 -> 1.45
   risetime_coef = 1.45
*/

float risetime_coef = 1.45;
float deep_risetimedep_saved;

float sd_rand = -1.0;
float fzero = 0.0;
float fone = 1.0;

float sigfac;
float slip_sigma = DEFAULT_SLIP_SIGMA;

float dtop, avgdip, mom, mag_med;
struct velmodel vmod;
char velfile[1024];

int read_erf = 1;

struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

int read_gsf = 0;
int write_gsf = 0;

int write_srf = 1;

struct generic_slip gslip;

float *stf, elon, elat;
int ntot;

struct pointsource *psrc;
struct stfpar2 stfparams;

int outbin = 0;
int writeout = 1;
int doslip = -1;
int dohypo = -1;

float slip_water_level = -1;
float slipmin;

float Calpha = 0.1;
float alphaT, fD, fR, avgrak;

velfile[0] = '\0';

init_slip_file[0] = '\0';

sprintf(infile,"stdin");
sprintf(outfile,"stdout");


sprintf(srf->version,"1.0");

rvfac_seg = NULL;
xseg = NULL;
gwid = NULL;

stfparams.dt = DEFAULT_DT;
stfparams.nt = NTMAX;
stfparams.trise = -1.0;
stfparams.rt_scalefac = DEFAULT_RT_SCALEFAC;
stfparams.risetimedep = DEFAULT_DEPTH_SCALING_LEVEL;
stfparams.risetimedep_range = DEFAULT_DEPTH_SCALING_RANGE;
stfparams.risetimefac = 2.0;
stfparams.deep_risetimedep = 16.5;
stfparams.deep_risetimedep_range = 1.5;
stfparams.deep_risetimefac = 1.5;
stfparams.rt_rand = 0.5;          /* ln(sigma) */
sprintf(stfparams.stype,"Mliu");  /* default is modified UCSB sincos */

/* RWG 2014-02-20 randomized hypocenter
   default probability tapering for randomized hypocenter */

/* along strike -> */
hpar_as.x0 = 0.2;	/* default tapering starts at 20% of fault length at end end */
hpar_as.x1 = 0.8;
hpar_as.f0 = 0.1;	/* default probability at edge is 10% of probability in middle of fault */
hpar_as.f1 = 0.1;

/* down dip -> */
hpar_dd.x0 = 0.4;	/* default tapering starts at 40% along top edge of fault */
hpar_dd.x1 = 0.8;	/* default tapering starts at 20% from fault bottom */
hpar_dd.f0 = 0.01;	/* default probability at top edge is 1% of probability in middle of fault */
hpar_dd.f1 = 0.1;	/* default probability at bottom edge is 10% of probability in middle of fault */

setpar(ac,av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);
getpar("doslip","d",&doslip);
getpar("dohypo","d",&dohypo);

getpar("velfile","s",velfile);

getpar("read_erf","d",&read_erf);
getpar("read_gsf","d",&read_gsf);

getpar("write_srf","d",&write_srf);

if(read_gsf == 1)
   {
   read_erf = 0;
   mstpar("mag","f",&mag);

   mstpar("nx","d",&nx);
   mstpar("ny","d",&ny_in);
   getpar("write_gsf","d",&write_gsf);
   }

if(read_erf == 0 && read_gsf == 0)
   {
   mstpar("mag","f",&mag);

   mstpar("nx","d",&nx);
   mstpar("ny","d",&ny_in);
   mstpar("dx","f",&dx);
   mstpar("dy","f",&dy);

   mstpar("dtop","f",&dtop);
   mstpar("strike","f",&avgstk);
   mstpar("dip","f",&avgdip);
   mstpar("rake","f",&rake);
   mstpar("elon","f",&elon);
   mstpar("elat","f",&elat);
   }

getpar("dt","f",&stfparams.dt);
getpar("nt","d",&stfparams.nt);
getpar("risetime","f",&stfparams.trise);
getpar("risetime_coef","f",&risetime_coef);
getpar("risetimefac","f",&stfparams.risetimefac);
getpar("risetimedep","f",&stfparams.risetimedep);
getpar("risetimedep_range","f",&stfparams.risetimedep_range);
getpar("rt_scalefac","f",&stfparams.rt_scalefac);
getpar("rt_rand","f",&stfparams.rt_rand);
getpar("stype","s",stfparams.stype);

getpar("deep_risetimefac","f",&stfparams.deep_risetimefac);
getpar("deep_risetimedep","f",&stfparams.deep_risetimedep);
getpar("deep_risetimedep_range","f",&stfparams.deep_risetimedep_range);

deep_risetimedep_saved = stfparams.deep_risetimedep;

getpar("shypo_step","f",&shypo_step);
getpar("shypo_min_off","f",&shypo_min_off);
getpar("dhypo_frac","f",&dhypo_frac);
getpar("slips_to_hypos","d",&slips_to_hypos);
getpar("ns","d",&ns);
getpar("nh","d",&nh);

getpar("shypo","f",&shypo);
getpar("dhypo","f",&dhypo);

/*XXXX*/
getpar("generate_seed","d",&generate_seed);

/* RWG 2014-04-24 uniform grid of hypocenters */
getpar("uniformgrid_hypo","d",&uniformgrid_hypo);

/* RWG 2014-02-20 randomized hypocenter */
getpar("random_hypo","d",&random_hypo);
getpar("uniform_prob4hypo","d",&uniform_prob4hypo);

getpar("nrup_min","d",&nrup_min);
/* RWG 2014-04-24 target spacing of hypocenters, replaces nrup_scale_fac */
getpar("target_hypo_spacing","f",&target_hypo_spacing);

getpar("hypo_taperperc_left","f",&hpar_as.x0);
getpar("hypo_taperval_left","f",&hpar_as.f0);
getpar("hypo_taperperc_right","f",&hpar_as.x1);
getpar("hypo_taperval_right","f",&hpar_as.f1);

getpar("hypo_taperperc_top","f",&hpar_dd.x0);
getpar("hypo_taperval_top","f",&hpar_dd.f0);
getpar("hypo_taperperc_bot","f",&hpar_dd.x1);
getpar("hypo_taperval_bot","f",&hpar_dd.f1);
/*XXXX*/

getpar("rvfrac","f",&rvfrac);
getpar("shal_vrup","f",&shal_vrup);
getpar("shal_vrup_dep","f",&shal_vrup_dep);
getpar("shal_vrup_deprange","f",&shal_vrup_deprange);

getpar("rup_delay","f",&rup_delay);

getpar("tsfac","f",&tsfac);
getpar("tsfac_factor","f",&tsfac_factor);
getpar("tsfac_dep","f",&tsfac_dep);
getpar("tsfac_deprange","f",&tsfac_deprange);
getpar("tsfac_rand","f",&tsfac_rand);
getpar("sd_rand","f",&sd_rand);

getpar("rand_rake_degs","f",&rand_rake_degs);
getpar("set_rake","f",&set_rake);

getpar("target_savg","f",&target_savg);

getpar("kmodel","d",&kmodel);

if(kmodel < 0)
   kmodel = INPUT_CORNERS_FLAG;
if(kmodel != MAI_FLAG && kmodel != INPUT_CORNERS_FLAG && kmodel < 100)
   kmodel = SOMERVILLE_FLAG;

if(kmodel == INPUT_CORNERS_FLAG)
   {
   mstpar("kx_corner","f",&kx_corner);
   mstpar("ky_corner","f",&ky_corner);

   xmag_exponent = 0.5;
   ymag_exponent = 0.5;
   getpar("xmag_exponent","f",&xmag_exponent);
   getpar("ymag_exponent","f",&ymag_exponent);
   }

getpar("mag_area_Acoef","f",&mag_area_Acoef);
getpar("mag_area_Bcoef","f",&mag_area_Bcoef);
getpar("use_median_mag","d",&use_median_mag);

getpar("modified_corners","d",&modified_corners);
getpar("circular_average","d",&circular_average);
getpar("flip_at_surface","d",&flip_at_surface);
getpar("stretch_kcorner","d",&stretch_kcorner);
getpar("seed","d",&seed);
getpar("side_taper","f",&side_taper);
getpar("bot_taper","f",&bot_taper);
getpar("top_taper","f",&top_taper);

getpar("truncate_zero_slip","d",&truncate_zero_slip);
getpar("slip_water_level","f",&slip_water_level);

getpar("slip_sigma","f",&slip_sigma);

getpar("init_slip_file","s",init_slip_file);

getpar("seg_delay","d",&seg_delay);
if(seg_delay == 1)
   {
   mstpar("nseg_bounds","d",&nseg_bounds);

   xseg = (float *)_check_malloc(nseg_bounds*sizeof(float));
   rvfac_seg = (float *)_check_malloc(nseg_bounds*sizeof(float));
   gwid = (float *)_check_malloc(nseg_bounds*sizeof(float));

   for(ig=0;ig<nseg_bounds;ig++)
      {
      rvfac_seg[ig] = 0.5;    /* default is 50% reduction of rupture speed at seg boundaries */
      gwid[ig] = 4.0;	/* default width of delay zone is 4 km */
      }

   mstpar("xseg","vf",xseg);
   getpar("rvfac_seg","vf",rvfac_seg);
   getpar("gwid","vf",gwid);

   /* change rvfac_seg to equivalent distance adjustment */
   for(ig=0;ig<nseg_bounds;ig++)
      rvfac_seg[ig] = 1.0/rvfac_seg[ig] - 1.0;
   }

endpar();

psrc = (struct pointsource *)NULL;
gslip.np = -1;
gslip.spar = (struct slippars *)NULL;

if(read_erf == 1)
#ifdef _USE_MEMCACHED
   psrc = _mc_read_ruppars(infile,psrc,&mag,&nx,&ny_in,&dx,&dy,&dtop,&avgstk,&avgdip,&elon,&elat,memcached_server);
#else

   psrc = _read_ruppars(infile,psrc,&mag,&nx,&ny_in,&dx,&dy,&dtop,&avgstk,&avgdip,&elon,&elat);
#endif
else if(read_gsf == 1)
   psrc = _read_gsfpars(infile,psrc,&gslip,&dx,&dy,&dtop,&avgdip);
else
   psrc = _set_ruppars(psrc,&mag,&nx,&ny_in,&dx,&dy,&dtop,&avgstk,&avgdip,&rake,&elon,&elat);

flen = nx*dx;
fwid = ny_in*dy;

bigM = log(10.0);
mom = exp(bigM*1.5*(mag + 10.7));
mag_med = 3.98 + log(flen*fwid)/bigM;
/* update 12/2005 */
mag_med = 3.87 + 1.05*log(flen*fwid)/bigM;

if(tsfac < -1.0e+10)
   tsfac = -tsfac_coef*1.0e-09*exp(log(mom)/3.0);

/* update 12/2007 */
if(mag_area_Acoef < 0.0)
   {
   mag_area_Acoef = 4.0;
   mag_area_Bcoef = 1.0;
   }
mag_med = mag_area_Acoef + mag_area_Bcoef*log(flen*fwid)/bigM;

if(mag < 0.0 || use_median_mag)
   mag = mag_med;

_init_plane_srf(srf,&gslip,&elon,&elat,nx,ny_in,&flen,&fwid,&dx,&dy,&avgstk,&avgdip,&dtop,&shypo,&dhypo);

if(velfile[0] != '\0')
   read_Fvelmodel(velfile,&vmod);
   /*
   read_velmodel(velfile,&vmod);
   */
else
   default_velmodel(&vmod);

/*
conv2vrup(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup);
*/
depmin1 = shal_vrup_dep - shal_vrup_deprange;
depmax1 = shal_vrup_dep + shal_vrup_deprange;
depmin2 = deep_vrup_dep - shal_vrup_deprange;
depmax2 = deep_vrup_dep + shal_vrup_deprange;
/*
conv2vrup_dd2(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup,&depmin1,&depmax1,&depmin2,&depmax2);
*/
conv2vrup_dd(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup,&depmin1,&depmax1);

/* generic fault length/width scaling */
xl = flen;
yl = fwid;

if(xl > 2.0*yl)
   xl = 2.0*yl;
if(yl > 2.0*xl)
   yl = 2.0*xl;

if(kmodel == MAI_FLAG) /* mai scaling */
   {
   /* add factor of log10(2*pi) */
   xl = exp(bigM*(0.5*mag - 2.50 + 0.79818));
   yl = exp(bigM*(0.3333*mag - 1.50 + 0.79818));

   if(circular_average)
      yl = xl;
   }

if(kmodel == SOMERVILLE_FLAG) /* somerville scaling */
   {
   xl = exp(bigM*(0.5*mag - 1.72));
   yl = exp(bigM*(0.5*mag - 1.93));

   if(circular_average)
      {
      xl = exp(bigM*(0.5*mag - 1.825));
      yl = xl;
      }
   }

if(modified_corners)
   {
   xl = exp(bigM*(0.5*mag - 2.00));
   yl = exp(bigM*(0.5*mag - 2.00));
   }

if(kmodel == INPUT_CORNERS_FLAG)
   {
   xl = exp(bigM*(xmag_exponent*mag - kx_corner));
   yl = exp(bigM*(ymag_exponent*mag - ky_corner));
   }

if(flip_at_surface < 0)
   {
   if(dtop <= 0.100001)
      flip_at_surface = 1;
   else
      flip_at_surface = 0;
   }

ny = ny_in;
if(flip_at_surface == 1)
   {
   ny = 2*ny_in;

   if(stretch_kcorner)
      {
      xl = sqrt(2.0)*xl;
      yl = sqrt(2.0)*yl;
      }
   }

dkx = 1.0/(nx*dx);
dky = 1.0/(ny*dy);

if(stfparams.trise < 0.0)
   stfparams.trise = risetime_coef*1.0e-09*exp(log(mom)/3.0);

/* 20131118: rise time modification for rake and dip => update of GP2010 alphaT parameter */

Calpha = 0.1;
         
fD = 0.0;
if(avgdip <= 90.0 && avgdip > 45.0)
   fD = 1.0 - (avgdip - 45.0)/45.0;
else if(avgdip <= 45.0 && avgdip >= 0.0)
   fD = 1.0;

avgrak = 0.0;
for(j=0;j<ny_in*nx;j++)
   avgrak = avgrak + psrc[j].rak;

avgrak = avgrak/((float)(nx*ny_in));
while(avgrak < -180.0)
   avgrak = avgrak + 360.0;
while(avgrak > 180.0)
   avgrak = avgrak - 360.0;
         
fR = 0.0;
if(avgrak <= 180.0 && avgrak >= 0.0)
   fR = 1.0 - sqrt((avgrak - 90.0)*(avgrak - 90.0))/90.0;
         
alphaT = 1.0/(1.0 + fD*fR*Calpha);
stfparams.trise = alphaT*stfparams.trise;

/* alphaT done */

cslip = (struct complex *) _check_malloc (nx*ny*sizeof(struct complex));

dkx_rk = 1.0/(nx*dx);
dky_rk = 1.0/(ny_in*dy);
crake = (struct complex *) _check_malloc (nx*ny_in*sizeof(struct complex));
sort_rake = (float *) _check_malloc (nx*ny_in*sizeof(float));
psrc_rake = (float *) _check_malloc (nx*ny_in*sizeof(float));

for(j=0;j<ny_in*nx;j++)
   psrc_rake[j] = psrc[j].rak;

/* XXXX */
/* RWG 2014-02-20 randomized hypocenter */
if(generate_seed == 1) {
   gseed(&seed,psrc,&flen,&fwid,&dtop,&mag);
}

/* RWG 2014-03-21 set starting seed */
starting_seed = seed;

if(random_hypo == 1)
   {
   hpar_as.x0 = hpar_as.x0*flen;
   hpar_as.x1 = hpar_as.x1*flen;
   hpar_as.xlen = flen;
   hpar_as.xshift = -0.5*flen;

   hpar_dd.x0 = hpar_dd.x0*fwid;
   hpar_dd.x1 = hpar_dd.x1*fwid;
   hpar_dd.xlen = fwid;
   hpar_dd.xshift = 0.0;
   }

if(read_erf == 1)
   {
   if(uniformgrid_hypo == 1) /* RWG 20140424 added option for uniform grid of hypocenters */
      {
      calc_shypo = 0;

      nhypo_s = (int)((float)(flen/target_hypo_spacing));
      hypo_ds = target_hypo_spacing;
      if(nhypo_s < 3)
         {
         nhypo_s = 3;
         hypo_ds = flen/(nhypo_s+1);
         }
      hypo_s0 = 0.5*(flen - (nhypo_s-1)*hypo_ds) - 0.5*flen;

      nhypo_d = (int)((float)(fwid/target_hypo_spacing));
      hypo_dd = target_hypo_spacing;
      if(nhypo_d < 3)
         {
         nhypo_d = 3;
         hypo_dd = fwid/(nhypo_d+1);
         }
      hypo_d0 = 0.5*(fwid - (nhypo_d-1)*hypo_dd);

      ns = nhypo_s*nhypo_d;
      ih_scnt = 0;
      ih_dcnt = 0;

      nh = 1;
      }
   else if(random_hypo == 1)
      {
      calc_shypo = 0;

      /* RWG 20140424 OLD way
        ns = (int)(0.1*nrup_scale_fac*flen*fwid + 0.5);
               if target_hypo_spacing ~= 4.5, then nrup_scale_fac ~= 0.5
     */
      ns = (int)(flen*fwid/(target_hypo_spacing*target_hypo_spacing) + 0.5);
      if(ns < nrup_min)
         ns = nrup_min;

      nh = 1;
      }
   else
      {
      calc_shypo = 1;
      shypo = -1.0e+15;
      dhypo = -1.0e+15;

      if(nh < 0)
         nh = (int)((flen-2.0*shypo_min_off)/shypo_step) + 1;
      if(ns < 0)
         {
         ns = slips_to_hypos*nh;
         if(nx == 1 && ny_in == 1)
            ns = 1;
         }

      sh0 = 0.5*(flen - (nh-1)*shypo_step);
      }
   }

if(shypo > -1.0e+14)
   calc_shypo = 0;

if(dhypo < -1.0e+14)
   dhypo = dhypo_frac*fwid;
/* XXXX */

//Added for RupGen API
stats->numslip = ns;
stats->numhypo = nh;

if (state==GET_STATS) {
   free(psrc);
   return 0;
}

if (dohypo>=nh) {
        fprintf(stderr,"Error:  You requested hypocenter ID %d but from infile %s genslip calculates %d hypocenters.\n", dohypo, infile, nh);
        exit(1);
}

if (doslip>=ns) {
        fprintf(stderr,"Error:  You requested slip ID %d but from infile %s genslip calculates %d slips.\n", doslip, infile, ns);
        exit(2);
}

if (ns<-1 || nh<-1) {
        fprintf(stderr, "Error: hypocenter and slip ID need to be >= -1.");
        exit(3);
}

printf("Calculating for slip %d, hypo %d.\n", doslip, dohypo);
//end add

fprintf(stderr,"mag= %.2f median mag= %.2f nslip= %d nhypo= %d\n",mag,mag_med,ns,nh);
fprintf(stderr,"nx= %d ny= %d dx= %10.4f dy= %10.4f\n",nx,ny_in,dx,dy);

for(js=0;js<ns;js++)    /* loop over slip/rupture realizations */
   {
   if (doslip!=-1 && js>doslip) {
	//We've already done the slip we need, break out of the loop
	break;
   }

   if (doslip!=-1 && js<doslip) {
	//Haven't gotten to the right slip yet, continue
	continue;
   }

/*
 *    RWG 2014-03-21 set initial seed using increments of starting seed.
 *       Allows reproducability of ruptures without having to generate entire set.
 *       */

   seed = starting_seed;
   for(k=0;k<10*js;k++)
      sval = _sfrand(&seed);

   fprintf(stderr,"js= %d seed= %d ran= %10.6f\n",js,seed,sval);

   if(uniformgrid_hypo == 1) /* RWG 20140424 added option for uniform grid of hypocenters */
      {
      /*if(ih_scnt == nhypo_s)
          {
          ih_scnt = 0;
          ih_dcnt++;
          }
       shypo = hypo_s0 + ih_scnt*hypo_ds;
       dhypo = hypo_d0 + ih_dcnt*hypo_dd;

       ih_scnt++;*/
      shypo = hypo_s0 + (js%nhypo_s)*hypo_ds;
      dhypo = hypo_d0 + (js/nhypo_s)*hypo_dd;
      }


/* RWG 2014-02-20 randomized hypocenter */
   if(random_hypo == 1)
      {
      if(uniform_prob4hypo == 1)
         rhypo_uniform(&seed,&shypo,&dhypo,&flen,&fwid);
      else
         {
         shypo = rhypo1_lintaper(&seed,&hpar_as);
         dhypo = rhypo1_lintaper(&seed,&hpar_dd);
	 }
      }

   printf("shypo=%f, dhypo=%f\n", shypo, dhypo);
   init_slip_IO_fftw(cslip,nx,ny,&dx,&dy,flip_at_surface,init_slip_file);

   //fftwf_execute(plan_fwd);
   fft2d_fftw(cslip,nx,ny,-1,&dx,&dy);
   /*for(j=0;j<nx*ny;j++) {
          printf("cslip[%d].re=%f, cslip[j].im=%f\n", j, cslip[j].re, cslip[j].im);
   }*/

   if(kmodel < 100)
      kfilt_lw(cslip,nx,ny,&dkx,&dky,&xl,&yl,&seed,kmodel,&flen,&fwid);
   fft2d_fftw(cslip,nx,ny,1,&dkx,&dky);
   //fftwf_execute(plan_inv);
   /*float normalize = 1.0/((float)(nx*ny));
   for(i=0; i<nx; i++) {
	for(j=0; j<ny; j++) {
		cslip[i*ny+j].re = normalize*cslip_fftwf[i*ny+j][0];
		cslip[i*ny+j].im = normalize*cslip_fftwf[i*ny+j][1];
	}
   }*/

/*
   truncate any negative slip values => should double check that spectra is
   not changed too much
*/

   ny_str = 0;
   ny_end = ny_in;

   sum = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         sum = sum + cslip[ip].re;
         //sum = sum + cslip[ip][0];
         }
      }
   sum = sum/(float)(nx*(ny_end-ny_str));

   sigma = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         sigma = sigma + (cslip[ip].re-sum)*(cslip[ip].re-sum);
         //sigma = sigma + (cslip[ip][0]-sum)*(cslip[ip][0]-sum);
         }
      }
   sigma = sqrt(sigma/((ny_end-ny_str)*nx))/sum;

   if(slip_sigma > 0.0)
      {
      sigfac = slip_sigma/sigma;
      for(j=ny_str;j<ny_end;j++)
         {
         for(i=0;i<nx;i++)
            {
	    ip = i + j*nx;
            cslip[ip].re = sigfac*(cslip[ip].re - sum) + sum;
            //cslip[ip][0] = sigfac*(cslip[ip][0] - sum) + sum;
            }
         }
      }

   neg_sum = 0.0;
   sum = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         if(cslip[ip].re < 0.0 && truncate_zero_slip)
         //if(cslip[ip][0] < 0.0 && truncate_zero_slip)
	    {
	    neg_sum = neg_sum - cslip[ip].re;
	    //neg_sum = neg_sum - cslip[ip][0];
            cslip[ip].re = 0.0;
            //cslip[ip][0] = 0.0;
	    /*
	    fprintf(stderr,"slip= %13.5e\n",cslip[ip].re);
	    */
	    }
	 sum = sum + cslip[ip].re;
	 //sum = sum + cslip[ip][0];
	 }
      }

   sum = sum/(float)(nx*(ny_end-ny_str));
   neg_sum = neg_sum/(float)(nx*(ny_end-ny_str));

   if(truncate_zero_slip)
      fprintf(stderr,"ratio (negative slip)/(positive slip)= %f\n",neg_sum/(sum+neg_sum));

   taper_slip_all(cslip,nx,ny_in,&side_taper,&bot_taper,&top_taper);

   if(slip_water_level > 0.0)
      {
      slipmin = sum*slip_water_level;
      for(j=ny_str;j<ny_end;j++)
         {
         for(i=0;i<nx;i++)
            {
            ip = i + j*nx;
            if(cslip[ip].re < slipmin)
               cslip[ip].re = slipmin;
            }
	    /*if(cslip[ip][0] < slipmin)
		cslip[ip][0] = slipmin;
	    }*/
         }
      }

/* check moment and scale slip */
   savg = target_savg;
   //fprintf(stderr,"savg=%f\n", savg);
   scale_slip_fftw(psrc,cslip,nx,ny_in,ny_str,&dx,&dy,&dtop,&avgdip,&mom,&vmod,&savg,&smax);
   fprintf(stderr,"mom= %13.5e avgslip= %.0f maxslip= %.0f\n",mom,savg,smax);

   fprintf(stderr,"orig_sigma= %f ... ",sigma);

/* recalculate just to check, plus normalize by savg */
   sigma = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         sigma = sigma + (psrc[ip].slip-savg)*(psrc[ip].slip-savg);
         }
      }
   sigma = sqrt(sigma/((ny_end-ny_str)*nx))/savg;

   fprintf(stderr,"new_sigma= %f\n",sigma);

   fprintf(stderr,"seed=%d\n", seed);

   savg = 0.0;
   snum = 0.0;
   sden = 0.0;
   for(j=0;j<ny_in*nx;j++)
      {
      savg = savg + psrc[j].slip;

      /* this uses slip-weighted averaging (does NOT include depth variation) */

      sabs = sqrt(psrc[j].slip*psrc[j].slip);
      snum = snum + sabs*exp(0.5*log(sabs));
      sden = sden + sabs;
      //fprintf(stderr,"j=%d, snum=%f, den=%f\n", j, snum, sden);
      }

   savg = savg/(float)(nx*ny_in);

   if(stfparams.rt_scalefac > 0)
      {
      stfparams.rt_scalefac = snum/sden;
      fprintf(stderr,"snum=%f, den=%f\n", snum, sden);
      fprintf(stderr,"rt_scalefac= %f\n",stfparams.rt_scalefac);
      }

/* now do rake */

   for(j=0;j<ny_in*nx;j++)
      {
      crake[j].re = 1.0;
      crake[j].im = 0.0;
      }

   fft2d_fftw(crake,nx,ny_in,-1,&dx,&dy);
   kfilt(crake,nx,ny_in,&dkx_rk,&dky_rk,&xl,&yl,&seed,kmodel);
   fft2d_fftw(crake,nx,ny_in,1,&dkx_rk,&dky_rk);

   printf("A: seed = %ld\n", seed);

   //for(j=0;j<ny_in*nx;j++)
   //   sort_rake[j] = crake[j].re;

   /*test = 1;
   while(test)
      {
      test = 0;
      for(i=0;i<nx*ny_in-1;i++)
         {
         if(sort_rake[i] > sort_rake[i+1])
            {
            ravg = sort_rake[i];
            sort_rake[i] = sort_rake[i+1];
            sort_rake[i+1] = ravg;
            test = 1;
            }
         }
      }*/
   //qsort(sort_rake, nx*ny_in, sizeof(float), compareFloat);

   //ip = (int)(0.5*nx*ny_in - 0.5);
   //rmed = sort_rake[ip];

   rmin = 1.0e+15;
   rmax = -1.0e+15;
   ravg = 0.0;
   for(j=0;j<ny_in;j++)
      {
      for(i=0;i<nx;i++)
         {
         ip = i + j*nx;

         if(crake[ip].re < rmin)
            rmin = crake[ip].re;
         if(crake[ip].re > rmax)
            rmax = crake[ip].re;

         ravg = ravg + crake[ip].re;
         }
      }

   ravg = ravg/(float)(nx*ny_in);

   fprintf(stderr,"ravg= %13.5e rmed= %13.5e rmin= %13.5e rmax= %13.5e\n",ravg,rmed,rmin,rmax);

   if(rmax-ravg >= ravg-rmin)
      rfac = rand_rake_degs/(rmax-ravg);
   else
      rfac = rand_rake_degs/(ravg-rmin);

   if(set_rake > -900.0)
      {
      for(j=0;j<ny_in*nx;j++)
         psrc_rake[j] = set_rake;
      }

   for(j=0;j<ny_in*nx;j++)
      psrc[j].rak = (crake[j].re - ravg)*rfac + psrc_rake[j];

/* XXXX */
   stfparams.deep_risetimedep = deep_risetimedep_saved;
   xhypo = dhypo*sin(avgdip*rperd) + dtop + stfparams.deep_risetimedep_range;
   if(xhypo > stfparams.deep_risetimedep)
      stfparams.deep_risetimedep = xhypo;
/* XXXX */

   printf("B: seed = %ld\n", seed);

   /*
   if (js<doslip) {
	int iseg, idip, istk;
	int nseg = srf->srf_prect.nseg;
	int ntot = 0;
        struct srf_prectsegments *prseg_ptr = srf->srf_prect.prectseg;
	for(iseg=0; iseg<nseg; iseg++) {
		ntot += prseg_ptr[iseg].nstk;
	}
	int ioff = 0;
	printf("nseg: %d\n", nseg);
	for(iseg=0; iseg<nseg; iseg++) {
		for (idip=0;idip<prseg_ptr[iseg].ndip;idip++) {
			for(istk=0;istk<prseg_ptr[iseg].nstk; istk++) {
				int ip0 = istk + ioff + idip*ntot;
				float sabs = sqrt(psrc[ip0].slip*psrc[ip0].slip);
				if (sabs > MINSLIP) {
					if (stfparams.rt_scalefac > 0.0) {
						if(stfparams.rt_rand > 0.0) {
				                	_gaus_rand(&(stfparams.rt_rand),&fzero,&seed);
						}
					}
				}
			}
		}
		ioff += prseg_ptr[iseg].nstk;
	}
	//Now do the hypos
	for (ih=0; ih<nh; ih++) {
		for(j=0; j<ny_in; j++) {
                        for(i=0; i<nx; i++) {
                               if(tsfac_rand > 0.0) {
                                     _gaus_rand(&tsfac_rand,&fzero,&seed);
                                }

                                if (sd_rand > 0.0) {
                                        _gaus_rand(&fone,&fzero,&seed);
                                        _gaus_rand(&fone,&fzero,&seed);
                                }
                        }
               }
	}
	continue;
   }*/

   if(write_srf) {
      _load_slip_srf_dd2(srf,&stfparams,psrc,&seed);
   }

   for(ih=0;ih<nh;ih++)    /* loop over hypocenter realizations */
      {
      /*if (((doslip >= 0) && (doslip != js)) ||
          ((dohypo >= 0) && (dohypo != ih))) {
        continue;
      }*/

	/*if (ih!=dohypo) {
	       for(j=0; j<ny_in; j++) {
			for(i=0; i<nx; i++) {
		               if(tsfac_rand > 0.0) {
        		             _gaus_rand(&tsfac_rand,&fzero,&seed);
                		}
	
        		        if (sd_rand > 0.0) {
        		                _gaus_rand(&fone,&fzero,&seed);
                		        _gaus_rand(&fone,&fzero,&seed);
                		}
			}
		}
                continue;
	}*/


      if(random_hypo != 1 && calc_shypo == 1)
         shypo = sh0 + ih*shypo_step - 0.5*flen;

/* calculate rupture time */

      if((smax/savg) != (float)(1.0))
         sfac0 = 1.0/log(smax/savg);
      else
         tsfac = sfac0 = 0.0;

      tsmin = 1.0e+15;
      float half_dy = 0.5*dy;
      float half_dx_minus_half_flen = 0.5*dx - 0.5*flen;
      float dep_const = sfac0*(tsfac_factor - 1.0)/(depmax1-depmin1);
      float dep_alpha = sfac0 + dep_const * depmax1;
      for(j=0;j<ny_in;j++)
         {
         //yy = (j + 0.5)*dy;
         yy = j*dy + half_dy;
         for(i=0;i<nx;i++)
            {
            //xx = (i+0.5)*dx - 0.5*flen;
            xx = i*dx + half_dx_minus_half_flen;
	    ip = i + j*nx;

            if(psrc[ip].dep >= depmax1)
               sf = sfac0;
            else if(psrc[ip].dep < depmax1 && psrc[ip].dep > depmin1) {
               //sf = sfac0*(1.0 + (tsfac_factor - 1.0)*(depmax1-psrc[ip].dep)/(depmax1-depmin1));
               sf = dep_alpha - dep_const*psrc[ip].dep;
            } else
               sf = sfac0*tsfac_factor;

	    shypo_mseg = shypo;
	    if(seg_delay == 1)
	       {
	       delh = 0.0;
	       hx = xx - shypo;

               for(ig=0;ig<nseg_bounds;ig++)
                  {
                  hg = xseg[ig]-shypo;
	          gwid2 = 0.5*gwid[ig];

                  if(hx >= 0.0 && hg >= 0.0 && shypo <= (xseg[ig]-gwid2))
                     {
                     if(hx > (hg - gwid2) && hx < (hg + gwid2))
                        delh = delh + rvfac_seg[ig]*(hx - (hg - gwid2));
                     else if(hx >= (hg + gwid2))
                        delh = delh + rvfac_seg[ig]*gwid[ig];
                     }

                  else if(hx < 0.0 && hg < 0.0 && shypo >= (xseg[ig]+gwid2))
                     {
                     if(hx < (hg + gwid2) && hx > (hg - gwid2))
                        delh = delh + rvfac_seg[ig]*(hx - (hg + gwid2));
                     else if(hx <= (hg - gwid2))
                        delh = delh - rvfac_seg[ig]*gwid[ig];
                     }

                  /*
                  else if(shypo >= (xseg[ig]-gwid2) && shypo <= (xseg[ig]+gwid2))
                     {
                     if(xx >= (xseg[ig] - gwid2) && xx <= (xseg[ig] + gwid2))
                        delh = delh + rvfac_seg[ig]*hx;
                     else if(xx < (xseg[ig]-gwid2))
                        delh = delh - rvfac_seg[ig]*(shypo-(xseg[ig]-gwid2));
                     else if(xx > (xseg[ig]+gwid2))
                        delh = delh + rvfac_seg[ig]*((xseg[ig]+gwid2)-shypo);
                     }
                  */
		  }

	       shypo_mseg = xx - (hx + delh);
	       }

            get_rupt(&rvmod,&htol,&dhypo,&yy,&shypo_mseg,&xx,&rayp,&rupt_rad,&rt);

/* 20131119 RWG:

   Add random perturbations to tsfac so that it is not 1:1 correlated with slip.
   Perturbations are log normal with ln(sigma)=tsfac_rand*tsfac.  Default for
   tsfac_rand=0.2.
   
*/

	    /*if (js!=doslip || ih!=dohypo) {
	       if(tsfac_rand > 0.0) {
          	     _gaus_rand(&tsfac_rand,&fzero,&seed);
		}

		if (sd_rand > 0.0) {
			_gaus_rand(&fone,&fzero,&seed);
			_gaus_rand(&fone,&fzero,&seed);
		}
		continue;
	     }*/


            if(tsfac_rand > 0.0)
               sf = sf*exp(_gaus_rand(&tsfac_rand,&fzero,&seed));

            if(psrc[ip].slip > 0.25*savg*savg/smax)
               psrc[ip].rupt = rt + sf*tsfac*log(psrc[ip].slip/savg);
            else
               psrc[ip].rupt = rt - sf*tsfac;

	    if(sd_rand > 0.0)
	       {
               psrc[ip].stk = psrc[ip].stk + sd_rand*_gaus_rand(&fone,&fzero,&seed);
               psrc[ip].dip = psrc[ip].dip + sd_rand*_gaus_rand(&fone,&fzero,&seed);
	       }

	    if(psrc[ip].rupt < tsmin)
	       tsmin = psrc[ip].rupt;
            }
         }

      printf("C: seed = %ld\n", seed);

      /* adjust to start at rt=0.0 */
      for(j=0;j<nx*ny_in;j++) {
         psrc[j].rupt = psrc[j].rupt - tsmin + rup_delay;
      }

      if (dohypo!=ih || doslip!=js) {
	continue;
      }

/* RWG 2014-02-20 randomized hypocenter */
      if(write_srf)
	 {
         _load_rupt_srf(srf,psrc,&shypo,&dhypo);

         if(strcmp(outfile,"stdout") == 0)
            sprintf(str,"stdout");

         else if(read_erf == 1)
	    {
            if(random_hypo == 1)
               sprintf(str,"%s-r%.6d.srf",outfile,js);

            else             /* old way */
               sprintf(str,"%s-s%.4d-h%.4d",outfile,js,ih);
	    }

         else
            sprintf(str,"%s.srf",outfile);

         write2srf(srf,str,outbin);
	 }

      else if(gslip.np > 0 && write_gsf)
         {
         if(strcmp(outfile,"stdout") == 0)
            sprintf(str,"stdout");
         else
            sprintf(str,"%s-s%.4d-h%.4d.gsf",outfile,js,ih);

         write2gsf(&gslip,psrc,infile,str);
	 }
      }

   //free_srf_stf(&srf);
   }
 
 //Free things 
 free(xseg);
 free(rvfac_seg);
 free(gwid);
 free(cslip);
 free(crake);
 free(sort_rake);
 free(psrc_rake);
 free(psrc);

 free(vmod.vp);
 free(vmod.vs);
 free(vmod.den);
 free(vmod.th);
 free(vmod.dep);
 free(vmod.mu);
 free(vmod.invb2);

 free(rvmod.vs);
 free(rvmod.th);
 free(rvmod.invb2);

 return(0);
}

void write2srf(struct standrupformat *srf,char *str,int outbin)
{
_write_srf(srf,str,outbin);
}
