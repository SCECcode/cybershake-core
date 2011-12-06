#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "genslip.h"
#include "getpar.h"
#include "ruptime.h"

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
*/

int genslip(int ac,char **av, rg_stats_t *stats)
{
FILE *fpr, *fpw;
struct complex *cslip, *crake;
float *aspec;
float flen, fwid, dx, dy, sval;
float dkx, dky, kx, ky, xx, yy, zz;
float cosA, sinA, cosD, sinD, dd, sn, se;
float sum, sigma, neg_sum;
int nx, ny;
int i, j, k, ip, nt6, it;
char infile[1024], slipfile[1024], specfile[1024], ruptfile[1024], str[1024];
char init_slip_file[1024], outfile[1024], logfile[1024];
FILE *lfile;

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

long seed = 0;
int kmodel = MAI_FLAG;   /* default is mai */
int circular_average = DEFAULT_CIRCULAR_AVERAGE;
int modified_corners = DEFAULT_MODIFIED_CORNERS;

float kx_corner, ky_corner, xmag_exponent, ymag_exponent;

float mag_area_Acoef = -1.0;
float mag_area_Bcoef = -1.0;
int use_median_mag = 0;

int calc_shypo = 1;
float shypo = -1.0e+15;
float dhypo = -1.0e+15;
float avgstk, rake, rt, tsmin;
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
float savg = -1.0;

float depmin, depmax;
float shal_vrup_dep = DEFAULT_DEPTH_SCALING_LEVEL;
float shal_vrup_deprange = DEFAULT_DEPTH_SCALING_RANGE;

float tsfac = -1.0e+15;;
float tsfac_coef = DEFAULT_TSFAC_COEF;
float tsfac_factor = DEFAULT_TSFAC_FACTOR;
float tsfac_dep = DEFAULT_DEPTH_SCALING_LEVEL;
float tsfac_deprange = DEFAULT_DEPTH_SCALING_RANGE;
float tsfac_rand = -1.0;
float sd_rand = -1.0;
float fzero = 0.0;
float fone = 1.0;

float sigfac;
float slip_sigma = DEFAULT_SLIP_SIGMA;

float dtop, avgdip, mom, mag_med;
struct velmodel vmod;
char velfile[1024];

int read_erf = 1;

struct standrupformat srf;
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
struct stfpar stfparams;

int outbin = 0;
int writeout = 1;
int doslip = -1;
int dohypo = -1;

slipfile[0] = '\0';
specfile[0] = '\0';
ruptfile[0] = '\0';

velfile[0] = '\0';

init_slip_file[0] = '\0';

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(logfile, "log");

sprintf(srf.version,"1.0");

 rvfac_seg = NULL;
 xseg = NULL;
 gwid = NULL;
 aspec = NULL;

stfparams.dt = DEFAULT_DT;
stfparams.nt = NTMAX;
stfparams.trise = -1.0;
stfparams.rt_scalefac = DEFAULT_RT_SCALEFAC;
stfparams.risetimedep = DEFAULT_DEPTH_SCALING_LEVEL;
stfparams.risetimedep_range = DEFAULT_DEPTH_SCALING_RANGE;
stfparams.risetimefac = 2.0;
sprintf(stfparams.stype,"ucsb");  /* default is UCSB sincos */

setpar(ac,av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("logfile","s",logfile);
getpar("outbin","d",&outbin);
getpar("writeout","d",&writeout);
getpar("doslip","d",&doslip);
getpar("dohypo","d",&dohypo);

getpar("velfile","s",velfile);

getpar("read_erf","d",&read_erf);
getpar("read_gsf","d",&read_gsf);

getpar("write_srf","d",&write_srf);

lfile = fopen(logfile, "w");
if (lfile == NULL) {
  fprintf(stderr, "Failed to open logfile %s\n", logfile);
  lfile = stderr;
}

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
getpar("risetimefac","f",&stfparams.risetimefac);
getpar("risetimedep","f",&stfparams.risetimedep);
getpar("risetimedep_range","f",&stfparams.risetimedep_range);
getpar("rt_scalefac","f",&stfparams.rt_scalefac);
getpar("stype","s",stfparams.stype);

getpar("shypo_step","f",&shypo_step);
getpar("shypo_min_off","f",&shypo_min_off);
getpar("dhypo_frac","f",&dhypo_frac);
getpar("slips_to_hypos","d",&slips_to_hypos);
getpar("ns","d",&ns);
getpar("nh","d",&nh);

getpar("shypo","f",&shypo);
getpar("dhypo","f",&dhypo);

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

getpar("savg","f",&savg);

getpar("slipfile","s",slipfile);
getpar("specfile","s",specfile);
getpar("ruptfile","s",ruptfile);

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

getpar("slip_sigma","f",&slip_sigma);

getpar("init_slip_file","s",init_slip_file);

getpar("seg_delay","d",&seg_delay);
if(seg_delay == 1)
   {
   mstpar("nseg_bounds","d",&nseg_bounds);

   xseg = (float *)check_malloc(nseg_bounds*sizeof(float));
   rvfac_seg = (float *)check_malloc(nseg_bounds*sizeof(float));
   gwid = (float *)check_malloc(nseg_bounds*sizeof(float));

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
   psrc = read_ruppars(infile,psrc,&mag,&nx,&ny_in,&dx,&dy,&dtop,&avgstk,&avgdip,&elon,&elat);
else if(read_gsf == 1)
   psrc = read_gsfpars(infile,psrc,&gslip,&dx,&dy,&dtop,&avgdip);
else
   psrc = set_ruppars(psrc,&mag,&nx,&ny_in,&dx,&dy,&dtop,&avgstk,&avgdip,&rake,&elon,&elat);

flen = nx*dx;
fwid = ny_in*dy;

if(shypo > -1.0e+14)
   calc_shypo = 0;

if(nh < 0)
   nh = (int)((flen-2.0*shypo_min_off)/shypo_step) + 1;
if(ns < 0)
   {
   ns = slips_to_hypos*nh;
   if(nx == 1 && ny_in == 1)
      ns = 1;
   }

sh0 = 0.5*(flen - (nh-1)*shypo_step);
if(dhypo < -1.0e+14)
   dhypo = dhypo_frac*fwid;

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

/* PES collect statistics */
stats->numslip = ns;
stats->numhypo = nh;

fprintf(lfile,"mag= %.2f median mag= %.2f nslip= %d nhypo= %d\n",mag,mag_med,ns,nh);
fprintf(lfile,"nx= %d ny= %d dx= %10.4f dy= %10.4f\n",nx,ny_in,dx,dy);

init_plane_srf(&srf,&gslip,&elon,&elat,nx,ny_in,&flen,&fwid,&dx,&dy,&avgstk,&avgdip,&dtop,&shypo,&dhypo);

if(velfile[0] != '\0')
   read_velmodel(velfile,&vmod);
else
   default_velmodel(&vmod);

/*
conv2vrup(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup);
*/
depmin = shal_vrup_dep - shal_vrup_deprange;
depmax = shal_vrup_dep + shal_vrup_deprange;
conv2vrup_dd(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup,&depmin,&depmax);

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

/* rise time from somerville */
if(stfparams.trise < 0.0)
   stfparams.trise = 1.6e-09*exp(log(mom)/3.0);

cslip = (struct complex *) check_malloc (nx*ny*sizeof(struct complex));

dkx_rk = 1.0/(nx*dx);
dky_rk = 1.0/(ny_in*dy);
crake = (struct complex *) check_malloc (nx*ny_in*sizeof(struct complex));
sort_rake = (float *) check_malloc (nx*ny_in*sizeof(float));
psrc_rake = (float *) check_malloc (nx*ny_in*sizeof(float));

for(j=0;j<ny_in*nx;j++)
   psrc_rake[j] = psrc[j].rak;

if(specfile[0] != '\0')
   {
   aspec = (float *) check_malloc (nx*ny*sizeof(float));
   for(j=0;j<nx*ny;j++)
      aspec[j] = 0.0;
   }

for(js=0;js<ns;js++)    /* loop over slip realizations */
   {
   /*
   init_slip(cslip,nx,ny,&side_taper,&bot_taper);

   for(j=ny/2-30;j<ny/2+30;j++)
      {
      fprintf(lfile,"j=%d\n",j);
      for(i=0;i<nx/2;i++)
         {
	 cslip[i+j*nx].re = 0.01*cslip[i+j*nx].re;
	 }
      }
   */

   init_slip_IO(cslip,nx,ny,&dx,&dy,flip_at_surface,init_slip_file);

   fft2d(cslip,nx,ny,-1,&dx,&dy);
   if(kmodel < 100)
      kfilt(cslip,nx,ny,&dkx,&dky,&xl,&yl,&seed,kmodel);
   fft2d(cslip,nx,ny,1,&dkx,&dky);

   taper_slip_all(cslip,nx,ny_in,&side_taper,&bot_taper,&top_taper);

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
            }
         }
      }

   fprintf(lfile,"orig_sigma= %f ... ",sigma);

/* recalculate just to check, plus normalize by sum */
   sigma = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         sigma = sigma + (cslip[ip].re-sum)*(cslip[ip].re-sum);
         }
      }
   sigma = sqrt(sigma/((ny_end-ny_str)*nx))/sum;

   fprintf(lfile,"new_sigma= %f\n",sigma);

   neg_sum = 0.0;
   for(j=ny_str;j<ny_end;j++)
      {
      for(i=0;i<nx;i++)
         {
	 ip = i + j*nx;
         if(cslip[ip].re < 0.0 && truncate_zero_slip)
	    {
	    neg_sum = neg_sum - cslip[ip].re;
            cslip[ip].re = 0.0;
	    /*
	    fprintf(lfile,"slip= %13.5e\n",cslip[ip].re);
	    */
	    }
	 }
      }

   neg_sum = neg_sum/(float)(nx*(ny_end-ny_str));

   if(truncate_zero_slip)
      fprintf(lfile,"ratio (negative slip)/(positive slip)= %f\n",neg_sum/(sum+neg_sum));

/* check moment and scale slip */

   scale_slip(psrc,cslip,nx,ny_in,ny_str,&dx,&dy,&dtop,&avgdip,&mom,&vmod,&savg,&smax);
   fprintf(lfile,"mom= %13.5e avgslip= %.0f maxslip= %.0f sigma= %.2f\n",mom,savg,smax,savg*sigma);

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
      }

   savg = savg/(float)(nx*ny_in);

   if(stfparams.rt_scalefac > 0)
      {
      stfparams.rt_scalefac = snum/sden;
      fprintf(lfile,"rt_scalefac= %f\n",stfparams.rt_scalefac);
      }

/* now do rake */

   for(j=0;j<ny_in*nx;j++)
      {
      crake[j].re = 1.0;
      crake[j].im = 0.0;
      }

   fft2d(crake,nx,ny_in,-1,&dx,&dy);
   kfilt(crake,nx,ny_in,&dkx_rk,&dky_rk,&xl,&yl,&seed,kmodel);
   fft2d(crake,nx,ny_in,1,&dkx_rk,&dky_rk);

   for(j=0;j<ny_in*nx;j++)
      sort_rake[j] = crake[j].re;

   test = 1;
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
      }

   ip = (int)(0.5*nx*ny_in - 0.5);
   rmed = sort_rake[ip];

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

   fprintf(lfile,"ravg= %13.5e rmed= %13.5e rmin= %13.5e rmax= %13.5e\n",ravg,rmed,rmin,rmax);

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

   if(slipfile[0] != '\0')
      {
      sprintf(str,"%s-s%.4d",slipfile,js);
      write_field(str,psrc,"slip",nx,ny_in,&dx,&dy);
      }

   if(specfile[0] != '\0')
      {
      sprintf(str,"%s-s%.4d",specfile,js);
      write_spec(str,aspec,cslip,nx,ny,&dx,&dy,&dkx,&dky,&xl,&yl,kmodel);
      }

   load_slip_srf(&srf,&stfparams,psrc);

   for(ih=0;ih<nh;ih++)    /* loop over hypocenter realizations */
      {
      if(calc_shypo == 1)
         shypo = sh0 + ih*shypo_step - 0.5*flen;

/* calculate rupture time */

      if((smax/savg) != (float)(1.0))
         sfac0 = 1.0/log(smax/savg);
      else
         tsfac = sfac0 = 0.0;

      tsmin = 1.0e+15;
      for(j=0;j<ny_in;j++)
         {
         yy = (j + 0.5)*dy;
         for(i=0;i<nx;i++)
            {
            xx = (i+0.5)*dx - 0.5*flen;
	    ip = i + j*nx;

            if(psrc[ip].dep >= depmax)
               sf = sfac0;
            else if(psrc[ip].dep < depmax && psrc[ip].dep > depmin)
               sf = sfac0*(1.0 + (tsfac_factor - 1.0)*(depmax-psrc[ip].dep)/(depmax-depmin));
            else
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

                  if(hx >= 0.0 && hg >= 0.0)
                     {
                     if(hx > (hg - gwid2) && hx < (hg + gwid2))
                        delh = delh + rvfac_seg[ig]*(hx - (hg - gwid2));
                     else if(hx >= (hg + gwid2))
                        delh = delh + rvfac_seg[ig]*gwid[ig];
                     }

                  else if(hx < 0.0 && hg < 0.0)
                     {
                     if(hx < (hg + gwid2) && hx > (hg - gwid2))
                        delh = delh + rvfac_seg[ig]*(hx - (hg + gwid2));
                     else if(hx <= (hg - gwid2))
                        delh = delh - rvfac_seg[ig]*gwid[ig];
                     }
		  }

	       shypo_mseg = xx - (hx + delh);
	       }

            get_rupt(&rvmod,&htol,&dhypo,&yy,&shypo_mseg,&xx,&rayp,&rupt_rad,&rt);

	    sval = 0.0;
            if(tsfac_rand > 0.0)
               sval = gaus_rand(&tsfac_rand,&fzero,&seed);

	    sval = 0.0;
            if(psrc[ip].slip > 0.25*savg*savg/smax)
               psrc[ip].rupt = rt + sf*tsfac*(1.0 + sval)*log(psrc[ip].slip/savg);
            else
               psrc[ip].rupt = rt - sf*tsfac;

	    if(sd_rand > 0.0)
	       {
               psrc[ip].stk = psrc[ip].stk + sd_rand*gaus_rand(&fone,&fzero,&seed);
               psrc[ip].dip = psrc[ip].dip + sd_rand*gaus_rand(&fone,&fzero,&seed);
	       }

/*
	    sigma = 0.15;
            sval = psrc[ip].slip*(1.0 + gaus_rand(&sigma,&fzero,&seed));

            if(sval > 0.25*savg*savg/smax)
               psrc[ip].rupt = rt + sf*tsfac*log(sval/savg);
            else
               psrc[ip].rupt = rt - sf*tsfac;
*/

	    if(psrc[ip].rupt < tsmin)
	       tsmin = psrc[ip].rupt;
            }
         }

      /* adjust to start at rt=0.0 */
      for(j=0;j<nx*ny_in;j++)
         psrc[j].rupt = psrc[j].rupt - tsmin + rup_delay;

      if(ruptfile[0] != '\0')
         {
         sprintf(str,"%s-s%.4d-h%.4d",ruptfile,js,ih);
         write_field(str,psrc,"rupt",nx,ny_in,&dx,&dy);
         }

      load_rupt_srf(&srf,psrc,&shypo,&dhypo);

      if(strcmp(outfile,"stdout") == 0)
         sprintf(str,"stdout");
      else
         sprintf(str,"%s-s%.4d-h%.4d",outfile,js,ih);

      if (((doslip < 0) || (doslip == js)) &&
	  ((dohypo < 0) || (dohypo == ih))) {

      if((write_srf) && (writeout))
         write2srf(&srf,str,outbin);

      }

      if(strcmp(outfile,"stdout") == 0)
         sprintf(str,"stdout");
      else
         sprintf(str,"%s-s%.4d-h%.4d.gsf",outfile,js,ih);

      if(gslip.np > 0 && write_gsf && writeout)
         write2gsf(&gslip,psrc,infile,str);
      }

   free_srf_stf(&srf);
   }

if(specfile[0] != '\0')
   {
   sprintf(str,"%s-s_avg",specfile);
   write_avgspec(str,aspec,ns,nx,ny,&dkx,&dky);
   }

 if (lfile != stderr) {
   fclose(lfile);
 }
 return(0);
}

void write2srf(struct standrupformat *srf,char *str,int outbin)
{
write_srf(srf,str,outbin);
}
