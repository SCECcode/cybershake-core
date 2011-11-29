#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

/*

   03/06/2000: version 2.1

   * Added option to read input fault description from Generic Slip Format (GSF)
     file.

   * Added option to specify "Brune" type STF.

*/

main(int ac,char **av)
{
FILE *fpr, *fpw;
struct complex *cslip;
float *aspec;
float flen, fwid, dx, dy, amp0, amp, lamp;
float dkx, dky, kx, ky, xx, yy, zz;
float cosA, sinA, cosD, sinD, dd, sn, se;
int nx, ny;
int i, j, k, ip, nt6, it;
char infile[256], slipfile[256], specfile[256], ruptfile[256], str[512];
char outfile[256];

float bigM;
float xl, yl;

int flip_at_surface = -1;
int stretch_kcorner = 0;
int ny_in, ny_str, ny_end;

float pi = 3.14159265;
double rperd = 0.017453293;

float mag;
float side_taper = DEFAULT_SIDE_TAP;
float bot_taper = DEFAULT_BOT_TAP;

long seed = 0;
int kmodel = 1;   /* default is somerville */
int circular_average = 1;   /* default is cicular average of correlation lengths (mai only) */
int modified_corners = 1;   /* default is modified xL and yL (both mai & somerville) */

int calc_shypo = 1;
float shypo = -1.0e+15;
float dhypo = -1.0e+15;
float avgstk, rake, rt, tsmin;
struct velmodel rvmod;
double rayp, rupt_rad;

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

float smax, sf;
float savg = -1.0;
float tsfac = DEFAULT_TSFAC;

float dtop, avgdip, mom, mag_med;
struct velmodel vmod;
char velfile[128];

int read_erf = 1;

struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

int read_gsf = 0;
int write_gsf = 0;

struct generic_slip gslip;

float *stf, elon, elat;
int ig, ntot;

struct pointsource *psrc;
struct stfpar stfparams;

int outbin = 0;

slipfile[0] = '\0';
specfile[0] = '\0';
ruptfile[0] = '\0';

velfile[0] = '\0';

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

sprintf(srf.version,"1.0");

stfparams.dt = DEFAULT_DT;
stfparams.nt = NTMAX;
stfparams.trise = -1.0;
sprintf(stfparams.stype,"urs");  /* default is URS 2tri STF */

setpar(ac,av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);

getpar("velfile","s",velfile);

getpar("read_erf","d",&read_erf);
getpar("read_gsf","d",&read_gsf);

if(read_gsf == 1)
   {
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
getpar("tsfac","f",&tsfac);

getpar("savg","f",&savg);

getpar("slipfile","s",slipfile);
getpar("specfile","s",specfile);
getpar("ruptfile","s",ruptfile);

getpar("kmodel","d",&kmodel);
getpar("modified_corners","d",&modified_corners);
getpar("circular_average","d",&circular_average);
getpar("flip_at_surface","d",&flip_at_surface);
getpar("stretch_kcorner","d",&stretch_kcorner);
getpar("seed","d",&seed);
getpar("side_taper","f",&side_taper);
getpar("bot_taper","f",&bot_taper);

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

fprintf(stderr,"mag= %.2f median mag= %.2f nslip= %d nhypo= %d\n",mag,mag_med,ns,nh);
fprintf(stderr,"nx= %d ny= %d dx= %10.4f dy= %10.4f\n",nx,ny_in,dx,dy);

init_plane_srf(&srf,&gslip,&elon,&elat,nx,ny_in,&flen,&fwid,&dx,&dy,&avgstk,&avgdip,&dtop,&shypo,&dhypo);

if(velfile[0] != '\0')
   read_velmodel(velfile,&vmod);
else
   default_velmodel(&vmod);

conv2vrup(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup);

if(kmodel != MAI_FLAG)
   kmodel = SOMERVILLE_FLAG;

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
   xl = exp(bigM*(0.5*mag_med - 2.50 + 0.79818));
   yl = exp(bigM*(0.3333*mag_med - 1.50 + 0.79818));

   if(circular_average)
      yl = xl;
   }

if(kmodel == SOMERVILLE_FLAG) /* somerville scaling */
   {
   xl = exp(bigM*(0.5*mag_med - 1.72));
   yl = exp(bigM*(0.5*mag_med - 1.93));

   if(circular_average)
      {
      xl = exp(bigM*(0.5*mag_med - 1.825));
      yl = xl;
      }
   }

if(modified_corners)
   {
   xl = exp(bigM*(0.5*mag_med - 2.00));
   yl = exp(bigM*(0.5*mag_med - 2.00));
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
   stfparams.trise = 1.8e-09*exp(log(mom)/3.0);

cslip = (struct complex *) check_malloc (nx*ny*sizeof(struct complex));

if(specfile[0] != '\0')
   {
   aspec = (float *) check_malloc (nx*ny*sizeof(float));
   for(j=0;j<nx*ny;j++)
      aspec[j] = 0.0;
   }

for(js=0;js<ns;js++)    /* loop over slip realizations */
   {
   init_slip(cslip,nx,ny,&side_taper,&bot_taper);

   fft2d(cslip,nx,ny,-1,&dx,&dy);
   kfilt(cslip,nx,ny,&dkx,&dky,&xl,&yl,&seed,kmodel);
   fft2d(cslip,nx,ny,1,&dkx,&dky);

   taper_slip(cslip,nx,ny,&side_taper,&bot_taper);

/*
   truncate any negative slip values => should double check that spectra is
   not changed too much
*/

   for(j=0;j<ny*nx;j++)
      {
      if(cslip[j].re < 0.0)
	 {
         cslip[j].re = 0.0;
	 /*
	 */
	 }
      }

   ny_str = 0;
   ny_end = ny;
   if(flip_at_surface == 1)
      ny_str = ny_in;

/* check moment and scale slip */

   scale_slip(psrc,cslip,nx,ny_in,ny_str,&dx,&dy,&dtop,&avgdip,&mom,&vmod,&savg,&smax);
   fprintf(stderr,"mom= %13.5e avgslip= %.0f maxslip= %.0f\n",mom,savg,smax);

   savg = 0.0;
   for(j=0;j<ny_in*nx;j++)
      savg = savg + psrc[j].slip;

   savg = savg/(float)(nx*ny_in);
   for(j=0;j<ny_in*nx;j++)
      {
      if(psrc[j].slip < 0.0)
	 {
	 /*
	 fprintf(stderr,"NEG: %12.5f\n",-100.0*psrc[j].slip/savg);
	 */
	 }
      }

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

      if((smax-savg) != (float)(0.0))
         sf = 1.0/(smax-savg);
      else
         tsfac = sf = 0.0;

      tsmin = 1.0e+15;
      for(j=0;j<ny_in;j++)
         {
         yy = (j + 0.5)*dy;
         for(i=0;i<nx;i++)
            {
            xx = (i+0.5)*dx - 0.5*flen;

            get_rupt(&rvmod,&htol,&dhypo,&yy,&shypo,&xx,&rayp,&rupt_rad,&rt);
            psrc[i+j*nx].rupt = rt + sf*tsfac*(psrc[i + j*nx].slip - savg);

	    if(psrc[i+j*nx].rupt < tsmin)
	       tsmin = psrc[i+j*nx].rupt;
            }
         }

      /* adjust to start at rt=0.0 */
      for(j=0;j<nx*ny_in;j++)
            psrc[j].rupt = psrc[j].rupt - tsmin;

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

      write_srf(&srf,str,outbin);

      if(strcmp(outfile,"stdout") == 0)
         sprintf(str,"stdout");
      else
         sprintf(str,"%s-s%.4d-h%.4d.gsf",outfile,js,ih);

      if(gslip.np > 0 && write_gsf)
         write2gsf(&gslip,psrc,infile,str);
      }

   free_srf_stf(&srf);
   }

if(specfile[0] != '\0')
   {
   sprintf(str,"%s-s_avg",specfile);
   write_avgspec(str,aspec,ns,nx,ny,&dkx,&dky);
   }
}
