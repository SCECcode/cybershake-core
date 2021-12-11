#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void *check_malloc(size_t len)
{
void *ptr;

ptr = (void *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
/*double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}*/

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
/*double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}*/

double gaussian_rand(float *sigma,float *mean,long *seed)
{
double r = 0.0;
double six = 6.0;
double one = 1.0;
double half = 0.5;
int i;

for(i=0;i<12;i++)
   r = r + half*(one + sfrand(seed));

return((double)((r - six)*(*sigma) + *mean));
}

/*zapit(s,n)
float *s;
int n;
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}*/

void set_ne(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*sn = ((*slat) - (*elat))*kperd_n;
*se = ((*slon) - (*elon))*kperd_e;
}

void set_ll(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*slat = (*sn)/kperd_n + *elat;
*slon = (*se)/kperd_e + *elon;
}

void swap_in_place(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[3];
   cbuf[3] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[2];
   cbuf[2] = cv;

   cbuf = cbuf + 4;
   }
}

int write_xyz(char *file,struct standrupformat *srf,char *type,int ig,int calc_xy,float *xoff,float *yoff,int tsstr,int tsend,int tsinc,float *svmin,float *slipmin,int keepsign,int dump_slip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, outval, slipval;
float sv, svmax, *stf1, *stf2, *stf3;
float ss, dd, dfac;
int it, maxnt, ts, itdel;
int nx, ny;
int nt1, nt2, nt3;
int i, j, k, npskip;
int igst, ignd, igp;
int nd, place;
char str[512], frmt[32];

int isgn = 1;

float rperd = 0.017453293;

if(strcmp(srf[0].type,"PLANE") != 0 || ig >= srf[0].srf_prect.nseg)
   return(-99);

if(ig < 0)
   {
   igst = 0;
   ignd = srf[0].srf_prect.nseg;
   }
else
   {
   igst = ig;
   ignd = igst + 1;
   }

if(strncmp(type,"rupture",7) == 0)
   {
   if(tsstr < 0)
      tsstr = 0;
   if(tsend < 0)
      {
      apval_ptr = srf[0].srf_apnts.apntvals;
      for(k=0;k<srf[0].srf_apnts.np;k++)
         {
	 itdel = (int)(apval_ptr[k].tinit/apval_ptr[k].dt + 0.5);

	 if(tsend < itdel + apval_ptr[k].nt1)
            tsend = itdel + apval_ptr[k].nt1;
	 if(tsend < itdel + apval_ptr[k].nt2)
            tsend = itdel + apval_ptr[k].nt2;
	 if(tsend < itdel + apval_ptr[k].nt3)
            tsend = itdel + apval_ptr[k].nt3;
	 }
      }

   nd = 0;
   place = 1;
   while((tsend-tsstr)/tsinc > (place-1))
      {
      place = place*10;
      nd++;
      }

   sprintf(frmt,"%s%%.%dd",file,nd);
   }
else
   {
   tsstr = 0;
   tsend = 1;
   }

for(ts=tsstr;ts<tsend;ts=ts+tsinc)
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      {
      if(strncmp(type,"rupture",7) == 0)
         {
         sprintf(str,frmt,(ts-tsstr)/tsinc);
         fpw = fopfile(str,"w");
         }
      else
         fpw = fopfile(file,"w");
      }

   for(igp=igst;igp<ignd;igp++)
      {
      prseg_ptr = srf[0].srf_prect.prectseg;

      elon = prseg_ptr[igp].elon;
      elat = prseg_ptr[igp].elat;
      nx = prseg_ptr[igp].nstk;
      ny = prseg_ptr[igp].ndip;
      len = prseg_ptr[igp].flen;
      wid = prseg_ptr[igp].fwid;
      stk = prseg_ptr[igp].stk;
      dip = prseg_ptr[igp].dip;
      dtop = prseg_ptr[igp].dtop;
      shypo = prseg_ptr[igp].shyp;
      dhypo = prseg_ptr[igp].dhyp;

      len2 = 0.5*len;
      dx = len/(float)(nx);
      dy = wid/(float)(ny);

      npskip = 0;
      for(i=0;i<igp;i++)
         npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

      apval_ptr = srf[0].srf_apnts.apntvals + npskip;

      for(j=0;j<ny;j++)
         {
         yy = (j+0.5)*dy;
         for(i=0;i<nx;i++)
            {
            xx = (i+0.5)*dx;
	    k = i + j*nx;

            if(calc_xy)
               {
               set_ne(&(prseg_ptr[igp].elon),&(prseg_ptr[igp].elat),&(apval_ptr[k].lon),&(apval_ptr[k].lat),&sn,&se);

               arg = (apval_ptr[k].stk)*rperd;
               xx = se*sin(arg) + sn*cos(arg) + 0.5*prseg_ptr[igp].flen;

               arg = (apval_ptr[k].dip)*rperd;
               yy = ((apval_ptr[k].dep)-(prseg_ptr[igp].dtop))/sin(arg);
               }

	    xx = xx + *xoff;
	    yy = yy + *yoff;

            if(strncmp(type,"slip",4) == 0 || dump_slip)
               {
	       slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                apval_ptr[k].slip2*apval_ptr[k].slip2 +
	                apval_ptr[k].slip3*apval_ptr[k].slip3);

	       if(keepsign)
		  {
		  isgn = 1;

		  if(apval_ptr[k].slip1 < 0.0)
		     isgn = -1;
		  if((apval_ptr[k].slip2*apval_ptr[k].slip2) > (apval_ptr[k].slip1*apval_ptr[k].slip1) && apval_ptr[k].slip2 < 0.0)
		     isgn = -1;
		  if((apval_ptr[k].slip3*apval_ptr[k].slip3) > (apval_ptr[k].slip1*apval_ptr[k].slip1) && (apval_ptr[k].slip3*apval_ptr[k].slip3) > (apval_ptr[k].slip2*apval_ptr[k].slip2) && apval_ptr[k].slip3 < 0.0)
		     isgn = -1;

	          slip = isgn*slip;
		  }

               outval = slip;
               slipval = slip;
	       }

            if(strncmp(type,"velocity",8) == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

	       maxnt = apval_ptr[k].nt1;
	       if(apval_ptr[k].nt2 > maxnt)
	          maxnt = apval_ptr[k].nt2;
	       if(apval_ptr[k].nt3 > maxnt)
	          maxnt = apval_ptr[k].nt3;

	       svmax = 0.0;
               for(it=0;it<maxnt;it++)
	          {
	          sv = 0.0;
	          if(it < apval_ptr[k].nt1)
	             sv = sv + stf1[it]*stf1[it];
	          if(it < apval_ptr[k].nt2)
	             sv = sv + stf2[it]*stf2[it];
	          if(it < apval_ptr[k].nt3)
	             sv = sv + stf3[it]*stf3[it];

                  sv = sqrt(sv);
	          if(sv > svmax)
	             svmax = sv;
	          }

               outval = svmax;
	       }
            else if(strncmp(type,"rupture",7) == 0)
               {
               stf1 = apval_ptr[k].stf1;
               stf2 = apval_ptr[k].stf2;
               stf3 = apval_ptr[k].stf3;

	       nt1 = apval_ptr[k].nt1;
	       while(stf1[nt1-1] < *svmin && nt1 > 0)
	          nt1--;

	       nt2 = apval_ptr[k].nt2;
	       while(stf2[nt2-1] < *svmin && nt2 > 0)
	          nt2--;

	       nt3 = apval_ptr[k].nt3;
	       while(stf3[nt3-1] < *svmin && nt3 > 0)
	          nt3--;

	       itdel = ts - (int)(apval_ptr[k].tinit/apval_ptr[k].dt + 0.5);

	       if(itdel < 1)
	          sv = 0.0;
	       else if(itdel < nt1 || itdel < nt2 || itdel < nt3)
	          {
	          sv = 0.0;
	          if(itdel < nt1)
	             sv = sv + stf1[itdel]*stf1[itdel];
	          if(itdel < nt2)
	             sv = sv + stf2[itdel]*stf2[itdel];
	          if(itdel < nt3)
	             sv = sv + stf3[itdel]*stf3[itdel];

                  sv = sqrt(sv);
		  }
	       else
	          {
		  maxnt = nt1;
		  if(nt2 > maxnt)
		     maxnt = nt2;
		  if(nt3 > maxnt)
		     maxnt = nt3;

		  dfac = 1.0;
		  if((itdel - maxnt) < maxnt)
		     dfac = (float)(itdel-maxnt)/(float)(maxnt);

	          sv = -dfac*sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                         apval_ptr[k].slip2*apval_ptr[k].slip2 +
	                         apval_ptr[k].slip3*apval_ptr[k].slip3);

		  if(sv > -(*slipmin))
		     sv = -(*slipmin);
		  }

               outval = sv;
	       }

            else if(strncmp(type,"tinit",5) == 0)
               outval = apval_ptr[k].tinit;

            else if(strncmp(type,"trise",5) == 0)
	       {
               outval = apval_ptr[k].nt1*apval_ptr[k].dt;
	       if(apval_ptr[k].nt2*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt2*apval_ptr[k].dt;
	       if(apval_ptr[k].nt3*apval_ptr[k].dt > outval)
                  outval = apval_ptr[k].nt3*apval_ptr[k].dt;
	       }

            else if(strncmp(type,"strike",6) == 0)
               outval = apval_ptr[k].stk;

            else if(strncmp(type,"dip",3) == 0)
               outval = apval_ptr[k].dip;

            else if(strncmp(type,"rake",4) == 0)
	       {
	       ss = cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 - sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;
	       dd = sin(rperd*apval_ptr[k].rake)*apval_ptr[k].slip1 + cos(rperd*apval_ptr[k].rake)*apval_ptr[k].slip2;

	       outval = 90.0;
	       if(ss != 0.0)
	          outval = atan(dd/ss)/rperd;
               if(ss < 0.0)
	          outval = outval + 180.0;

	       if(ss == 0.0 && dd == 0.0)
                  outval = apval_ptr[k].rake;

	       while(outval >= 360.0)
	          outval = outval - 360.0;
	       while(outval < 0.0)
	          outval = outval + 360.0;
	       }

            if(dump_slip)
               fprintf(fpw,"%12.5e %12.5e %12.5e %12.5e\n",xx,yy,outval,slipval);
	    else
               fprintf(fpw,"%12.5e %12.5e %12.5e\n",xx,yy,outval);
            }
         }

      *xoff = *xoff + len;
      }
   fclose(fpw);
   }
}

void write_maxsvfX(char *file,struct standrupformat *srf,char *type,int ig,float *maxslip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy;
float *stf;
int it, kp, nx, ny;
int i, j, k, npskip, pflag;

float rperd = 0.017453293;

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

*maxslip = 0.0;
if(pflag)
   {
   apval_ptr = srf[0].srf_apnts.apntvals + npskip;

fprintf(stderr,"nx= %5d ny= %d\n",nx,ny);

   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         if(strncmp(type,"slip",4) == 0)
            {
	    slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

	    if(slip > *maxslip)
	       {
	       *maxslip = slip;
	       kp = k;
	       }
	    }
         }
      }

fprintf(stderr,"kp= %5d\n",kp);

   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"svf svf\n");
   fprintf(fpw,"%d %13.5e\n",apval_ptr[kp].nt1,apval_ptr[kp].dt);

   stf = apval_ptr[kp].stf1;
   for(it=0;it<apval_ptr[kp].nt1;it++)
      fprintf(fpw,"%13.5e\n",stf[it]);

   fclose(fpw);
   }
}

void write_maxsvf(char *file,struct standrupformat *srf,char *type,int ig,float *maxslip)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, dels;
float *stfp, *stfout;
int it, kp, nx, ny, ntout;
int i, j, k, npskip, pflag;

float rperd = 0.017453293;

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

dels = 1.0e+15;
if(pflag)
   {
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         apval_ptr = &(srf[0].srf_apnts.apntvals[k+npskip]);

         if(strncmp(type,"slip",4) == 0)
            {
	    slip = sqrt(apval_ptr->slip1*apval_ptr->slip1 +
	             apval_ptr->slip2*apval_ptr->slip2 +
	             apval_ptr->slip3*apval_ptr->slip3);

	    if(*maxslip < 0.0)
	       {
	       if((slip + *maxslip)*(slip + *maxslip) < dels)
	          {
	          kp = k;
		  dels = (slip + *maxslip)*(slip + *maxslip);

fprintf(stderr,"kp= %d dels= %13.5e\n",kp,dels);
		  }
	       }
	    else if(slip > *maxslip)
	       {
	       *maxslip = slip;
	       kp = k;
	       }
	    }
         }
      }

   apval_ptr = &(srf[0].srf_apnts.apntvals[kp+npskip]);

   ntout = 0;
   if(apval_ptr->nt1 > ntout)
      ntout = apval_ptr->nt1;
   if(apval_ptr->nt2 > ntout)
      ntout = apval_ptr->nt2;
   if(apval_ptr->nt3 > ntout)
      ntout = apval_ptr->nt3;

   stfout = (float *)check_malloc(ntout*sizeof(float));

   for(it=0;it<ntout;it++)
      stfout[it] = 0.0;

   if(apval_ptr->nt1)
      {
      stfp = apval_ptr->stf1;
      for(it=0;it<apval_ptr->nt1;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(apval_ptr->nt2)
      {
      stfp = apval_ptr->stf2;
      for(it=0;it<apval_ptr->nt2;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(apval_ptr->nt3)
      {
      stfp = apval_ptr->stf3;
      for(it=0;it<apval_ptr->nt3;it++)
         stfout[it] = stfout[it] + stfp[it]*stfp[it];
      }

   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"svf svf\n");
   fprintf(fpw,"%d %13.5e\n",ntout,apval_ptr->dt);

   for(it=0;it<ntout;it++)
      fprintf(fpw,"%13.5e\n",sqrt(stfout[it]));

   fclose(fpw);
   free(stfout);
   if(*maxslip < 0.0)
      *maxslip = -(*maxslip);
   }
}

void get_vmax2slip(char *file,struct standrupformat *srf,char *type,int ig)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float elon, elat;
float sn, se, arg;
float len, wid, dx, dy, stk, dip, dtop, shypo, dhypo;
float len2, slip, xx, yy, dels;
float *stfp, *stfout;
int it, kp, nx, ny, ntout;
int i, j, k, npskip, pflag;
float vos, maxv;
int nc;

float rperd = 0.017453293;

if(strcmp(file,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(file,"w");

if(strcmp(srf[0].type,"PLANE") == 0 && ig < srf[0].srf_prect.nseg)
   {
   prseg_ptr = srf[0].srf_prect.prectseg;

   elon = prseg_ptr[ig].elon;
   elat = prseg_ptr[ig].elat;
   nx = prseg_ptr[ig].nstk;
   ny = prseg_ptr[ig].ndip;
   len = prseg_ptr[ig].flen;
   wid = prseg_ptr[ig].fwid;
   stk = prseg_ptr[ig].stk;
   dip = prseg_ptr[ig].dip;
   dtop = prseg_ptr[ig].dtop;
   shypo = prseg_ptr[ig].shyp;
   dhypo = prseg_ptr[ig].dhyp;

   len2 = 0.5*len;
   dx = len/(float)(nx);
   dy = wid/(float)(ny);

   npskip = 0;
   for(i=0;i<ig;i++)
      npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;

   pflag = 1;
   }
else
   {
   elon = -99;
   elat = -99;
   pflag = 0;

   pflag = 1;
   npskip = 0;

   nx = srf[0].srf_apnts.np;
   ny = 1;
   }

dels = 1.0e+15;
if(pflag)
   {
   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
	 k = i + j*nx;

         apval_ptr = &(srf[0].srf_apnts.apntvals[k+npskip]);

         nc = 0;
         vos = 0.0;
         if(apval_ptr->nt1)
            {
            stfp = apval_ptr->stf1;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt1;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip1*apval_ptr->slip1));
            nc++;
            }

         if(apval_ptr->nt2)
            {
            stfp = apval_ptr->stf2;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt2;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip2*apval_ptr->slip2));
            nc++;
            }

         if(apval_ptr->nt3)
            {
            stfp = apval_ptr->stf3;
            maxv = 0.0;
            for(it=0;it<apval_ptr->nt3;it++)
               {
               if(stfp[it]*stfp[it] > maxv)
                  maxv = stfp[it]*stfp[it];
	       }
            vos = vos + sqrt(maxv/(apval_ptr->slip3*apval_ptr->slip3));
            nc++;
            }

         if(nc)
            fprintf(fpw,"%13.5e\n",vos/(float)(nc));

         }
      }

   fclose(fpw);
   }
}

void read_Fvelmodel(char *vfile,struct velmodel *vm)
{
FILE *fpr, *fopfile();
int i, nr;
char str[512];

int nblock = 50;

fpr = fopfile(vfile,"r");

fgets(str,512,fpr);

if(str[0] == '#')
   {
   while(str[0] == '#')
      fgets(str,512,fpr);

   vm->nlay = nblock;
   vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
   vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

   i = 0;
   sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

   vm->dep[i] = vm->th[i];
   vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */

   while(fgets(str,512,fpr) != NULL)
      {
      i++;

      if(i == vm->nlay)
         {
	 vm->nlay = vm->nlay + nblock;
	 vm->vp = (float *)check_realloc(vm->invb2,vm->nlay*sizeof(float));
	 vm->vs = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
	 vm->den = (float *)check_realloc(vm->invb2,vm->nlay*sizeof(float));
	 vm->th = (float *)check_realloc(vm->invb2,vm->nlay*sizeof(float));
	 vm->dep = (float *)check_realloc(vm->invb2,vm->nlay*sizeof(float));
	 vm->mu = (float *)check_realloc(vm->invb2,vm->nlay*sizeof(float));
	 vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
	 }

      sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

      vm->dep[i] = vm->dep[i-1] + vm->th[i];
      vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
      }

   vm->nlay = i+1;
   vm->vp = (float *)check_realloc(vm->vp,vm->nlay*sizeof(float));
   vm->vs = (double *)check_realloc(vm->vs,vm->nlay*sizeof(double));
   vm->den = (float *)check_realloc(vm->den,vm->nlay*sizeof(float));
   vm->th = (float *)check_realloc(vm->th,vm->nlay*sizeof(float));
   vm->dep = (float *)check_realloc(vm->dep,vm->nlay*sizeof(float));
   vm->mu = (float *)check_realloc(vm->mu,vm->nlay*sizeof(float));
   vm->invb2 = (double *)check_realloc(vm->invb2,vm->nlay*sizeof(double));
   }
else
   {
   sscanf(str,"%d",&vm->nlay);

   vm->vp = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->vs = (double *)check_malloc(vm->nlay*sizeof(double));
   vm->den = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->th = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->dep = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->mu = (float *)check_malloc(vm->nlay*sizeof(float));
   vm->invb2 = (double *)check_malloc(vm->nlay*sizeof(double));

   for(i=0;i<vm->nlay;i++)
      {
      fgets(str,512,fpr);
      sscanf(str,"%f %f %lf %f",&vm->th[i],&vm->vp[i],&vm->vs[i],&vm->den[i]);

      if(i==0)
         vm->dep[i] = vm->th[i];
      else
         vm->dep[i] = vm->dep[i-1] + vm->th[i];

      vm->mu[i] = vm->vs[i]*vm->vs[i]*vm->den[i]*1.0e+10;  /* in CMS units */
      }
   }

fclose(fpr);
}

void get_moment(struct standrupformat *srf,struct velmodel *vm)
{
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float sigma, aslip, tarea, tmom, slip;
float amax, amin;
int np;
int j, k;

amax = -1.0e+15;
amin = 1.0e+15;
aslip = 0.0;
tarea = 0.0;
tmom = 0.0;
np = srf[0].srf_apnts.np;
apval_ptr = srf[0].srf_apnts.apntvals;

fprintf(stderr,"np= %d\n",np);

for(k=0;k<np;k++)
   {
   slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

   j = 0;
   while(vm->dep[j] < apval_ptr[k].dep)
      j++;

   aslip = aslip + slip;
   tarea = tarea + apval_ptr[k].area*1.0e-10;
   tmom = tmom + slip*apval_ptr[k].area*vm->mu[j];

   /*
fprintf(stderr,"dep= %13.5e mu= %13.5e\n",apval_ptr[k].dep,vm->mu[j]);

   if(apval_ptr[k].nt1 == 0 && apval_ptr[k].nt2 == 0 && apval_ptr[k].nt3 == 0)
      fprintf(stderr,"ZERO %5d %13.5e %13.5e %13.5e\n",k,apval_ptr[k].slip1,apval_ptr[k].slip2,apval_ptr[k].slip3);
*/

   if(apval_ptr[k].area < amin)
      amin = apval_ptr[k].area;
   if(apval_ptr[k].area > amax)
      amax = apval_ptr[k].area;
   }

aslip = aslip/(float)(np);
sigma = 0.0;
for(k=0;k<np;k++)
   {
   slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	             apval_ptr[k].slip2*apval_ptr[k].slip2 +
	             apval_ptr[k].slip3*apval_ptr[k].slip3);

   sigma = sigma + (aslip - slip)*(aslip - slip);
   }

fprintf(stderr,"Total area= %13.5e km*km\n",tarea);
fprintf(stderr,"Average area= %13.5e km*km\n",tarea/(float)(np));
fprintf(stderr,"Min. area= %13.5e km*km\n",amin*1.0e-10);
fprintf(stderr,"Max. area= %13.5e km*km\n",amax*1.0e-10);
fprintf(stderr,"Average slip= %13.5e cm (sigma/aslip= %.3f)\n",aslip,sqrt(sigma/(np-1))/aslip);
fprintf(stderr,"Total moment= %13.5e dyne-cm (Mw= %.2f)\n",tmom,0.66667*(log(tmom)/log(10))-10.7);
}

int write_lld(char *file,struct standrupformat *srf,int ig,float *dmin,float *dmax)
{
FILE *fpw, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
float slip;
int nx, ny;
int i, j, k, npskip;
int igst, ignd, igp;
char str[512], frmt[32];

igst = 0;
ignd = 1;
if(strcmp(srf[0].type,"PLANE") == 0)
   {
   if(ig < 0)
      {
      igst = 0;
      ignd = srf[0].srf_prect.nseg;
      }
   else
      {
      igst = ig;
      ignd = igst + 1;
      }
   }

if(strcmp(file,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(file,"w");

for(igp=igst;igp<ignd;igp++)
   {
   npskip = 0;
   nx = srf[0].srf_apnts.np;
   ny = 1;

   if(strcmp(srf[0].type,"PLANE") == 0)
      {
      prseg_ptr = srf[0].srf_prect.prectseg;

      nx = prseg_ptr[igp].nstk;
      ny = prseg_ptr[igp].ndip;

      npskip = 0;
      for(i=0;i<igp;i++)
         npskip = npskip + prseg_ptr[i].nstk*prseg_ptr[i].ndip;
      }

   apval_ptr = srf[0].srf_apnts.apntvals + npskip;

   for(j=0;j<ny;j++)
      {
      for(i=0;i<nx;i++)
         {
         k = i + j*nx;

         if(apval_ptr[k].dep >= *dmin && apval_ptr[k].dep <= *dmax)
	    {
	    slip = sqrt(apval_ptr[k].slip1*apval_ptr[k].slip1 +
	                            apval_ptr[k].slip2*apval_ptr[k].slip2 +
				                            apval_ptr[k].slip3*apval_ptr[k].slip3);
            fprintf(fpw,"%12.5f %12.5f %12.5e %12.5e\n",apval_ptr[k].lon,apval_ptr[k].lat,apval_ptr[k].dep,slip);
	    }
         }
      }
   }

fclose(fpw);
}
