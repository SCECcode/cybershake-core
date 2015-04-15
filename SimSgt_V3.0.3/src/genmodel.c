/*
   genmodel.c contains the following functions:

      adj_medsten_fs()
      avgmedslice()
      copybasinprof()
      do_digit()
      do_elp()
      do_geom()
      do_init()
      do_open_elp()
      do_prof()
      do_topo()
      getlineRWG()
      genmod3d()
      getmedslice()
      init_name()
      modsten2mplot()
      set_max_depths()
      set_vminvmax_genmod()
*/

#include "include.h"

genmod3d(nprof,mfile,xorg,yorg,nx,ny,nz,h,medsten,medprofs,fs,medbnd)
char *mfile;
struct medstencil *medsten;
struct medprof *medprofs;
float *xorg, *yorg, *h;
int nx, ny, nz, fs, medbnd, *nprof;
{
FILE *fpr, *fopfile();
char *words[128], *line, *getlineRWG();
int nword;
int ip = 0;
char buf[256];

int linenum = 0;
int ungetflag = 0;

do_init(nx*ny,medsten,nz);

fpr = fopfile(mfile,"r");

while((line = getlineRWG(fpr,&linenum,&ungetflag,buf)) != NULL)
   {
   if(line[0] == '\0') continue;     /* blank line */
   if(line[0] == '#') continue;      /* comment */
   nword = getwords(line,words);
   if(nword == 0) continue;          /* nothing on line */

   if(strcmp("DEF",words[0]) == 0)
      {
      do_prof(fpr,&medprofs[ip],&words[1],nz,h,fs,&linenum,&ungetflag,buf,medbnd);
      ip++;
      continue;
      }
   if(strcmp("GEOM",words[0]) == 0)
      {
      do_geom(nx,ny,nz,h,medsten,nword,words);
      continue;
      }
   if(strcmp("DIG",words[0]) == 0)
      {
      do_digit(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   if(strcmp("ELP",words[0]) == 0)
      {
      do_elp(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   if(strcmp("OPELP",words[0]) == 0)
      {
      do_open_elp(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   if(strncmp("TOP",words[0],3) == 0)
      {
      do_topo(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   if(strncmp("CUBE",words[0],4) == 0)
      {
      do_cube(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   if(strncmp("DIPPING_LAYER",words[0],13) == 0)
      {
      do_diplay(nx,ny,nz,h,medsten,xorg,yorg,nword,words);
      continue;
      }
   }

fclose(fpr);

if(fs) /* Adjust media stencil to add extra row at model top for free-surface */
   adj_medsten_fs(nx,ny,nz,medsten,fs);

/*
   Set max. depths for basin profiles and check to ensure that input
   profiles are deep enough.
*/

set_max_depths(medsten,medprofs,ip,nx,ny);

*nprof = ip;
}

adj_medsten_fs(nx,ny,nz,medsten,fs)
struct medstencil *medsten;
int nx, ny, nz, fs;
{
int ix, iy;

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      if(medsten[ix + iy*nx].idep0 != 0)
         medsten[ix + iy*nx].idep0 = medsten[ix + iy*nx].idep0 + fs;

      medsten[ix + iy*nx].idep1 = medsten[ix + iy*nx].idep1 + fs;
      if(medsten[ix + iy*nx].idep1 > nz)
         medsten[ix + iy*nx].idep1 = nz;
      }
   }
}

do_digit(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
FILE *fpr, *fopfile();
int ix, iy, iz, nw;
float *dval;
char str[512];
char fname[128];
char name[4];
char units[32];
float unitconv = 1.0;

dval = (float *) check_malloc (nx*ny*sizeof(float));

init_name(name,4);

if(nword != 4)
   {
   fprintf(stdout,"nword = %d, should be 4\n",nword);
   exit(-1);
   }
   
strcpy(fname,words[1]);
strcpy(units,words[2]);
strncpy(name,words[3],3);

if(units[0] == 'm' || units[0] == 'M')       /* meters */
   unitconv = 0.001;
else if(units[0] == 'k' || units[0] == 'K')  /* kilometers */
   unitconv = 1.0;

fprintf(stderr,"*** do_digit(): unitconv= %f\n",unitconv);

fpr = fopfile(fname,"r");

/* First, skip comment lines */
 
fgets(str,512,fpr);
while(str[0] == '#')
   {
   fprintf(stderr,"                %s",str);
   fgets(str,512,fpr);
   }
 
/*
   If there are no comment lines, get depth values stored in the character
   string "str", and then get remaining depth values for iy=0 from
   first part of file,
*/
 
nw = 0;                /* remove NEWLINE character from "str" */
while(str[nw] != '\n')
   nw++;
str[nw] = '\0';
 
nw = getwords(str,words);
for(ix=0;ix<nw;ix++)
   dval[ix] = atof(words[ix]);
 
for(ix=nw;ix<nx;ix++)
   fscanf(fpr,"%f",&dval[ix]);
 
/* Next, get remaining depth values (iy>0) from input file */
 
for(iy=1;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      fscanf(fpr,"%f",&dval[ix + iy*nx]);
   }
fclose(fpr);

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      if(dval[ix + iy*nx] > 0.0)   /*  within basin  */
         {
         iz = (int)(unitconv*dval[ix + iy*nx]/(*dh) + 0.5);
	 if(iz > nz)
	    iz = nz;
         strncpy(medsten[iy*nx + ix].med_name,name,4);
         medsten[iy*nx + ix].idep1 = iz;
         }
      }
   }
free(dval);
}

do_geom(nx,ny,nz,dh,medsten,nword,words)
float *dh;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
float x0, y0, z0, x, y, z, dep, rad, r, x2, y2;
int ix, iy, iz, iflag, idep;
char type[16];
char name[4];

init_name(name,4);
strcpy(type,words[1]);

if(strncmp("cylind",type,6) == 0)
   {
   if(nword != 7)
      {
      printf("nword = %d, should be 7\n",nword);
      exit(-1);
      }

   x0 = atof(words[2]);
   y0 = atof(words[3]);
   dep = atof(words[4]);
   rad = atof(words[5]);
   strncpy(name,words[6],3);

   idep = (int)(dep/(*dh) + 0.5);

   for(iy=0;iy<ny;iy++)
      {
      y = iy*(*dh)-y0;
      y2 = y*y;
      for(ix=0;ix<nx;ix++)
         {
         x = ix*(*dh)-x0;
         r = sqrt(x*x + y2);
   
         if(r <= rad)   /*  within cylinder  */
	    {
            strncpy(medsten[iy*nx + ix].med_name,name,4);
	    medsten[iy*nx + ix].idep1 = idep;
	    }
         }
      }
   }
else if(strncmp("sphere",type,6) == 0)
   {
   if(nword != 7)
      {
      printf("nword = %d, should be 7\n",nword);
      exit(-1);
      }

   x0 = atof(words[2]);
   y0 = atof(words[3]);
   z0 = atof(words[4]);
   rad = atof(words[5]);
   strncpy(name,words[6],3);

   for(iy=0;iy<ny;iy++)
      {
      y = iy*(*dh) - y0;
      y2 = y*y;
      for(ix=0;ix<nx;ix++)
         {
         x = ix*(*dh) - x0;
	 x2 = x*x;

	 iflag = 0;
         for(iz=0;iz<nz;iz++)
            {
            z = iz*(*dh) - z0;
            r = sqrt(x2 + y2 + z*z);
   
            if(r <= rad)   /*  within sphere  */
	       {
	       if(iflag == 0) /* just entering sphere */
		  {
		  iflag = 1;
                  strncpy(medsten[iy*nx + ix].med_name,name,4);
	          medsten[iy*nx + ix].idep0 = iz;
		  }
	       }
            else          /*  outside sphere  */
	       {
	       if(iflag == 1) /* just exiting sphere */
		  {
		  iflag = 0;
	          medsten[iy*nx + ix].idep1 = iz-1;
		  }
	       }
	    }
         }
      }
   }
else
   {
   printf("*** geometric type %s is not a valid choice, exiting...\n",type);
   exit(-1);
   }
}

do_elp(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
double arg;
float cosT, sinT;
float a, b, c, h, k, l, theta;
float inva2, invb2, c2, x, y, z, xprm, yprm;
float sum, check;
int ix, iy, iz;
char name[4];

init_name(name,4);

if(nword != 9)
   {
   fprintf(stdout,"nword = %d, should be 9\n",nword);
   exit(-1);
   }
   
a = atof(words[1]);
b = atof(words[2]);
c = atof(words[3]);
h = atof(words[4]);
k = atof(words[5]);
l = atof(words[6]);
theta = atof(words[7]);
strncpy(name,words[8],3);

if(l > 0.0)
   {
   fprintf(stderr,"\n           **********************\n");
   fprintf(stderr,"Program assumes l to be negative (ie. above the free\n");
   fprintf(stderr,"surface).  Changing l to -l.  Check model to be sure\n");
   fprintf(stderr,"it is what you want.  Basin name=%s.\n\n",name);
   l = -l;
   }

arg = theta*RPERD;
cosT = cos(arg);
sinT = sin(arg);

inva2 = 1.0/(a*a);
invb2 = 1.0/(b*b);
c2 = c*c;

for(iy=0;iy<ny;iy++)
   {
   y = iy*(*dh) + *yorg;
   for(ix=0;ix<nx;ix++)
      {
      x = ix*(*dh) + *xorg;

      xprm = -(x-h)*sinT - (y-k)*cosT;
      yprm = (x-h)*cosT - (y-k)*sinT;

      sum = c2*(1.0 - inva2*xprm*xprm - invb2*yprm*yprm);
      if(sum >= 0.0)
         {
         check = sqrt(sum) + l;
         for(iz=0;iz<nz;iz++)
            {
            z = iz*(*dh);
            if(z <= check)   /*  within ellipsoid  */
	       continue;
	    else
	       {
               strncpy(medsten[iy*nx + ix].med_name,name,4);
	       medsten[iy*nx + ix].idep1 = iz;
	       break;
	       }
            }
         }
      }
   }
}

do_open_elp(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
double arg;
float cosT, sinT;
float a, b, c, h, k, l, theta;
float inva2, invb2, c2, x, y, z, xprm, yprm;
float sum, check;
int ix, iy, iz;
char name[4];

init_name(name,4);

if(nword != 9)
   {
   fprintf(stdout,"nword = %d, should be 9\n",nword);
   exit(-1);
   }
   
a = atof(words[1]);
b = atof(words[2]);
c = atof(words[3]);
h = atof(words[4]);
k = atof(words[5]);
l = atof(words[6]);
theta = atof(words[7]);
strncpy(name,words[8],3);

if(l > 0.0)
   {
   fprintf(stderr,"\n           **********************\n");
   fprintf(stderr,"Program assumes l to be negative (ie. above the free\n");
   fprintf(stderr,"surface).  Changing l to -l.  Check model to be sure\n");
   fprintf(stderr,"it is what you want.  Basin name=%s.\n\n",name);
   l = -l;
   }

arg = theta*RPERD;
cosT = cos(arg);
sinT = sin(arg);

inva2 = 1.0/(a*a);
invb2 = 1.0/(b*b);
c2 = c*c;

for(iy=0;iy<ny;iy++)
   {
   y = iy*(*dh) + *yorg;
   for(ix=0;ix<nx;ix++)
      {
      x = ix*(*dh) + *xorg;

      xprm = -(x-h)*sinT - (y-k)*cosT;
      if(xprm > 0.0)
	 xprm = 0.0;
      yprm = (x-h)*cosT - (y-k)*sinT;

      sum = c2*(1.0 - inva2*xprm*xprm - invb2*yprm*yprm);
      if(sum >= 0.0)
         {
         check = sqrt(sum) + l;
         for(iz=0;iz<nz;iz++)
            {
            z = iz*(*dh);
            if(z <= check)   /*  within ellipsoid  */
	       continue;
	    else
	       {
               strncpy(medsten[iy*nx + ix].med_name,name,4);
	       medsten[iy*nx + ix].idep1 = iz;
	       break;
	       }
            }
         }
      }
   }
}

do_topo(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
double arg;
float x0, y0, z0, a0;
float a2, x, y, z, sum;
float one = 1.0;
int ix, iy, iz;
char name[4];

init_name(name,4);

if(nword != 6)
   {
   fprintf(stdout,"nword = %d, should be 6\n",nword);
   exit(-1);
   }
   
x0 = atof(words[1]);            /* x location of center of mound */
y0 = atof(words[2]);            /* y location of center of mound */
z0 = atof(words[3]);            /* height of mound */
a0 = atof(words[4]);            /* half-width radius of mound */
strncpy(name,words[5],3);

a2 = -one/(a0*a0);

for(iy=0;iy<ny;iy++)
   {
   y = iy*(*dh) + *yorg - y0;
   for(ix=0;ix<nx;ix++)
      {
      x = ix*(*dh) + *xorg - x0;

      arg = a2*(x*x + y*y);
      sum = z0*(one - exp(arg)) + 4.0*(*dh); /* at least 4*h below model top */

      for(iz=0;iz<nz+1;iz++)
         {
         z = iz*(*dh);
         if(z < sum && iz < nz)   /*  above topography  */
            continue;
	 else
	    {
            strncpy(medsten[iy*nx + ix].med_name,name,4);
	    medsten[iy*nx + ix].idep1 = iz;
	    break;
	    }
         }
      }
   }
}

do_cube(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
double arg;
int x0, y0, z0, x1, y1, z1;
int ix, iy;
char name[4];

init_name(name,4);

if(nword != 8)
   {
   fprintf(stdout,"nword = %d, should be 8\n",nword);
   exit(-1);
   }
   
x0 = atoi(words[1]);
y0 = atoi(words[2]);
z0 = atoi(words[3]);
x1 = atoi(words[4]);
y1 = atoi(words[5]);
z1 = atoi(words[6]);
strncpy(name,words[7],3);

if(x0 < 0)
   x0 = 0;
if(y0 < 0)
   y0 = 0;
if(z0 < 0)
   z0 = 0;
if(x1 > nx)
   x1 = nx;
if(y1 > ny)
   y1 = ny;
if(z1 > nz)
   z1 = nz;

for(iy=y0;iy<y1;iy++)
   {
   for(ix=x0;ix<x1;ix++)
      {
      strncpy(medsten[iy*nx + ix].med_name,name,4);
      medsten[iy*nx + ix].idep0 = z0;
      medsten[iy*nx + ix].idep1 = z1;
      }
   }
}

do_diplay(nx,ny,nz,dh,medsten,xorg,yorg,nword,words)
float *dh, *xorg, *yorg;
struct medstencil *medsten;
int nword, nx, ny, nz;
char *words[];
{
float theta, beta, b;
double c0, c1;
int ix, iy, z1;
char name[4];

double pi = 3.141592654;

init_name(name,4);

if(nword != 5)
   {
   fprintf(stdout,"nword = %d, should be 5\n",nword);
   exit(-1);
   }
   
theta = atof(words[1]);
beta = atof(words[2]);
b = atoi(words[3]);
strncpy(name,words[4],3);

if(theta >= 90.0)
   theta = 89.99;
if(theta <= -90.0)
   theta = -89.99;

c0 = cos(beta*pi/180.0)*tan(theta*pi/180.0);
c1 = sin(beta*pi/180.0)*tan(theta*pi/180.0);

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      z1 = 1 + (int)(ix*c0 + iy*c1 + b + 0.5);
fprintf(stderr,"%5d  %5d  %5d\n",ix,iy,z1);
      if(z1 > nz)
	 z1 = nz;

      if(z1 >= 0)
	 {
         strncpy(medsten[iy*nx + ix].med_name,name,4);
         medsten[iy*nx + ix].idep0 = 0;
         medsten[iy*nx + ix].idep1 = z1;
	 }
      }
   }
}

do_init(n,msten,ndep)
struct medstencil *msten;
int n, ndep;
{
int i;
char symb[4];

symb[0] = 'H';
symb[1] = 'S';
symb[2] = 'T';
symb[3] = '\0';

for(i=0;i<n;i++)
   {
   strncpy(msten[i].med_name,symb,4);
   msten[i].idep0 = 0;
   msten[i].idep1 = ndep;
   }
}

do_prof(fp,medprofs,words,n,h,fs,linenum,ungetflag,buf,medbnd)
struct medprof *medprofs;
float *h;
int n, fs, *linenum, *ungetflag, medbnd;
FILE *fp;
char *words[], *buf;
{
double atof();
float *lam2mu, *lam, *invmu, *rho, *qp, *qs;
float alpha, beta, den, qp0, qs0;
char *line, *getlineRWG();
int ist, iend, i;
float dep;

float one = 1.0;
float two = 2.0;
int ndeps = 0;

lam2mu = (medprofs->medptr);
lam    = (medprofs->medptr) +   n;
invmu  = (medprofs->medptr) + 2*n; /* store as 1/mu */
rho    = (medprofs->medptr) + 3*n;
qp     = (medprofs->medptr) + 4*n;
qs     = (medprofs->medptr) + 5*n;

init_name(medprofs->prof_name,4);
strncpy(medprofs->prof_name,words[0],3);

ist = 0;
while((line = getlineRWG(fp,linenum,ungetflag,buf)) != NULL)
   {
   if(line[0] == '\0') continue;     /* blank line */
   if(line[0] == '#') continue;      /* comment */
   if(line[0] != ' ' && line[0] != '\t')      /* end of profile */
      {
      *ungetflag = 1;
      break;
      }
   i = getwords(line,words);
   if(i == 0) continue;          /* nothing on line */
   if(i != 6)
      {
      fprintf(stderr,"\n***********  ERROR  ***********\n");
      fprintf(stderr,"Need P velocity, S velocity, density, Q and depth value for each\n");
      fprintf(stderr,"entry in profile.  Line no. = %d.\n", *linenum);
      exit(-1);
      }

/*
      words[0] = P velocity
      words[1] = S velocity
      words[2] = density
      words[3] = Q for p-waves
      words[4] = Q for s-waves
      words[5] = depth to bottom of layer
*/

   /* round to nearest depth increment */

   if(medbnd == 0)
      dep = atof(words[5])/(*h) + 1.0;
   else
      dep = atof(words[5])/(*h) + 0.5;

   iend = (int)(dep);
   if(iend > n)
      iend = n;

   alpha = atof(words[0]);
   beta = atof(words[1]);
   den = atof(words[2]);
   qp0 = atof(words[3]);
   qs0 = atof(words[4]);

/*
fprintf(stderr,"%7.4f %7.4f %7.4f %7.4f %7.4f %3d\n",alpha,beta,den,qp0,qs0,iend-1);
*/

   if(qp0 <= 0.0)
      qp0 = 1.0e+10;
   if(qs0 <= 0.0)
      qs0 = 1.0e+10;

   for(i=ist;i<iend;i++)
      {
      lam2mu[i] = alpha*alpha*den;
      invmu[i] = beta*beta*den;

      lam[i] = lam2mu[i] - two*invmu[i];
      invmu[i] = one/invmu[i];

      rho[i] = den;
      qp[i] = qp0;
      qs[i] = qs0;

      ndeps++;
      }
   ist = iend;
   }

/*
     Add dummy row to top of model for free-surface.  The rest of the model
     is shifted down by one grid point.
*/
if(fs)
   {
   ndeps = ndeps + fs;
   if(ndeps > n)
      ndeps = n;

   for(i=ndeps-1;i>=fs;i--)
      {
      lam2mu[i] = lam2mu[i-fs];
      lam[i] = lam[i-fs];
      invmu[i] = invmu[i-fs];
      rho[i] = rho[i-fs];
      qp[i] = qp[i-fs];
      qs[i] = qs[i-fs];
      }
   }

medprofs->ndep = ndeps;
}

set_max_depths(ms,mdp,nprof,nx,ny)
struct medstencil *ms;
struct medprof *mdp;
int nprof, nx, ny;
{
int ispot, iprof;
int ix, iy;

for(iprof=0;iprof<nprof;iprof++)
   mdp[iprof].mxdp = -1;

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ispot = iy*nx + ix;
   
      iprof = 0;
      while(strcmp(ms[ispot].med_name,mdp[iprof].prof_name) != 0)
         iprof++;

      if(ms[ispot].idep1 > mdp[iprof].mxdp)
         mdp[iprof].mxdp = ms[ispot].idep1;
      }
   }

for(iprof=0;iprof<nprof;iprof++)
   {
   if(mdp[iprof].mxdp > mdp[iprof].ndep)
      {
      fprintf(stderr,"\n***********  ERROR  ***********\n");
      fprintf(stderr,"Profile for media name '%s' ",mdp[iprof].prof_name);
      fprintf(stderr,"does not go deep enough.\n");
      fprintf(stderr,"**** maxdep= %d, ndep= %d\n",mdp[iprof].mxdp,mdp[iprof].ndep);
      exit(-1);
      }
   }
}

set_vminvmax_gen(ms,mdp,nprof,fdc,nx,ny,nz)
struct medstencil *ms;
struct medprof *mdp;
struct fdcoefs *fdc;
int nprof, nx, ny, nz;
{
int *prof_list;
int ispot, iprof;
int ix, iy;
float *lam2mu, *invmu, *rho;
float a2, b2, den, vmin, vmin2, vmax;
int iz, izmin;
float one = 1.0;

/*
   Set lowest valid shear velocity to be 50 m/sec (0.05 km/sec).
   This avoids problems in regions where mu -> 0.  The variable 'b2floor'
   is set to 0.05*0.05 = 0.0025 as a cutoff to find 'vmin'.
*/

float b2floor = 0.0025;

/*
   Find which profiles are used for current model space.
*/

prof_list = (int *) check_malloc (sizeof(int)*nprof);

for(iprof=0;iprof<nprof;iprof++)
   prof_list[iprof] = -1;           /* initialize all profiles unactive */

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ispot = iy*nx + ix;
   
      iprof = 0;
      while(strcmp(ms[ispot].med_name,mdp[iprof].prof_name) != 0)
         iprof++;

      prof_list[iprof] = 1;         /* profile is used */
      }
   }

/*
   Find vmin and vmax only for active profiles.
*/

vmin = 1.0e+15;
vmax = -1.0;
for(iprof=0;iprof<nprof;iprof++)
   {
   if(prof_list[iprof] == 1)
      {
      lam2mu = (mdp[iprof].medptr);
      invmu  = (mdp[iprof].medptr) + 2*nz;
      rho    = (mdp[iprof].medptr) + 3*nz;

      for(iz=0;iz<mdp[iprof].mxdp;iz++)
	 {
	 den = one/rho[iz];
	 a2 = lam2mu[iz]*den;
	 b2 = den/invmu[iz];

	 if(b2 > b2floor && b2 < vmin)  /* store as velocity squared for now */
	    vmin = b2;
	 if(a2 > vmax)    /* store as velocity squared for now */
	    vmax = a2;
	 }
      }
   }

/*
   Find 'izord2'.  This is the depth index below which all velocities are
   at least two times larger than 'vmin' (remember, 'vmin' is stored as
   velocity squared at this point).  For this part of the model,
   second-order FD operators (ord2) can be used without substantial loss
   of accuracy, but with added computational efficiency.

   To be tested: 05 May 1994, RWG
*/

vmin2 = 4.0*vmin;   /* 2*vmin squared */
izmin = 0;
for(iprof=0;iprof<nprof;iprof++)
   {
   if(prof_list[iprof] == 1)
      {
      lam2mu = (mdp[iprof].medptr);
      invmu  = (mdp[iprof].medptr) + 2*nz;
      rho    = (mdp[iprof].medptr) + 3*nz;

      iz = (mdp[iprof].mxdp) - 1;
      b2 = one/(rho[iz]*invmu[iz]);
      while(b2 > vmin2 && iz > 0)    /* store as velocity squared for now */
         {
         iz--;
         b2 = one/(rho[iz]*invmu[iz]);
	 }

      if((iz+1) > izmin)
         izmin = iz + 1;
      }
   }

fdc->vmin = sqrt(vmin);
fdc->vmax = sqrt(vmax);
fdc->izord2 = izmin + 1;  /* add 1 just to be sure */

free(prof_list);
}

/*
   reedmedslice() reads the media parameters into an xz slice of the model.

   These values are loaded into the array medf as follows:

      medf           -> lambda + 2*mu
      medf +   nx*nz -> lambda
      medf + 4*nx*nz -> 1/mu
      medf + 7*nx*nz -> rho
      medf + 9*nx*nz -> qp  (later turned into pk coefs.)
      medf +10*nx*nz -> qs  (later turned into sk coefs.)
      medf +11 nx*nz -> lambda (raw, never modified by avg. or atten. )
      medf +12*nx*nz -> mu     (raw, never modified by avg. or atten. )

   After reedmedslice() returns to main(), a call to avgmedslice() computes the
   smoothed values of rigidity (mu) and bouyancy (1/rho) which are then
   stored in the array medf and are ready to be used in the differencing
   operations.

*/

void reedmedslice(int fdp,int fds,int fdd,float *medf,int nx,int ny,int nz,int fs,struct qvalues *qv,int iy,int dw,float *qbmax,float *qbmin,int swapb)
{
float *a, *b;
float *lam2mu, *mu, *lam, *rho, *qp, *qs;
float *lamraw, *muraw;
int i, ip, k, ix, iz, nn, shft;
int powr, xpowr, ypowr, zpowr;
float rfac, qbndf, qfac0, qp0, qs0, qfac;
float aa, bb, rr;

float one = 1.0;
float two = 2.0;

rfac = one/(dw - NBND_PAD - 1);
qbndf = (*qbmin)/(*qbmax);
qfac0 = exp(rfac*log(qbndf));

nn = nx*nz;

a = medf + 5*nn;  /* temporary, not used until avgmedslice */
b = medf + 6*nn;  /* temporary, not used until avgmedslice */

lam2mu = medf;
lam    = medf + nn;
mu     = medf + 4*nn;
rho    = medf + 7*nn;
qp     = medf + 9*nn;
qs     = medf + 10*nn;

/*
   Store raw values of lambda and mu into medf array (11*nx*nz and
   12*nx*nz offsets).  These values will never be modified by averaging
   or attenuation (phase velocity adjustments).
*/

lamraw = medf + 11*nn;
muraw  = medf + 12*nn;

reed(fdp,a,nn*sizeof(float)); /* p velocity file */
reed(fds,b,nn*sizeof(float)); /* s velocity file */
reed(fdd,rho,nn*sizeof(float)); /* density file */

if(swapb == 1) /* need to swap bytes */
   {
   swap_in_place(nn,((char *)(a)));
   swap_in_place(nn,((char *)(b)));
   swap_in_place(nn,((char *)(rho)));
   }

for(ix=0;ix<nx;ix++)
   {
   lam2mu[ix] = a[ix]*a[ix]*rho[ix];
   mu[ix] = b[ix]*b[ix]*rho[ix];

   lam[ix] = lam2mu[ix] - two*mu[ix];
   mu[ix] = one/mu[ix];
   lamraw[ix] = lam[ix];
   muraw[ix] = one/mu[ix];

   get_qvals(&qp[ix],&qs[ix],&a[ix],&b[ix],qv);

   if(qp[ix] <= 0.0 || qs[ix] <= 0.0)
      {
      fprintf(stderr,"Problem with zero/negative Q, exiting...\n");
      exit(-99);
      }
   }

shft = 0;
if(fs)      /* shift model down one grid point for free surface */
   shft = 1;

for(iz=nz-1;iz>=1;iz--)
   {
   for(ix=0;ix<nx;ix++)
      {
      i = iz*nx + ix;
      ip = (iz-shft)*nx + ix;

      lam2mu[i] = a[ip]*a[ip]*rho[ip];
      mu[i] = b[ip]*b[ip]*rho[ip];

      lam[i] = lam2mu[i] - two*mu[i];
      mu[i] = one/mu[i];
      lamraw[i] = lam[i];
      muraw[i] = one/mu[i];
      rho[i] = rho[ip];

      get_qvals(&qp[i],&qs[i],&a[ip],&b[ip],qv);

/* damping region for absorbing boundary */

      xpowr = 0;
      if(ix < dw)
	 xpowr = dw - ix;
      else if(ix >= nx-dw)
	 xpowr = ix - (nx - 1 - dw);

      ypowr = 0;
      if(iy < dw)
	 ypowr = dw - iy;
      else if(iy >= ny-dw)
	 ypowr = iy - (ny - 1 - dw);

      zpowr = 0;
      if(iz < dw && fs == 0)
	 zpowr = dw - iz;
      else if(iz >= nz-dw)
	 zpowr = iz - (nz - 1 - dw);

      if(xpowr || ypowr || zpowr)
         {
         powr = xpowr + ypowr + zpowr - 1;
	 qfac = one;
         while(powr--)
            qfac = qfac*qfac0;

         qp0 = qp[i];
         if(qp0 > 2*(*qbmax))
            qp0 = 2*(*qbmax);

         qs0 = qs[i];
         if(qs0 > (*qbmax))
            qs0 = (*qbmax);

         qp[i] = qfac*qp0;
         qs[i] = qfac*qs0;

         if(qp[i] < 5.0)
            qp[i] = 5.0;
         if(qs[i] < 5.0)
            qs[i] = 5.0;

         qp[i] = -qp[i];
         qs[i] = -qs[i];
         }
      }
   }
}

get_qvals(qp,qs,vp,vs,qv)
float *qp, *qs, *vp, *vs;
struct qvalues *qv;
{
float bdel, babs;
int i;

float del = 0.0005;

if(qv->qpfrac > 0.0)
   *qp = (*vp)*qv->qpfrac;
else
   {
   *qp = -1;

   for(i=0;i<qv->n;i++)
      {
      if(qv->vp[i] > ((*vp)-del) && qv->vp[i] < ((*vp)+del))
         (*qp) = qv->qp[i];
      }

   if((*qp) < 0.0) /* if vp is not found in table, use closest value */
      {
      bdel = 1.0e+15;
      for(i=0;i<qv->n;i++)
         {
         babs = qv->vp[i] - (*vp);
         if(babs < 0.0)
            babs = -babs;
         if(babs < bdel)
            {
            bdel = babs;
            *(qp) = qv->qp[i];
            }
         }
      }
   }

if((*qp) < 0.0)
   {
   fprintf(stderr,"Problem with negative Qp= %f Vp= %f, exiting...\n",(*qp),(*vp));
   exit(-99);
   }

if(qv->qsfrac > 0.0)
   *qs = (*vs)*qv->qsfrac;
else
   {
   *qs = -1;

   for(i=0;i<qv->n;i++)
      {
      if(qv->vs[i] > ((*vs)-del) && qv->vs[i] < ((*vs)+del))
         (*qs) = qv->qs[i];
      }

   if((*qs) < 0.0) /* if vs is not found in table, use closest value */
      {
      bdel = 1.0e+15;
      for(i=0;i<qv->n;i++)
         {
         babs = qv->vs[i] - (*vs);
         if(babs < 0.0)
            babs = -babs;
         if(babs < bdel)
            {
            bdel = babs;
            *(qs) = qv->qs[i];
            }
         }
      }
   }

if((*qs) < 0.0)
   {
   fprintf(stderr,"Problem with negative Qs, exiting...\n");
   exit(-99);
   }

if(qv->qpqs_factor > 0.0)
   *qp = (*qs)*qv->qpqs_factor;
}

/*
   getmedslice() loads the media parameters into an xz slice of the model.
   First, the 'HOST' (HST) media profiles are loaded into the array medf,
   then the appropriate basin profile values are loaded over these.  The media
   profiles stored in the array medptr in the structure mdp are as follows:

      mdp[iprof].medptr        -> lambda + 2*mu
      mdp[iprof].medptr +   nz -> lambda
      mdp[iprof].medptr + 2*nz -> 1/mu
      mdp[iprof].medptr + 3*nz -> rho
      mdp[iprof].medptr + 4*nz -> qp factor
      mdp[iprof].medptr + 5*nz -> qs factor

   These values are loaded into the array medf as follows:

      medf           -> lambda + 2*mu
      medf +   nx*nz -> lambda
      medf + 4*nx*nz -> 1/mu
      medf + 7*nx*nz -> rho
      medf + 9*nx*nz -> qp  (later turned into pk coefs.)
      medf +10*nx*nz -> qs  (later turned into sk coefs.)
      medf +11 nx*nz -> lambda (raw, never modified by avg. or atten. )
      medf +12*nx*nz -> mu     (raw, never modified by avg. or atten. )

   After getmedslice() returns to main(), a call to avgmedslice() computes the
   smoothed values of rigidity (mu) and bouyancy (1/rho) which are then
   stored in the array medf and are ready to be used in the differencing
   operations.

*/

getmedslice(ms,mdp,nprof,medf,nx,ny,nz,fs,iy,dw,qbmax,qbmin)
struct medstencil *ms;
struct medprof *mdp;
float *medf, *qbmax, *qbmin;
int nprof, nx, ny, nz, fs, iy, dw;
{
float *mptr, *m2ptr;
float *mptr1, *mptr2, *mptr3, *mptr4;
int iprof, ix, izst, izend;

float *lam2mu, *mu, *lam;
float *qp, *qs, qp0, qs0, qfac;
float rfac, qbndf, qfac0;
int i, ip, powr, xpowr, ypowr, zpowr;
int iz, iend;

float one = 1.0;
float two = 2.0;

	     /*  initialize host structure  */

iprof = 0;
while(strcmp("HST",mdp[iprof].prof_name) != 0)
   {
   iprof++;
   if(iprof == nprof)
      {
      fprintf(stderr,"\n***********  ERROR  ***********\n");
      fprintf(stderr,"Host media, 'HST', does not have an DEF ");
      fprintf(stderr,"entry in model file.\n");
      exit(-1);
      }
   }

for(ix=0;ix<nx;ix++)
   {
   mptr = medf + ix;
   m2ptr = mdp[iprof].medptr;
   store(nz,m2ptr,     1,mptr,         nx);  /* lambda + 2*mu */
   store(nz,m2ptr+  nz,1,mptr+   nx*nz,nx);  /* lambda */
   store(nz,m2ptr+2*nz,1,mptr+ 4*nx*nz,nx);  /* 1/mu */
   store(nz,m2ptr+3*nz,1,mptr+ 7*nx*nz,nx);  /* rho */
   store(nz,m2ptr+4*nz,1,mptr+ 9*nx*nz,nx);  /* qp */
   store(nz,m2ptr+5*nz,1,mptr+10*nx*nz,nx);  /* qs */
   }

	     /*  get basin structure  */

for(ix=0;ix<nx;ix++)
   {
   if(strcmp("HST",ms[ix].med_name) != 0)
      {
      iprof = 0;
      izst = ms[ix].idep0;
      izend = ms[ix].idep1;
      if(izend > nz)
         izend = nz;

      while(strcmp(ms[ix].med_name,mdp[iprof].prof_name) != 0)
         {
         iprof++;
         if(iprof == nprof)
            {
            fprintf(stderr,"\n***********  ERROR  ***********\n");
            fprintf(stderr,"Media name '%s' ",ms[ix].med_name);
            fprintf(stderr,"does not have a DEF ");
            fprintf(stderr,"entry in model file.\n"); 
            exit(-1);
            }
         }

      if(izend > fs)
	 {
         mptr = medf + izst*nx + ix;
	 if(BASINPROF && izst == 0)        /* assume izst=0 */
            copybasinprof(mptr,nx*nz,&mdp[iprof],nz,izend,nx,fs);
	 else
	    {
	    m2ptr = mdp[iprof].medptr;
            store(izend-izst,m2ptr+izst,     1,mptr,         nx);
            store(izend-izst,m2ptr+izst+  nz,1,mptr+   nx*nz,nx);
            store(izend-izst,m2ptr+izst+2*nz,1,mptr+ 4*nx*nz,nx);
            store(izend-izst,m2ptr+izst+3*nz,1,mptr+ 7*nx*nz,nx);
            store(izend-izst,m2ptr+izst+4*nz,1,mptr+ 9*nx*nz,nx);
            store(izend-izst,m2ptr+izst+5*nz,1,mptr+10*nx*nz,nx);
	    }
	 }
      }
   }

/* Apply damping function near boundaries */

qp = medf + 9*nx*nz;
qs = medf + 10*nx*nz;

rfac = one/(dw - NBND_PAD - 1);
qbndf = (*qbmin)/(*qbmax);
qfac0 = exp(rfac*log(qbndf));

ypowr = 0;
if(iy < dw)
   ypowr = dw - iy;
else if(iy >= ny-dw)
   ypowr = iy - (ny - 1 - dw);

for(iz=0;iz<nz;iz++)
   {
   zpowr = 0;
   if(iz < dw && fs == 0)
      zpowr = dw - iz;
   else if(iz >= nz-dw)
      zpowr = iz - (nz - 1 - dw);

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iz*nx;

      xpowr = 0;
      if(ix < dw)
	 xpowr = dw - ix;
      else if(ix >= nx-dw)
	 xpowr = ix - (nx - 1 - dw);

      if(xpowr || ypowr || zpowr)
         {
         powr = xpowr + ypowr + zpowr - 1;
         qfac = one;
         while(powr--)
            qfac = qfac*qfac0;

	 qp0 = qp[ip];
	 if(qp0 > 2*(*qbmax))
	    qp0 = 2*(*qbmax);

	 qs0 = qs[ip];
	 if(qs0 > (*qbmax))
	    qs0 = (*qbmax);

	 qp[ip] = qfac*qp0;
	 qs[ip] = qfac*qs0;

	 if(qp[ip] < 5.0)
	    qp[ip] = 5.0;
	 if(qs[ip] < 5.0)
	    qs[ip] = 5.0;

         qp[ip] = -qp[ip];
         qs[ip] = -qs[ip];
         }
      }
   }

/*
   Store raw values of lambda and mu into medf array (11*nx*nz and
   12*nx*nz offsets).  These values will never be modified by averaging
   or attenuation (phase velocity adjustments).
*/

mptr1 = medf +    nx*nz;     /* lambda */
mptr2 = medf + 11*nx*nz;     /* lambda raw */
mptr3 = medf +  4*nx*nz;     /* 1/mu */
mptr4 = medf + 12*nx*nz;     /* mu raw */

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = ix+iz*nx;
      mptr2[ip] = mptr1[ip];
      mptr4[ip] = one/mptr3[ip];
      }
   }
}

getmedslice_modfac(ms,mdp,nprof,medf,nx,ny,nz,fs,iy,dw,qbnd,nt,dt,f0)
struct medstencil *ms;
struct medprof *mdp;
float *medf, *qbnd, *dt, *f0;
int nprof, nx, ny, nz, fs, iy, dw, nt;
{
float *mptr, *m2ptr;
int iprof, ix, izst, izend;

float *lam2mu, *mu, *lam;
float *qp, *qs, qp0, qs0, qfac;
float xfac, yfac, zfac, invdw;
int i, ip, powr, xpowr, ypowr, zpowr;
int iz, iend;

float taumin, taumax, w0, wmax, fmax, pi2;
float lg0, pfac, sfac;
float fre, fim, xre, xim;

float pi = 3.1415927;

float one = 1.0;
float two = 2.0;
float rfac = 40.0;  /* Q is reduced to Q/rfac over width of 'dw' */

pi2 = 2.0*pi;

taumin = (*dt)/pi;
taumax = 5*nt*taumin;

w0 = pi2*(*f0);
if((*f0) <= 0.0)
   {
   w0 = 1.0/sqrt(taumax*taumin);
   *f0 = w0/pi2;
   }

/*
fprintf(stderr,"f0=%13.5e\n",*f0);
*/

wmax = 1.0/taumin;
fmax = wmax/pi2;

lg0 = log(wmax/w0);

fre = w0*w0*taumin*taumin*taumax*taumax;
fre = log((taumax*taumax + fre)/(taumin*taumin + fre));

fim = w0*(taumax - taumin)/(1.0 + w0*w0*taumin*taumax);
fim = atan(fim);

	     /*  initialize host structure  */

iprof = 0;
while(strcmp("HST",mdp[iprof].prof_name) != 0)
   {
   iprof++;
   if(iprof == nprof)
      {
      fprintf(stderr,"\n***********  ERROR  ***********\n");
      fprintf(stderr,"Host media, 'HST', does not have an DEF ");
      fprintf(stderr,"entry in model file.\n");
      exit(-1);
      }
   }

for(ix=0;ix<nx;ix++)
   {
   mptr = medf + ix;
   m2ptr = mdp[iprof].medptr;
   store(nz,m2ptr,     1,mptr,         nx);  /* lambda + 2*mu */
   store(nz,m2ptr+  nz,1,mptr+   nx*nz,nx);  /* lambda */
   store(nz,m2ptr+2*nz,1,mptr+ 4*nx*nz,nx);  /* 1/mu */
   store(nz,m2ptr+3*nz,1,mptr+ 7*nx*nz,nx);  /* rho */
   store(nz,m2ptr+4*nz,1,mptr+ 9*nx*nz,nx);  /* qp */
   store(nz,m2ptr+5*nz,1,mptr+10*nx*nz,nx);  /* qs */
   }

	     /*  get basin structure  */

for(ix=0;ix<nx;ix++)
   {
   if(strcmp("HST",ms[ix].med_name) != 0)
      {
      iprof = 0;
      izst = ms[ix].idep0;
      izend = ms[ix].idep1;
      if(izend > nz)
         izend = nz;

      while(strcmp(ms[ix].med_name,mdp[iprof].prof_name) != 0)
         {
         iprof++;
         if(iprof == nprof)
            {
            fprintf(stderr,"\n***********  ERROR  ***********\n");
            fprintf(stderr,"Media name '%s' ",ms[ix].med_name);
            fprintf(stderr,"does not have a DEF ");
            fprintf(stderr,"entry in model file.\n"); 
            exit(-1);
            }
         }

      if(izend > fs)
	 {
         mptr = medf + izst*nx + ix;
	 if(BASINPROF && izst == 0)        /* assume izst=0 */
            copybasinprof(mptr,nx*nz,&mdp[iprof],nz,izend,nx,fs);
	 else
	    {
	    m2ptr = mdp[iprof].medptr;
            store(izend-izst,m2ptr+izst,     1,mptr,         nx);
            store(izend-izst,m2ptr+izst+  nz,1,mptr+   nx*nz,nx);
            store(izend-izst,m2ptr+izst+2*nz,1,mptr+ 4*nx*nz,nx);
            store(izend-izst,m2ptr+izst+3*nz,1,mptr+ 7*nx*nz,nx);
            store(izend-izst,m2ptr+izst+4*nz,1,mptr+ 9*nx*nz,nx);
            store(izend-izst,m2ptr+izst+5*nz,1,mptr+10*nx*nz,nx);
	    }
	 }
      }
   }

/* Apply damping function near boundaries */

qp = medf + 9*nx*nz;
qs = medf + 10*nx*nz;

powr = 1;                 /* linear taper */
powr = 2;                 /* parabolic taper */

invdw = one/(float)(dw);
rfac = one - one/rfac;

ypowr = 0;
if(iy < dw)
   ypowr = dw - iy;
else if(iy >= ny-dw)
   ypowr = iy - (ny - 1 - dw);

yfac = one;
if(ypowr)
   {
   for(i=0;i<powr;i++)
      yfac = yfac*ypowr*invdw;

   yfac = one - yfac*rfac;
   }

for(iz=0;iz<nz;iz++)
   {
   zpowr = 0;
   if(iz < dw && fs == 0)
      zpowr = dw - iz;
   else if(iz >= nz-dw)
      zpowr = iz - (nz - 1 - dw);

   zfac = one;
   if(zpowr)
      {
      for(i=0;i<powr;i++)
         zfac = zfac*zpowr*invdw;

      zfac = one - zfac*rfac;
      }

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iz*nx;

      xpowr = 0;
      if(ix < dw)
	 xpowr = dw - ix;
      else if(ix >= nx-dw)
	 xpowr = ix - (nx - 1 - dw);

      xfac = one;
      if(xpowr)
	 {
         for(i=0;i<powr;i++)
            xfac = xfac*xpowr*invdw;

         xfac = one - xfac*rfac;
	 }

      if(xfac < one || yfac < one || zfac < one)
	 {
	 qfac = xfac*yfac*zfac;

	 qp0 = qp[ip];
	 if(qp0 > 2*(*qbnd))
	    qp0 = 2*(*qbnd);

	 qs0 = qs[ip];
	 if(qs0 > (*qbnd))
	    qs0 = (*qbnd);

	 qp[ip] = qfac*qp0;
	 qs[ip] = qfac*qs0;
	 }

      if(qp[ip] <= 0.0 || qs[ip] <= 0.0)
         {
         fprintf(stderr,"Problem with zero/negative Q, exiting...\n");
         exit(-99);
         }
      }
   }

/* Apply phase velocity correction for attenuation dispersion */
/* Remember that mu needs to go to 1/mu */

lam2mu = medf;
lam    = medf + nx*nz;
mu     = medf + 4*nx*nz;
qp     = medf + 9*nx*nz;
qs     = medf + 10*nx*nz;

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = iz*nx + ix;

      pfac = exp(2.0*lg0*(atan(1.0/qp[ip])/pi));
      sfac = exp(2.0*lg0*(atan(1.0/qs[ip])/pi));

      xre = 1.0 - fre/(qp[ip]*pi);
      xim = 2.0*fim/(qs[ip]*pi);
      pfac = 1.0/sqrt(xre*xre + xim*xim);

      xre = 1.0 - fre/(qs[ip]*pi);
      xim = 2.0*fim/(qs[ip]*pi);
      sfac = 1.0/sqrt(xre*xre + xim*xim);

/*
fprintf(stderr,"%3d %3d   pf=%7.4f sf=%7.4f\n",ix,iz,pfac,sfac);
*/
pfac = 1.0;
sfac = 1.0;

      lam2mu[ip] = pfac*lam2mu[ip];
      mu[ip] = sfac*mu[ip];

      lam[ip] = lam2mu[ip] - two*mu[ip];
      mu[ip] = one/mu[ip];
      }
   }
}

copybasinprof(medf,off1,mdp,off2,iend,stride,fs)
struct medprof *mdp;
float *medf;
int off1, off2, iend, stride, fs;
{
float fj, *medp;
int nend, i, j, n1;
float fn1, fe1, inve1, sfac;

float half = 0.5;
float one = 1.0;

nend = mdp->mxdp;
medp = mdp->medptr;

   /* set first (fs+1) values to input depth profile without scaling depths

for(i=0;i<fs+1;i++)
   {
   medf[i*stride]        = medp[i];
   medf[i*stride+  off1] = medp[i+  off2];
   medf[i*stride+4*off1] = medp[i+2*off2];
   medf[i*stride+7*off1] = medp[i+3*off2];
   medf[i*stride+8*off1] = medp[i+4*off2];
   }
   */

   /* scale the remainder of the profile depths */
   /* to give sloping layer interfaces          */

fn1 = (float)(nend - 1);
fe1 = (float)(iend - 1);
sfac = sqrt(fn1/fe1);

for(i=0;i<iend;i++)
   {
   fj = ((float)(i))*sfac + half;
   j = (int)(fj);

   medf[i*stride]        = medp[j];
   medf[i*stride+  off1] = medp[j+  off2];
   medf[i*stride+4*off1] = medp[j+2*off2];
   medf[i*stride+7*off1] = medp[j+3*off2];
   medf[i*stride+8*off1] = medp[j+4*off2];
   }
}

/*
   avgmedslice() computes the effective media parameters for lambda+2*mu,
   lambda, mu, and bouyancy (1/rho) for use in the differencing
   function diff().  All averages are harmonic, see Zahradnik et al. (1993).

   The averaged values are stored in the array medf as follows:

      l2m = medf           -> for txx,tyy,tzz update
      lam = medf + 1*nx*nz -> for txx,tyy,tzz update
      mxy = medf + 2*nx*nz -> for txy update
      mxz = medf + 3*nx*nz -> for txz update
      myz = medf + 4*nx*nz -> for tyz update
      bx  = medf + 5*nx*nz -> for vx update
      by  = medf + 6*nx*nz -> for vy update
      bz  = medf + 7*nx*nz -> for vz update

   In the following, plane 2 is centered at the current location of
   iy, with plane 1 at iy-1 and plane 3 at iy+1.

   CASE 1)

   When the media values are defined to be centered at the normal stress
   node (iflag == NORMAL_STRESS_CENTERED), then the media boundaries go
   through the velocity nodes (see Graves, 1996).  The averaging is
   given by:

      mxy[i] = 4/(1/mu2[i] + 1/mu2[i+1]  + 1/mu3[i]    + 1/mu2[i+1])
      mxz[i] = 4/(1/mu2[i] + 1/mu2[i+1]  + 1/mu2[i+nx] + 1/mu2[i+1+nx])
      myz[i] = 4/(1/mu2[i] + 1/mu2[i+nx] + 1/mu3[i]    + 1/mu3[i+nx])

      bx[i] = 2/(rho2[i] + rho2[i+1])
      by[i] = 2/(rho2[i] + rho3[i])
      bz[i] = 2/(rho2[i] + rho2[i+nx])

   No averaging is needed for lambda+2*mu and lambda.  Also note that
   all averages are computed in the forward sense.

   CASE 2)

   When the media values are defined to be centered at the NULL node, that
   is the node on the opposite corner of the unit cell from the normal stress
   node (iflag == NULL_NODE_CENTERED), then the media boundaries go
   through the normal stress node.  The averaging is given by:

      x[i] = 8/(1/x2[i] + 1/x2[i-1] + 1/x1[i] + 1/x2[i-nx]
                + 1/x1[i-1] + 1/x2[i-1-nx] + 1/x1[i-nx] + 1/x1[i-1-nx])

      where x=l2m, lam.

      mxy[i] = 2/(1/mu2[i] + 1/mu2[i-nx])
      mxz[i] = 2/(1/mu2[i] + 1/mu1[i])
      myz[i] = 2/(1/mu2[i] + 1/mu2[i-1])

      bx[i] = 4/(rho2[i] + rho1[i]   + rho2[i-nx] + rho1[i-nx])
      by[i] = 4/(rho2[i] + rho2[i-1] + rho2[i-nx] + rho2[i-1-nx])
      bz[i] = 4/(rho2[i] + rho2[i-1] + rho1[i]    + rho1[i-1])

   Also note that all averages are computed in the reverse sense.

*/

avgmedslice(mdf0,mdf1,mdf2,mdf3,nx,nz,iflag)
float *mdf0, *mdf1, *mdf2, *mdf3;
int nx, nz, iflag;
{
float *l2m0, *lam0, *ak0, *qp0, *qs0;
float *mxy0, *mxz0, *myz0;
float *bx0, *by0, *bz0;
float *l2m1, *l2m2, *l2m3, *lam1, *lam2, *lam3, *ak2, *qp2, *qs2;
float *invmu1, *invmu2, *invmu3, *rho1, *rho2, *rho3;
int i, ix, iz, i1, inx, inx1;

int m1, mnx;
float inwgt, outwgt, sumwgt, fxp, fxi;

double one = 1.0;
double two = 2.0;
double four = 4.0;
double eight = 8.0;

l2m1   = mdf1;
lam1   = mdf1 +   nx*nz;
invmu1 = mdf1 + 4*nx*nz;
rho1   = mdf1 + 7*nx*nz;

l2m2   = mdf2;
lam2   = mdf2 +   nx*nz;
invmu2 = mdf2 + 4*nx*nz;
rho2   = mdf2 + 7*nx*nz;
ak2    = mdf2 + 8*nx*nz;
qp2    = mdf2 + 9*nx*nz;
qs2    = mdf2 + 10*nx*nz;

l2m3   = mdf3;
lam3   = mdf3 +   nx*nz;
invmu3 = mdf3 + 4*nx*nz;
rho3   = mdf3 + 7*nx*nz;

l2m0 = mdf0;
lam0 = mdf0 +   nx*nz;
mxy0 = mdf0 + 2*nx*nz;
mxz0 = mdf0 + 3*nx*nz;
myz0 = mdf0 + 4*nx*nz;
bx0  = mdf0 + 5*nx*nz;
by0  = mdf0 + 6*nx*nz;
bz0  = mdf0 + 7*nx*nz;
ak0  = mdf0 + 8*nx*nz;
qp0  = mdf0 + 9*nx*nz;
qs0  = mdf0 + 10*nx*nz;

inwgt = 1.0;
outwgt = 1.0;
fxp = 5.0;
fxi = 1.0/fxp;

sumwgt = inwgt + 4*outwgt;
inwgt = inwgt/sumwgt;
outwgt = outwgt/sumwgt;

if(iflag == NORMAL_STRESS_CENTERED)  /* do forward averaging */
   {
   for(iz=0;iz<nz-1;iz++)
      {
      for(ix=0;ix<nx-1;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i + nx;
         m1 = i - 1;
         mnx = i - nx;

         l2m0[i] = l2m2[i];
         lam0[i] = lam2[i];

/*
   RWG, 05/09/2006:

   The following is to ensure stability at the normal stress free-surface
   when there is an extreme lateral velocity contrast right at the
   surface.  ESG2006, Grenoble.
*/
	    /*
	       */
	 if(iz == 1 && ix > 0)
	    {
	    if(l2m2[i] > fxp*l2m2[m1] || l2m2[i] < fxi*l2m2[m1] ||
	       l2m2[i] > fxp*l2m2[i1] || l2m2[i] < fxi*l2m2[i1] ||
	       l2m2[i] > fxp*l2m1[ i] || l2m2[i] < fxi*l2m1[ i] ||
	       l2m2[i] > fxp*l2m3[ i] || l2m2[i] < fxi*l2m3[ i] ||
	       lam2[i] > fxp*lam2[m1] || lam2[i] < fxi*lam2[m1] ||
	       lam2[i] > fxp*lam2[i1] || lam2[i] < fxi*lam2[i1] ||
	       lam2[i] > fxp*lam1[ i] || lam2[i] < fxi*lam1[ i] ||
	       lam2[i] > fxp*lam3[ i] || lam2[i] < fxi*lam3[ i] ||
	       invmu2[i] > fxp*invmu2[m1] || invmu2[i] < fxi*invmu2[m1] ||
	       invmu2[i] > fxp*invmu2[i1] || invmu2[i] < fxi*invmu2[i1] ||
	       invmu2[i] > fxp*invmu1[ i] || invmu2[i] < fxi*invmu1[ i] ||
	       invmu2[i] > fxp*invmu3[ i] || invmu2[i] < fxi*invmu3[ i])
	       {
               l2m0[i] = inwgt/l2m2[i] + outwgt*(one/l2m2[m1] + one/l2m2[i1]
                            + one/l2m1[i] + one/l2m3[i]);
               l2m0[i] = one/l2m0[i];

               lam0[i] = inwgt/lam2[i] + outwgt*(one/lam2[m1] + one/lam2[i1]
                            + one/lam1[i] + one/lam3[i]);
               lam0[i] = one/lam0[i];

	       /* flag for 2nd order operators in tstepp() */
	       l2m0[nx] = -1.0;

/*
	       l2m0[nx] = 1.0;
*/
	       }
	    }

         mxy0[i] = four/(invmu2[i] + invmu2[i1] + invmu3[i] + invmu3[i1]);
         mxz0[i] = four/(invmu2[i] + invmu2[i1] + invmu2[inx] + invmu2[inx+1]);
         myz0[i] = four/(invmu2[i] + invmu2[inx] + invmu3[i] + invmu3[inx]);

         bx0[i] = two/(rho2[i] + rho2[i1]);
         by0[i] = two/(rho2[i] + rho3[i]);
         bz0[i] = two/(rho2[i] + rho2[inx]);

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }

      i = nx-1 + iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }

   for(ix=0;ix<nx;ix++)
      {
      i = ix + (nz-1)*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }
   }

else if(iflag == NULL_NODE_CENTERED)  /* do reverse averaging */
   {
   for(ix=0;ix<nx;ix++)
      {
      i = ix;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }

   for(iz=1;iz<nz;iz++)
      {
      i = iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];

/* NULL NODE 
*/
      for(ix=1;ix<nx;ix++)
         {
         i = ix + iz*nx;
         i1 = i - 1;
         inx = i - nx;
         inx1 = i - nx - 1;

         l2m0[i] = eight/(one/l2m2[i] + one/l2m2[i1] + one/l2m1[i]
		    + one/l2m2[inx] + one/l2m1[i1] + one/l2m2[inx1]
		       + one/l2m1[inx] + one/l2m1[inx1]);

         lam0[i] = eight/(one/lam2[i] + one/lam2[i1] + one/lam1[i]
		    + one/lam2[inx] + one/lam1[i1] + one/lam2[inx1]
		       + one/lam1[inx] + one/lam1[inx1]);

         mxy0[i] = two/(invmu2[i] + invmu2[inx]);
         mxz0[i] = two/(invmu2[i] + invmu1[i]);
         myz0[i] = two/(invmu2[i] + invmu2[i1]);

         bx0[i] = four/(rho2[i] + rho1[i] + rho2[inx] + rho1[inx]);
         by0[i] = four/(rho2[i] + rho2[i1] + rho2[inx] + rho2[inx1]);
         bz0[i] = four/(rho2[i] + rho2[i1] + rho1[i] + rho1[i1]);

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }

/* vz NODE
      for(ix=1;ix<nx;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i - nx;
         inx1 = i - nx + 1;

         l2m0[i] = two/(one/l2m2[i] + one/l2m2[inx]);
         lam0[i] = two/(one/lam2[i] + one/lam2[inx]);

         mxy0[i] = eight/(invmu2[i] + invmu2[inx] + invmu3[i]
		    + invmu3[inx] + invmu2[i1] + invmu2[inx1]
		       + invmu3[i1] + invmu3[inx1]);

         mxz0[i] = two/(invmu2[i] + invmu2[i1]);
         myz0[i] = two/(invmu2[i] + invmu3[i]);

         bx0[i] = four/(rho2[i] + rho2[inx] + rho2[i1] + rho2[inx1]);
         by0[i] = four/(rho2[i] + rho2[inx] + rho3[inx] + rho3[i]);
         bz0[i] = one/rho2[i];

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }
*/

      }
   }

else if(iflag == KBO_AT_VXNODE)  /* do forward arithmetic averaging centerd at vx */
   {
   for(iz=0;iz<nz-1;iz++)
      {
      for(ix=0;ix<nx-1;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i + nx;
         m1 = i - 1;
         mnx = i - nx;

         l2m0[i] = 0.5*(l2m2[i] + l2m2[m1]);
         lam0[i] = 0.5*(lam2[i] + lam2[m1]);

         mxy0[i] = 0.5*(one/invmu2[i] + one/invmu3[i]);
         mxz0[i] = 0.5*(one/invmu2[i] + one/invmu2[inx]);
         myz0[i] = 0.125*(one/invmu2[i] + one/invmu2[inx] + one/invmu2[m1] + one/invmu2[inx-1] + one/invmu3[i] + one/invmu3[inx] + one/invmu3[m1] + one/invmu3[inx-1]);

         bx0[i] = one/rho2[i];
         by0[i] = four/(rho2[i] + rho2[m1] + rho3[i] + rho3[m1]);
         bz0[i] = four/(rho2[i] + rho2[m1] + rho2[inx] + rho2[inx-1]);

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }

      i = nx-1 + iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }

   for(ix=0;ix<nx;ix++)
      {
      i = ix + (nz-1)*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }
   }

else if(iflag == KBO_AT_VZNODE)  /* do forward arithmetic averaging centerd at vz */
   {
   for(ix=0;ix<nx;ix++)
      {
      i = ix;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }

   for(iz=1;iz<nz;iz++)
      {
      for(ix=0;ix<nx-1;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i + nx;
         m1 = i - 1;
         mnx = i - nx;

         l2m0[i] = 0.5*(l2m2[i] + l2m2[mnx]);
         lam0[i] = 0.5*(lam2[i] + lam2[mnx]);

         mxy0[i] = 0.125*(one/invmu2[i] + one/invmu2[mnx] + one/invmu2[i1] + one/invmu2[mnx+1] + one/invmu3[i] + one/invmu3[mnx] + one/invmu3[i1] + one/invmu3[mnx+1]);
         mxz0[i] = 0.5*(one/invmu2[i] + one/invmu2[i1]);
         myz0[i] = 0.5*(one/invmu2[i] + one/invmu3[i]);

         bx0[i] = four/(rho2[i] + rho2[mnx] + rho2[i1] + rho2[mnx+1]);
         by0[i] = four/(rho2[i] + rho2[mnx] + rho3[i] + rho3[mnx]);
         bx0[i] = 0.25*(one/rho2[i] + one/rho2[mnx] + one/rho2[i1] + one/rho2[mnx+1]);
         by0[i] = 0.25*(one/rho2[i] + one/rho2[mnx] + one/rho3[i] + one/rho3[mnx]);
         bz0[i] = one/rho2[i];

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }

      i = nx-1 + iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }
   }

else if(iflag == HARMONIC_AT_VZNODE)  /* do harmonic averaging centerd at vz */
   {
   for(ix=0;ix<nx;ix++)
      {
      i = ix;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }

   for(iz=1;iz<nz;iz++)
      {
      for(ix=0;ix<nx-1;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i + nx;
         m1 = i - 1;
         mnx = i - nx;

         l2m0[i] = 2.0/(one/l2m2[i] + one/l2m2[mnx]);
         lam0[i] = 2.0/(one/lam2[i] + one/lam2[mnx]);

         mxy0[i] = 8.000/(invmu2[i] + invmu2[mnx] + invmu2[i1] + invmu2[mnx+1] + invmu3[i] + invmu3[mnx] + invmu3[i1] + invmu3[mnx+1]);
         mxz0[i] = 2.0/(invmu2[i] + invmu2[i1]);
         myz0[i] = 2.0/(invmu2[i] + invmu3[i]);

         bx0[i] = four/(rho2[i] + rho2[mnx] + rho2[i1] + rho2[mnx+1]);
         by0[i] = four/(rho2[i] + rho2[mnx] + rho3[i] + rho3[mnx]);
         bz0[i] = one/rho2[i];

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }

      i = nx-1 + iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      ak0[i] = ak2[i];
      qp0[i] = qp2[i];
      qs0[i] = qs2[i];
      }
   }

else if(iflag == -999)  /* don't do averaging */
   {
   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
         i = ix + iz*nx;

         l2m0[i] = l2m2[i];
         lam0[i] = lam2[i];

         mxy0[i] = one/invmu2[i];
         mxz0[i] = one/invmu2[i];
         myz0[i] = one/invmu2[i];

         bx0[i] = one/rho2[i];
         by0[i] = one/rho2[i];
         bz0[i] = one/rho2[i];

         ak0[i] = ak2[i];
         qp0[i] = qp2[i];
         qs0[i] = qs2[i];
         }
      }
   }
}

avgmedsliceBBBB(mdf0,mdf1,mdf2,mdf3,nx,nz,iflag)
float *mdf0, *mdf1, *mdf2, *mdf3;
int nx, nz, iflag;
{
float *l2m0, *lam0, *qs0;
float *mxy0, *mxz0, *myz0;
float *bx0, *by0, *bz0;
float *l2m1, *l2m2, *l2m3, *lam1, *lam2, *lam3, *qs2;
float *invmu1, *invmu2, *invmu3, *rho1, *rho2, *rho3;
int i, ix, iz, i1, inx, inx1;

double half = 0.5;
double one = 1.0;
double two = 2.0;
double three = 3.0;
double four = 4.0;
double eight = 8.0;

l2m1   = mdf1;
lam1   = mdf1 +   nx*nz;
invmu1 = mdf1 + 4*nx*nz;
rho1   = mdf1 + 7*nx*nz;

l2m2   = mdf2;
lam2   = mdf2 +   nx*nz;
invmu2 = mdf2 + 4*nx*nz;
rho2   = mdf2 + 7*nx*nz;
qs2    = mdf2 + 8*nx*nz;

l2m3   = mdf3;
lam3   = mdf3 +   nx*nz;
invmu3 = mdf3 + 4*nx*nz;
rho3   = mdf3 + 7*nx*nz;

l2m0 = mdf0;
lam0 = mdf0 +   nx*nz;
mxy0 = mdf0 + 2*nx*nz;
mxz0 = mdf0 + 3*nx*nz;
myz0 = mdf0 + 4*nx*nz;
bx0  = mdf0 + 5*nx*nz;
by0  = mdf0 + 6*nx*nz;
bz0  = mdf0 + 7*nx*nz;
qs0  = mdf0 + 8*nx*nz;

if(iflag == NORMAL_STRESS_CENTERED)  /* do forward averaging */
   {
   for(iz=0;iz<nz-1;iz++)
      {
      for(ix=0;ix<nx-1;ix++)
         {
         i = ix + iz*nx;
         i1 = i + 1;
         inx = i + nx;

         l2m0[i] = l2m2[i];
         lam0[i] = lam2[i];

         mxy0[i] = four/(invmu2[i] + invmu2[i1] + invmu3[i] + invmu3[i1]);
         mxz0[i] = four/(invmu2[i] + invmu2[i1] + invmu2[inx] + invmu2[inx+1]);
         myz0[i] = four/(invmu2[i] + invmu2[inx] + invmu3[i] + invmu3[inx]);

         bx0[i] = two/(rho2[i] + rho2[i1]);
         by0[i] = two/(rho2[i] + rho3[i]);
         bz0[i] = two/(rho2[i] + rho2[inx]);

         qs0[i] = qs2[i];
         }

      i = nx-1 + iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      qs0[i] = qs2[i];
      }

   for(ix=0;ix<nx;ix++)
      {
      i = ix + (nz-1)*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      qs0[i] = qs2[i];
      }

/* GGGGGGGGGGGGGg */

   for(iz=1;iz<nz-2;iz++)
      {
      for(ix=1;ix<nx-2;ix++)
         {
         i = ix + iz*nx;

         mxz0[i] = three/(0.5*invmu2[i-nx]+invmu2[i]+invmu2[i+nx]+0.5*invmu2[i+2*nx]);
         myz0[i] = three/(0.5*invmu2[i-nx]+invmu2[i]+invmu2[i+nx]+0.5*invmu2[i+2*nx]);

         bz0[i] = three/(0.5*rho2[i-nx]+rho2[i]+rho2[i+nx]+0.5*rho2[i+2*nx]);
         }
      }
   }

else if(iflag == NULL_NODE_CENTERED)  /* do reverse averaging */
   {
   for(ix=0;ix<nx;ix++)
      {
      i = ix;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      qs0[i] = qs2[i];
      }

   for(iz=1;iz<nz;iz++)
      {
      i = iz*nx;

      l2m0[i] = l2m2[i];
      lam0[i] = lam2[i];

      mxy0[i] = one/invmu2[i];
      mxz0[i] = one/invmu2[i];
      myz0[i] = one/invmu2[i];

      bx0[i] = one/rho2[i];
      by0[i] = one/rho2[i];
      bz0[i] = one/rho2[i];

      qs0[i] = qs2[i];

      for(ix=1;ix<nx;ix++)
         {
         i = ix + iz*nx;
         i1 = i - 1;
         inx = i - nx;
         inx1 = i - nx - 1;

         l2m0[i] = eight/(one/l2m2[i] + one/l2m2[i1] + one/l2m1[i]
		    + one/l2m2[inx] + one/l2m1[i1] + one/l2m2[inx1]
		       + one/l2m1[inx] + one/l2m1[inx1]);

         lam0[i] = eight/(one/lam2[i] + one/lam2[i1] + one/lam1[i]
		    + one/lam2[inx] + one/lam1[i1] + one/lam2[inx1]
		       + one/lam1[inx] + one/lam1[inx1]);

         mxy0[i] = two/(invmu2[i] + invmu2[inx]);
         mxz0[i] = two/(invmu2[i] + invmu1[i]);
         myz0[i] = two/(invmu2[i] + invmu2[i1]);

         bx0[i] = four/(rho2[i] + rho1[i] + rho2[inx] + rho1[inx]);
         by0[i] = four/(rho2[i] + rho2[i1] + rho2[inx] + rho2[inx1]);
         bz0[i] = four/(rho2[i] + rho2[i1] + rho1[i] + rho1[i1]);

         qs0[i] = qs2[i];
         }
      }

/* GGGGGGGGGGGGGg */

   for(iz=2;iz<nz-1;iz++)
      {
      for(ix=1;ix<nx-1;ix++)
         {
         i = ix + iz*nx;

         l2m0[i] = three/(0.5/l2m2[i-2*nx] + one/l2m2[i-nx]
		   + one/l2m2[i] + 0.5/l2m2[i+nx]);

         lam0[i] = three/(0.5/lam2[i-2*nx] + one/lam2[i-nx]
		   + one/lam2[i] + 0.5/lam2[i+nx]);

         mxy0[i] = three/(0.5*invmu2[i-2*nx] + invmu2[i-nx]
		   + invmu2[i] + 0.5*invmu2[i+nx]);

         bx0[i] = three/(0.5*rho2[i-2*nx] + rho2[i-nx]
		   + rho2[i] + 0.5*rho2[i+nx]);
         by0[i] = three/(0.5*rho2[i-2*nx] + rho2[i-nx]
		   + rho2[i] + 0.5*rho2[i+nx]);
         }
      }
   }
}

copy_medslice(mdf1,mdf2,nx,nz)
float *mdf1, *mdf2;
int nx, nz;
{
float *l2m1, *lam1, *mu1, *boy1, *qp1, *qs1;
float *l2m2, *lam2, *mu2, *boy2, *qp2, *qs2;
float *lraw1, *lraw2;
float *mraw1, *mraw2;
int ix, iz;

l2m1 = mdf1;
lam1 = mdf1 +   nx*nz;
mu1  = mdf1 + 4*nx*nz;
boy1 = mdf1 + 7*nx*nz;
qp1  = mdf1 + 9*nx*nz;
qs1  = mdf1 + 10*nx*nz;
lraw1  = mdf1 + 11*nx*nz;
mraw1  = mdf1 + 12*nx*nz;

l2m2 = mdf2;
lam2 = mdf2 +   nx*nz;
mu2  = mdf2 + 4*nx*nz;
boy2 = mdf2 + 7*nx*nz;
qp2  = mdf2 + 9*nx*nz;
qs2  = mdf2 + 10*nx*nz;
lraw2  = mdf2 + 11*nx*nz;
mraw2  = mdf2 + 12*nx*nz;

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      l2m2[iz*nx + ix] = l2m1[iz*nx + ix];
      lam2[iz*nx + ix] = lam1[iz*nx + ix];
      mu2[iz*nx + ix]  = mu1[iz*nx + ix];
      boy2[iz*nx + ix] = boy1[iz*nx + ix];
      qp2[iz*nx + ix]  = qp1[iz*nx + ix];
      qs2[iz*nx + ix]  = qs1[iz*nx + ix];
      lraw2[iz*nx + ix]  = lraw1[iz*nx + ix];
      mraw2[iz*nx + ix]  = mraw1[iz*nx + ix];
      }
   }
}

copy_avgmedslice(mdf1,mdf2,nx,nz)
float *mdf1, *mdf2;
int nx, nz;
{
float *l2m1, *lam1, *mu1x, *mu1y, *mu1z, *boy1x, *boy1y, *boy1z, *ak1, *qp1, *qs1;
float *l2m2, *lam2, *mu2x, *mu2y, *mu2z, *boy2x, *boy2y, *boy2z, *ak2, *qp2, *qs2;
float *lraw1, *lraw2;
float *mraw1, *mraw2;
int ix, iz;

l2m1  = mdf1;
lam1  = mdf1 +   nx*nz;
mu1x  = mdf1 + 2*nx*nz;
mu1y  = mdf1 + 3*nx*nz;
mu1z  = mdf1 + 4*nx*nz;
boy1x = mdf1 + 5*nx*nz;
boy1y = mdf1 + 6*nx*nz;
boy1z = mdf1 + 7*nx*nz;
ak1   = mdf1 + 8*nx*nz;
qp1   = mdf1 + 9*nx*nz;
qs1   = mdf1 + 10*nx*nz;
lraw1   = mdf1 + 11*nx*nz;
mraw1   = mdf1 + 12*nx*nz;

l2m2  = mdf2;
lam2  = mdf2 +   nx*nz;
mu2x  = mdf2 + 2*nx*nz;
mu2y  = mdf2 + 3*nx*nz;
mu2z  = mdf2 + 4*nx*nz;
boy2x = mdf2 + 5*nx*nz;
boy2y = mdf2 + 6*nx*nz;
boy2z = mdf2 + 7*nx*nz;
ak2   = mdf2 + 8*nx*nz;
qp2   = mdf2 + 9*nx*nz;
qs2   = mdf2 + 10*nx*nz;
lraw2   = mdf2 + 11*nx*nz;
mraw2   = mdf2 + 12*nx*nz;

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      l2m2[iz*nx + ix]  = l2m1[iz*nx + ix];
      lam2[iz*nx + ix]  = lam1[iz*nx + ix];
      mu2x[iz*nx + ix]  = mu1x[iz*nx + ix];
      mu2y[iz*nx + ix]  = mu1y[iz*nx + ix];
      mu2z[iz*nx + ix]  = mu1z[iz*nx + ix];
      boy2x[iz*nx + ix] = boy1x[iz*nx + ix];
      boy2y[iz*nx + ix] = boy1y[iz*nx + ix];
      boy2z[iz*nx + ix] = boy1z[iz*nx + ix];
      ak2[iz*nx + ix]   = ak1[iz*nx + ix];
      qp2[iz*nx + ix]   = qp1[iz*nx + ix];
      qs2[iz*nx + ix]   = qs1[iz*nx + ix];
      lraw2[iz*nx + ix]   = lraw1[iz*nx + ix];
      mraw2[iz*nx + ix]   = mraw1[iz*nx + ix];
      }
   }
}

xzpad_medslice(mdf,nx,nz,npad,fs)
float *mdf;
int nx, nz, npad, fs;
{
float *l2m, *lam, *mu, *boy, *qp, *qs;
float *lraw, *mraw;
int ix, iz;

l2m = mdf;
lam = mdf +   nx*nz;
mu  = mdf + 4*nx*nz;
boy = mdf + 7*nx*nz;
qp  = mdf + 9*nx*nz;
qs  = mdf + 10*nx*nz;
lraw = mdf + 11*nx*nz;
mraw = mdf + 12*nx*nz;

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<npad;ix++)
      {
      l2m[iz*nx + ix] = l2m[iz*nx + npad];
      lam[iz*nx + ix] = lam[iz*nx + npad];
      mu[iz*nx + ix]  = mu[iz*nx + npad];
      boy[iz*nx + ix] = boy[iz*nx + npad];
      qp[iz*nx + ix]  = qp[iz*nx + npad];
      qs[iz*nx + ix]  = qs[iz*nx + npad];
      lraw[iz*nx + ix]  = lraw[iz*nx + npad];
      mraw[iz*nx + ix]  = mraw[iz*nx + npad];
      }
   for(ix=nx-npad;ix<nx;ix++)
      {
      l2m[iz*nx + ix] = l2m[iz*nx + (nx-npad-1)];
      lam[iz*nx + ix] = lam[iz*nx + (nx-npad-1)];
      mu[iz*nx + ix]  = mu[iz*nx + (nx-npad-1)];
      boy[iz*nx + ix] = boy[iz*nx + (nx-npad-1)];
      qp[iz*nx + ix]  = qp[iz*nx + (nx-npad-1)];
      qs[iz*nx + ix]  = qs[iz*nx + (nx-npad-1)];
      lraw[iz*nx + ix]  = lraw[iz*nx + (nx-npad-1)];
      mraw[iz*nx + ix]  = mraw[iz*nx + (nx-npad-1)];
      }
   }

if(!fs)
   {
   for(iz=0;iz<npad;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
         l2m[iz*nx + ix] = l2m[npad*nx + ix];
         lam[iz*nx + ix] = lam[npad*nx + ix];
         mu[iz*nx + ix]  = mu[npad*nx + ix];
         boy[iz*nx + ix] = boy[npad*nx + ix];
         qp[iz*nx + ix]  = qp[npad*nx + ix];
         qs[iz*nx + ix]  = qs[npad*nx + ix];
         lraw[iz*nx + ix]  = lraw[npad*nx + ix];
         mraw[iz*nx + ix]  = mraw[npad*nx + ix];
         }
      }
   }
for(iz=nz-npad;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      l2m[iz*nx + ix] = l2m[(nz-npad-1)*nx + ix];
      lam[iz*nx + ix] = lam[(nz-npad-1)*nx + ix];
      mu[iz*nx + ix]  = mu[(nz-npad-1)*nx + ix];
      boy[iz*nx + ix] = boy[(nz-npad-1)*nx + ix];
      qp[iz*nx + ix]  = qp[(nz-npad-1)*nx + ix];
      qs[iz*nx + ix]  = qs[(nz-npad-1)*nx + ix];
      lraw[iz*nx + ix]  = lraw[(nz-npad-1)*nx + ix];
      mraw[iz*nx + ix]  = mraw[(nz-npad-1)*nx + ix];
      }
   }
}

char *getlineRWG(fp,linenum,ungetflag,buf)
FILE *fp;
char *buf;
int *linenum, *ungetflag;
{
char *ptr, c;

if((*ungetflag))
   {
   *ungetflag = 0;
   return(buf);
   }
ptr = buf;
while((c = getc(fp)) != EOF && !feof(fp))
   {
   if(c == '\n')
      {
      *ptr = '\0';
      (*linenum)++;
      return(buf);
      }
   *ptr++ = c;
   }
return(NULL);
}

getwords(str,w)
char *str, *w[];
{
int nw;

nw = 0;
while(*str != '\0')
   {
   while(*str == ' ' || *str == '\t')
      str++;
   if(*str == '\0')
      break;

   w[nw] = str;
   nw++;
   while(*str != ' ' && *str != '\t' && *str != '\0')
      str++;
   if(*str != '\0')
      *str++ = '\0';
   }
return(nw);
}

init_name(name,n)
char *name;
int n;
{
int i;

for(i=0;i<n;i++)
   name[i] = '\0';
}

modsten2mplot(nx,ny,nz,ms,name,fs,h)
struct medstencil *ms;
char *name;
float h;
int nx, ny, nz, fs;
{
FILE *fpw, *fopfile();
char string[128];
int ix, iy, ifs;
float dep_val;

ifs = 0;
if(fs)
   ifs = 1;

sprintf(string,"%s.mplot",name);
fpw = fopfile(string,"w");
fprintf(fpw,"%d %d\n",nx,ny);

for(iy=0;iy<ny;iy++)
   {
   for(ix=0;ix<nx;ix++)
      {
      if(ms[(ny-1-iy)*nx + ix].idep1 == nz &&
         (strcmp("HST",ms[(ny-1-iy)*nx + ix].med_name) == 0))
         dep_val = 0.0;
      else
         dep_val = (float)((ms[(ny-1-iy)*nx + ix].idep1-ifs)*h);
   
      fprintf(fpw," %10.5f\n",dep_val);
      }
   }
fclose(fpw);
}

diffmedslice(iy,med1,med2,nx,nz)
float *med1, *med2;
int iy, nx, nz;
{
float *lam2mu1, *mu1, *lam1, *rho1, *qs1;
float *lam2mu2, *mu2, *lam2, *rho2, *qs2;
float diff;
int i, ix, iz, nn;

float del = 0.01;
float hun = 100.0;

nn = nx*nz;

lam2mu1 = med1;
lam1    = med1 + nn;
mu1     = med1 + 4*nn;
rho1    = med1 + 7*nn;
qs1     = med1 + 8*nn;

lam2mu2 = med2;
lam2    = med2 + nn;
mu2     = med2 + 4*nn;
rho2    = med2 + 7*nn;
qs2     = med2 + 8*nn;

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      i = iz*nx + ix;

      diff = hun*(lam2mu1[i] - lam2mu2[i])/lam2mu1[i];
      printdiff(ix,iy,iz,&diff,&del,"lam2mu");

      diff = hun*(lam1[i] - lam2[i])/lam1[i];
      printdiff(ix,iy,iz,&diff,&del,"lam");

      diff = hun*(mu1[i] - mu2[i])/mu1[i];
      printdiff(ix,iy,iz,&diff,&del,"mu");

      diff = hun*(rho1[i] - rho2[i])/rho1[i];
      printdiff(ix,iy,iz,&diff,&del,"rho");

      diff = hun*(log(qs1[i]) - log(qs2[i]))/log(qs1[i]);
      printdiff(ix,iy,iz,&diff,&del,"qs");
      }
   }
}

printdiff(ix,iy,iz,diff,del,name)
float *diff, *del;
int ix, iy, iz;
char *name;
{
float zap = 0.0;

if(*diff < zap)
   *diff = -(*diff);
      
if(*diff > *del)
   {
   fprintf(stderr,"*** Media model mismatch (%s):\n",name);
   fprintf(stderr,"    ix=%3d iy=%3d iz=%3d diff=%8.3g tolerance=%8.3g\n",ix,iy,iz,*diff,*del);
   }
}

read_qvals(qfile,qv)
char *qfile;
struct qvalues *qv;
{
FILE *fpr, *fopfile();
int i;

fpr = fopfile(qfile,"r");

fscanf(fpr,"%d",&qv->n);

qv->vp = (float *) check_malloc ((qv->n)*sizeof(float));
qv->qp = (float *) check_malloc ((qv->n)*sizeof(float));
qv->vs = (float *) check_malloc ((qv->n)*sizeof(float));
qv->qs = (float *) check_malloc ((qv->n)*sizeof(float));
for(i=0;i<qv->n;i++)
   {
   fscanf(fpr,"%f %f %f %f",&qv->vp[i],&qv->qp[i],&qv->vs[i],&qv->qs[i]);

   if(qv->qp[i] <= 0.0)
      qv->qp[i] = 1.0e+10;
   if(qv->qs[i] <= 0.0)
      qv->qs[i] = 1.0e+10;
   }

fclose(fpr);
}

void init_media(struct media_input *mi,float *mfield,struct runparams *rpars,struct qvalues *qv,struct modelstorage *mstor,int swapb)
{
int nodeType, nx, ny, ny1, ny2, nz, fs;
float h;

struct medstencil *ms;
struct medprof *mp;
float *medf;
int i, iy, iymin, iymax, print_tol, ptol;
int ix, iz;
int fdp, fds, fdd;
int nbpad = NBND_PAD;
off_t boff;

nx = rpars->nx;
ny = rpars->globny;
nz = rpars->nz;
ny1 = rpars->ny1;
ny2 = rpars->ny2;
fs = rpars->freesurf;
nodeType = rpars->nodeType;
h = rpars->h;

print_tol = 25;
ptol = print_tol;

init_model_seek(mstor,MEDIA_FIELD);

if(mi->model_style == 0)  /* generate model from stencil */
   {
   mi->medprofs = (struct medprof *) check_malloc (NPROFS*sizeof(struct medprof));
   mi->medsten = (struct medstencil *) check_malloc (nx*ny*sizeof(struct medstencil));

   ms = mi->medsten;
   mp = mi->medprofs;

   for(i=0;i<NPROFS;i++)
      mp[i].medptr = (float *) check_malloc (N_MEDPROF_VARS*nz*sizeof(float));

   genmod3d(&mi->nprof,mi->modfile,&mi->xorg,&mi->yorg,nx,ny,nz,&h,ms,mp,fs,mi->media_boundary);
   }

else if(mi->model_style == 1) /* read in binary velocity model from disk */
   {
   fdp = opfile_ro(mi->pmodfile);
   fds = opfile_ro(mi->smodfile);
   fdd = opfile_ro(mi->dmodfile);

   if(qv->qpfrac <= 0.0 || qv->qsfrac <= 0.0)
      read_qvals(mi->qmodfile,qv);
   }

/*
   Add padding to iy=0,ny-1 ends of model to ensure stability of absorbing
   boundaries.  There should be no media variations in the direction
   perpendicular to the absorbing boundary.  See also, xzpad_medslice()
   for padding of x and z boundaries.
*/

medf = mfield;

if(nodeType == PARL_HEAD)
   {
   if(mi->model_style == 0)
      getmedslice(ms+nbpad*nx,mp,mi->nprof,medf,nx,ny,nz,fs,nbpad,mi->dampwidth,&mi->qbndmax,&mi->qbndmin);

   else if(mi->model_style == 1)
      {
      boff = (off_t)(nbpad)*(off_t)(nx)*(off_t)(nz*sizeof(float));
      lseek(fdp,boff,SEEK_SET);
      lseek(fds,boff,SEEK_SET);
      lseek(fdd,boff,SEEK_SET);

      reedmedslice(fdp,fds,fdd,medf,nx,ny,nz,fs,qv,nbpad,mi->dampwidth,&mi->qbndmax,&mi->qbndmin,swapb);
      }

   xzpad_medslice(medf,nx,nz,nbpad,fs);

   for(iy=0;iy<=nbpad;iy++)
      {
      if(mstor->intmem)
         {
         if(iy)
            copy_medslice(medf,mfield+iy*N_MED_VARS*nx*nz,nx,nz);
         }
      else
         rite_model(mstor,medf,MEDIA_FIELD);
      }
   iymin = nbpad + 1;
   }
else
   {
   iymin = ny1;

   if(mi->model_style == 1)
      {
      boff = (off_t)(iymin)*(off_t)(nx)*(off_t)(nz*sizeof(float));
      lseek(fdp,boff,SEEK_SET);
      lseek(fds,boff,SEEK_SET);
      lseek(fdd,boff,SEEK_SET);
      }
   }

iymax = ny - nbpad;
if(iymax > ny2)
   iymax = ny2;

for(iy=iymin;iy<iymax;iy++)
   {
   /*SN: if iy>=ny1, then it is model space for this machine */

   if(mstor->intmem)
      medf = mfield + (iy-ny1)*N_MED_VARS*nx*nz;

   if(mi->model_style == 0)
      getmedslice(ms+iy*nx,mp,mi->nprof,medf,nx,ny,nz,fs,iy,mi->dampwidth,&mi->qbndmax,&mi->qbndmin);
   else if(mi->model_style == 1)
      reedmedslice(fdp,fds,fdd,medf,nx,ny,nz,fs,qv,iy,mi->dampwidth,&mi->qbndmax,&mi->qbndmin,swapb);

   xzpad_medslice(medf,nx,nz,nbpad,fs);
   rite_model(mstor,medf,MEDIA_FIELD);

   if((float)(100.0*(iy-iymin))/(float)(iymax-1-iymin) >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy-iymin,iymax-1-iymin);
      ptol = ptol + print_tol;
      }
   }

if(nodeType == PARL_TAIL)
   {
   for(iy=(ny-nbpad);iy<ny;iy++)
      {
      if(mstor->intmem)
         copy_medslice(medf,mfield+(iy-ny1)*N_MED_VARS*nx*nz,nx,nz);
      else
         rite_model(mstor,medf,MEDIA_FIELD);
      }
   }

if(mi->model_style == 1)
   {
   close(fdp);
   close(fds);
   close(fdd);
   }

init_model_seek(mstor,MEDIA_FIELD);
}

int eff_media(struct media_input *mi,float *mfield,int nx,int ny,int ny1,int ny2,int nz,struct modelstorage *mstor,int nodeType)
{
float *medbuf, *medf[4], *tmp_ptr;
int iy, print_tol, ptol;
int ix, iz;
int extraPlane;

ny = ny2 - ny1;

init_model_seek(mstor,MEDIA_FIELD);

medf[0] = mfield;

medbuf = (float *) check_malloc (3*N_MED_VARS*nx*nz*sizeof(float));
medf[1] = medbuf;
medf[2] = medbuf + N_MED_VARS*nx*nz;
medf[3] = medbuf + 2*N_MED_VARS*nx*nz;

if(mstor->intmem)
   {
   copy_medslice(mfield,medf[1],nx,nz);
   copy_medslice(mfield+N_MED_VARS*nx*nz,medf[2],nx,nz);
   }
else
   {
   medf[1] = reed_model(mstor,medf[1],MEDIA_FIELD);
   medf[2] = reed_model(mstor,medf[2],MEDIA_FIELD);
   }

/* replicate first plane */
copy_medslice(medf[2],medf[3],nx,nz);
avgmedslice(medf[0],medf[1],medf[2],medf[3],nx,nz,mi->media_boundary);

rite_model(mstor,medf[0],MEDIA_FIELD);

print_tol = 25;
ptol = print_tol;

for(iy=1;iy<ny-1;iy++)
   {
   if(mstor->intmem)
      {
      medf[0] = mfield + iy*N_MED_VARS*nx*nz;
      copy_medslice(mfield+(iy+1)*N_MED_VARS*nx*nz,medf[3],nx,nz);
      }
   else
      medf[3] = reed_model(mstor,medf[3],MEDIA_FIELD);

   avgmedslice(medf[0],medf[1],medf[2],medf[3],nx,nz,mi->media_boundary);

   rite_model(mstor,medf[0],MEDIA_FIELD);

   tmp_ptr = medf[1];
   medf[1] = medf[2];
   medf[2] = medf[3];
   medf[3] = tmp_ptr;

   if((float)(100.0*(iy-1))/(float)(ny-3) >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy-1,ny-3);
      ptol = ptol + print_tol;
      }
   }

/* replicate last plane */
if(mstor->intmem)
   {
   medf[0] = mfield + (ny-1)*N_MED_VARS*nx*nz;
   copy_avgmedslice(mfield+(ny-2)*N_MED_VARS*nx*nz,medf[0],nx,nz);
   }
else
   rite_model(mstor,medf[0],MEDIA_FIELD);

free(medbuf);
init_model_seek(mstor,MEDIA_FIELD);
return(1);
}

int eff_mediaP3(struct media_input *mi,float *mfield,struct runparamsP3 *rpars)
{
float *medbuf, *medf[4], *tmp_ptr;
int iy, print_tol, ptol;
int nx, nz, ny;

struct nodeinfo *ni;

ni = &(rpars->ni);

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

medf[0] = mfield;

medbuf = (float *) check_malloc (3*N_MED_VARS*nx*nz*sizeof(float));
medf[1] = medbuf;
medf[2] = medbuf + N_MED_VARS*nx*nz;
medf[3] = medbuf + 2*N_MED_VARS*nx*nz;

copy_medslice(mfield,medf[1],nx,nz);
copy_medslice(mfield+N_MED_VARS*nx*nz,medf[2],nx,nz);

print_tol = 25;
ptol = print_tol;

for(iy=1;iy<ny-1;iy++)
   {
   /* get next plane */
   copy_medslice(mfield+(iy+1)*N_MED_VARS*nx*nz,medf[3],nx,nz);

   medf[0] = mfield + iy*N_MED_VARS*nx*nz;
   avgmedslice(medf[0],medf[1],medf[2],medf[3],nx,nz,mi->media_boundary);

   /* rotate pointers */
   tmp_ptr = medf[1];
   medf[1] = medf[2];
   medf[2] = medf[3];
   medf[3] = tmp_ptr;

   if((float)(100.0*(iy-1))/(float)(ny-3) >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy-1,ny-3);
      ptol = ptol + print_tol;
      }
   }

/* replicate first and last plane */

medf[0] = mfield;
copy_avgmedslice(mfield+N_MED_VARS*nx*nz,medf[0],nx,nz);

medf[0] = mfield + (ny-1)*N_MED_VARS*nx*nz;
copy_avgmedslice(mfield+(ny-2)*N_MED_VARS*nx*nz,medf[0],nx,nz);

free(medbuf);
return(1);
}

void init_mediaP3(struct media_input *mi,float *mfield,struct runparamsP3 *rpars,struct qvalues *qv,int swapb)
{
int nx, ny, nz, globnz, fs;
float h;

struct nodeinfo *ni;
struct medstencil *ms;
struct medprof *mp;
float *medf;
int i, iy, iymin, iymax, print_tol, ptol;
int ix, iz;
int fdp, fds, fdd;
int nbpad = NBND_PAD;
off_t boff;

ni = &(rpars->ni);

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;
globnz = ni->globnz;

fs = rpars->freesurf;
h = rpars->h;

print_tol = 25;
ptol = print_tol;

if(mi->model_style == 0)  /* generate model from stencil */
   {
   mi->medprofs = (struct medprof *) check_malloc (NPROFS*sizeof(struct medprof));
   mi->medsten = (struct medstencil *) check_malloc (nx*ny*sizeof(struct medstencil));

   ms = mi->medsten;
   mp = mi->medprofs;

   for(i=0;i<NPROFS;i++)
      mp[i].medptr = (float *) check_malloc (N_MEDPROF_VARS*globnz*sizeof(float));

   genmod3dP3(&mi->nprof,mi->modfile,rpars,ms,mp,mi->media_boundary);
   }
else if(mi->model_style == 1) /* read in binary velocity model from disk */
   {
   fdp = opfile_ro(mi->pmodfile);
   fds = opfile_ro(mi->smodfile);
   fdd = opfile_ro(mi->dmodfile);

   if(qv->qpfrac <= 0.0 || qv->qsfrac <= 0.0)
      read_qvals(mi->qmodfile,qv);

   boff = (off_t)(ni->ny1)*(off_t)(ni->globnx)*(off_t)(ni->globnz*sizeof(float));
   lseek(fdp,boff,SEEK_SET);
   lseek(fds,boff,SEEK_SET);
   lseek(fdd,boff,SEEK_SET);
   }

for(iy=0;iy<ny;iy++)
   {
   medf = mfield + iy*N_MED_VARS*nx*nz;

   if(mi->model_style == 0)
      getmedsliceP3(ms+iy*nx,mp,mi->nprof,medf,rpars,iy,mi->dampwidth,&mi->qbndmax,&mi->qbndmin);
   else if(mi->model_style == 1)
      reedmedsliceP3(fdp,fds,fdd,medf,fs,qv,mi->dampwidth,&mi->qbndmax,&mi->qbndmin,swapb,(ni->ny1)+iy,ni);

   if((float)(100.0*(iy+1))/(float)(ny) >= (float)(ptol))
      {
      fprintf(stderr,"\t %3d percent done (%5d of %5d)\n",ptol,iy+1,ny);
      ptol = ptol + print_tol;
      }
   }

if(mi->model_style == 1)
   {
   close(fdp);
   close(fds);
   close(fdd);
   }

/*
   Add padding to iy=0,ny-1 ends of model to ensure stability of absorbing
   boundaries.  There should be no media variations in the direction
   perpendicular to the absorbing boundary.  See also, xzpad_medslice()
   for padding of x and z boundaries.
*/

pad_medslice(mfield,ni,fs);
}

/*
   reedmedsliceP3() reads the media parameters into an xz slice of the model.

   These values are loaded into the array medf as follows:

      medf           -> lambda + 2*mu
      medf +   nx*nz -> lambda
      medf + 4*nx*nz -> 1/mu
      medf + 7*nx*nz -> rho
      medf + 9*nx*nz -> qp  (later turned into pk coefs.)
      medf +10*nx*nz -> qs  (later turned into sk coefs.)
      medf +11 nx*nz -> lambda (raw, never modified by avg. or atten. )
      medf +12*nx*nz -> mu     (raw, never modified by avg. or atten. )

   After reedmedslice() returns to main(), a call to avgmedslice() computes the
   smoothed values of rigidity (mu) and bouyancy (1/rho) which are then
   stored in the array medf and are ready to be used in the differencing
   operations.

*/

void reedmedsliceP3(int fdp,int fds,int fdd,float *medf,int fs,struct qvalues *qv,int dw,float *qbmax,float *qbmin,int swapb,int glob_iy,struct nodeinfo *ni)
{
float *a, *b;
float *lam2mu, *mu, *lam, *rho, *qp, *qs;
float *lamraw, *muraw;
int i, ip, k, ix, iz, nn, shft;
int powr, xpowr, ypowr, zpowr;
float rfac, qbndf, qfac0, qp0, qs0, qfac;
off_t boff;

float *fbuf;
int jp, kp;

int glob_ix, glob_iz;

float one = 1.0;
float two = 2.0;

rfac = one/(dw - NBND_PAD - 1);
qbndf = (*qbmin)/(*qbmax);
qfac0 = exp(rfac*log(qbndf));

nn = (ni->loc_nx)*(ni->loc_nz);

a = medf + 5*nn;  /* temporary, not used until avgmedslice */
b = medf + 6*nn;  /* temporary, not used until avgmedslice */

lam2mu = medf;
lam    = medf + nn;
mu     = medf + 4*nn;
rho    = medf + 7*nn;
qp     = medf + 9*nn;
qs     = medf + 10*nn;

/*
   Store raw values of lambda and mu into medf array (11*nx*nz and
   12*nx*nz offsets).  These values will never be modified by averaging
   or attenuation (phase velocity adjustments).
*/

lamraw = medf + 11*nn;
muraw  = medf + 12*nn;

boff = (off_t)(ni->globnx)*(off_t)((ni->nz1)*sizeof(float));
lseek(fdp,boff,SEEK_CUR);
lseek(fds,boff,SEEK_CUR);
lseek(fdd,boff,SEEK_CUR);

/*  old way with a bunch of lseeks, pretty slow?
for(iz=0;iz<(ni->loc_nz);iz++)
   {
   boff = (off_t)((ni->nx1)*sizeof(float));
   lseek(fdp,boff,SEEK_CUR);
   lseek(fds,boff,SEEK_CUR);
   lseek(fdd,boff,SEEK_CUR);

   reed(fdp,a + iz*(ni->loc_nx),(ni->loc_nx)*sizeof(float));
   reed(fds,b + iz*(ni->loc_nx),(ni->loc_nx)*sizeof(float));
   reed(fdd,rho + iz*(ni->loc_nx),(ni->loc_nx)*sizeof(float));

   boff = (off_t)((ni->globnx - ni->nx2)*sizeof(float));
   lseek(fdp,boff,SEEK_CUR);
   lseek(fds,boff,SEEK_CUR);
   lseek(fdd,boff,SEEK_CUR);
   }
*/

/*  new way: read large chunk and then copy out subset into arrays */
fbuf = (float *) check_malloc (ni->globnx*ni->loc_nz*sizeof(float));

reed(fdp,fbuf,ni->globnx*ni->loc_nz*sizeof(float));
for(iz=0;iz<(ni->loc_nz);iz++)
   {
   for(ix=0;ix<(ni->loc_nx);ix++)
      {
      jp = ix + iz*ni->loc_nx;
      kp = (ni->nx1+ix) + iz*ni->globnx;
      a[jp] = fbuf[kp];
      }
   }

reed(fds,fbuf,ni->globnx*ni->loc_nz*sizeof(float));
for(iz=0;iz<(ni->loc_nz);iz++)
   {
   for(ix=0;ix<(ni->loc_nx);ix++)
      {
      jp = ix + iz*ni->loc_nx;
      kp = (ni->nx1+ix) + iz*ni->globnx;
      b[jp] = fbuf[kp];
      }
   }

reed(fdd,fbuf,ni->globnx*ni->loc_nz*sizeof(float));
for(iz=0;iz<(ni->loc_nz);iz++)
   {
   for(ix=0;ix<(ni->loc_nx);ix++)
      {
      jp = ix + iz*ni->loc_nx;
      kp = (ni->nx1+ix) + iz*ni->globnx;
      rho[jp] = fbuf[kp];
      }
   }

free(fbuf);
/* end of new way */

/* seek to next global iy plane */
boff = (off_t)(ni->globnx)*(off_t)((ni->globnz - ni->nz2)*sizeof(float));
lseek(fdp,boff,SEEK_CUR);
lseek(fds,boff,SEEK_CUR);
lseek(fdd,boff,SEEK_CUR);

if(swapb == 1) /* need to swap bytes */
   {
   swap_in_place(nn,((char *)(a)));
   swap_in_place(nn,((char *)(b)));
   swap_in_place(nn,((char *)(rho)));
   }

for(ix=0;ix<(ni->loc_nx);ix++)
   {
   lam2mu[ix] = a[ix]*a[ix]*rho[ix];
   mu[ix] = b[ix]*b[ix]*rho[ix];

   lam[ix] = lam2mu[ix] - two*mu[ix];
   mu[ix] = one/mu[ix];
   lamraw[ix] = lam[ix];
   muraw[ix] = one/mu[ix];

   get_qvals(&qp[ix],&qs[ix],&a[ix],&b[ix],qv);

   if(qp[ix] <= 0.0 || qs[ix] <= 0.0)
      {
      fprintf(stderr,"Problem with zero/negative Q, exiting...\n");
      exit(-99);
      }
   }

shft = 0;
if(fs)      /* shift model down one grid point for free surface */
   shft = 1;

for(iz=(ni->loc_nz)-1;iz>=1;iz--)
   {
   glob_iz = iz + (ni->nz1);
   for(ix=0;ix<(ni->loc_nx);ix++)
      {
      glob_ix = ix + (ni->nx1);

      i = iz*(ni->loc_nx) + ix;
      ip = (iz-shft)*(ni->loc_nx) + ix;

      lam2mu[i] = a[ip]*a[ip]*rho[ip];
      mu[i] = b[ip]*b[ip]*rho[ip];

      lam[i] = lam2mu[i] - two*mu[i];
      mu[i] = one/mu[i];
      lamraw[i] = lam[i];
      muraw[i] = one/mu[i];
      rho[i] = rho[ip];

      get_qvals(&qp[i],&qs[i],&a[ip],&b[ip],qv);

/* damping region for absorbing boundary */

      xpowr = 0;
      if(glob_ix < dw)
	 xpowr = dw - glob_ix;
      else if(glob_ix >= (ni->globnx)-dw)
	 xpowr = glob_ix - ((ni->globnx) - 1 - dw);

      ypowr = 0;
      if(glob_iy < dw)
	 ypowr = dw - glob_iy;
      else if(glob_iy >= (ni->globny)-dw)
	 ypowr = glob_iy - ((ni->globny) - 1 - dw);

      zpowr = 0;
      if(glob_iz < dw && fs == 0)
	 zpowr = dw - glob_iz;
      else if(glob_iz >= (ni->globnz)-dw)
	 zpowr = glob_iz - ((ni->globnz) - 1 - dw);

      if(xpowr || ypowr || zpowr)
         {
         powr = xpowr + ypowr + zpowr - 1;
	 qfac = one;
         while(powr--)
            qfac = qfac*qfac0;

         qp0 = qp[i];
         if(qp0 > 2*(*qbmax))
            qp0 = 2*(*qbmax);

         qs0 = qs[i];
         if(qs0 > (*qbmax))
            qs0 = (*qbmax);

         qp[i] = qfac*qp0;
         qs[i] = qfac*qs0;

         if(qp[i] < 5.0)
            qp[i] = 5.0;
         if(qs[i] < 5.0)
            qs[i] = 5.0;

         qp[i] = -qp[i];
         qs[i] = -qs[i];
         }
      }
   }
}

void pad_medslice(float *mfield,struct nodeinfo *ni,int fs)
{
float *mdf1;
float *l2m, *lam, *mu, *boy, *qp, *qs;
float *lraw, *mraw;
int ix, iy, iz, ip, off1;
int nx, ny, nz;

int npad = NBND_PAD;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

for(iy=0;iy<ny;iy++)
   {
   mdf1 = mfield + N_MED_VARS*iy*nx*nz;

   l2m = mdf1;
   lam = mdf1 +   nx*nz;
   mu  = mdf1 + 4*nx*nz;
   boy = mdf1 + 7*nx*nz;
   qp  = mdf1 + 9*nx*nz;
   qs  = mdf1 + 10*nx*nz;
   lraw = mdf1 + 11*nx*nz;
   mraw = mdf1 + 12*nx*nz;

   for(iz=0;iz<nz;iz++)
      {
      if(ni->minusId_x < 0)
         {
         for(ix=0;ix<npad;ix++)
            {
            l2m[iz*nx + ix] = l2m[iz*nx + npad];
            lam[iz*nx + ix] = lam[iz*nx + npad];
            mu[iz*nx + ix]  = mu[iz*nx + npad];
            boy[iz*nx + ix] = boy[iz*nx + npad];
            qp[iz*nx + ix]  = qp[iz*nx + npad];
            qs[iz*nx + ix]  = qs[iz*nx + npad];
            lraw[iz*nx + ix]  = lraw[iz*nx + npad];
            mraw[iz*nx + ix]  = mraw[iz*nx + npad];
            }
	 }
      if(ni->plusId_x < 0)
         {
         for(ix=nx-npad;ix<nx;ix++)
            {
            l2m[iz*nx + ix] = l2m[iz*nx + (nx-npad-1)];
            lam[iz*nx + ix] = lam[iz*nx + (nx-npad-1)];
            mu[iz*nx + ix]  = mu[iz*nx + (nx-npad-1)];
            boy[iz*nx + ix] = boy[iz*nx + (nx-npad-1)];
            qp[iz*nx + ix]  = qp[iz*nx + (nx-npad-1)];
            qs[iz*nx + ix]  = qs[iz*nx + (nx-npad-1)];
            lraw[iz*nx + ix]  = lraw[iz*nx + (nx-npad-1)];
            mraw[iz*nx + ix]  = mraw[iz*nx + (nx-npad-1)];
            }
         }
      }

   if(!fs && ni->minusId_z < 0)
      {
      for(iz=0;iz<npad;iz++)
         {
         for(ix=0;ix<nx;ix++)
            {
            l2m[iz*nx + ix] = l2m[npad*nx + ix];
            lam[iz*nx + ix] = lam[npad*nx + ix];
            mu[iz*nx + ix]  = mu[npad*nx + ix];
            boy[iz*nx + ix] = boy[npad*nx + ix];
            qp[iz*nx + ix]  = qp[npad*nx + ix];
            qs[iz*nx + ix]  = qs[npad*nx + ix];
            lraw[iz*nx + ix]  = lraw[npad*nx + ix];
            mraw[iz*nx + ix]  = mraw[npad*nx + ix];
            }
         }
      }
   if(ni->plusId_z < 0)
      {
      for(iz=nz-npad;iz<nz;iz++)
         {
         for(ix=0;ix<nx;ix++)
            {
            l2m[iz*nx + ix] = l2m[(nz-npad-1)*nx + ix];
            lam[iz*nx + ix] = lam[(nz-npad-1)*nx + ix];
            mu[iz*nx + ix]  = mu[(nz-npad-1)*nx + ix];
            boy[iz*nx + ix] = boy[(nz-npad-1)*nx + ix];
            qp[iz*nx + ix]  = qp[(nz-npad-1)*nx + ix];
            qs[iz*nx + ix]  = qs[(nz-npad-1)*nx + ix];
            lraw[iz*nx + ix]  = lraw[(nz-npad-1)*nx + ix];
            mraw[iz*nx + ix]  = mraw[(nz-npad-1)*nx + ix];
            }
         }
      }
   }

if(ni->minusId_y < 0)
   {
   mdf1 = mfield + N_MED_VARS*npad*nx*nz;

   l2m = mdf1;
   lam = mdf1 +   nx*nz;
   mu  = mdf1 + 4*nx*nz;
   boy = mdf1 + 7*nx*nz;
   qp  = mdf1 + 9*nx*nz;
   qs  = mdf1 + 10*nx*nz;
   lraw = mdf1 + 11*nx*nz;
   mraw = mdf1 + 12*nx*nz;

   for(iy=0;iy<npad;iy++)
      {
      off1 = N_MED_VARS*(iy+1)*nx*nz;

      for(ip=0;ip<nx*nz;ip++)
         {
         l2m[ip-off1] = l2m[ip];
         lam[ip-off1] = lam[ip];
         mu[ip-off1]  = mu[ip];
         boy[ip-off1] = boy[ip];
         qp[ip-off1]  = qp[ip];
         qs[ip-off1]  = qs[ip];
         lraw[ip-off1]  = lraw[ip];
         mraw[ip-off1]  = mraw[ip];
         }
      }
   }

if(ni->plusId_y < 0)
   {
   mdf1 = mfield + N_MED_VARS*(ny-npad-1)*nx*nz;

   l2m = mdf1;
   lam = mdf1 +   nx*nz;
   mu  = mdf1 + 4*nx*nz;
   boy = mdf1 + 7*nx*nz;
   qp  = mdf1 + 9*nx*nz;
   qs  = mdf1 + 10*nx*nz;
   lraw = mdf1 + 11*nx*nz;
   mraw = mdf1 + 12*nx*nz;

   for(iy=0;iy<npad;iy++)
      {
      off1 = N_MED_VARS*(iy+1)*nx*nz;

      for(ip=0;ip<nx*nz;ip++)
         {
         l2m[ip+off1] = l2m[ip];
         lam[ip+off1] = lam[ip];
         mu[ip+off1]  = mu[ip];
         boy[ip+off1] = boy[ip];
         qp[ip+off1]  = qp[ip];
         qs[ip+off1]  = qs[ip];
         lraw[ip+off1]  = lraw[ip];
         mraw[ip+off1]  = mraw[ip];
         }
      }
   }
}

void genmod3dP3(int *nprof,char *mfile,struct runparamsP3 *rpars,struct medstencil *medsten,struct medprof *medprofs,int medbnd)
{
FILE *fpr, *fopfile();
char *words[128], *line, *getlineRWG();
int nword;
char buf[1024];

int ip = 0;
int linenum = 0;
int ungetflag = 0;

struct nodeinfo *ni;
int nx, ny, nz, globnz, fs;
float h;

ni = &(rpars->ni);
fs = rpars->freesurf;
h = rpars->h;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;
globnz = ni->globnz;

do_init(nx*ny,medsten,globnz);

fpr = fopfile(mfile,"r");

while((line = getlineRWG(fpr,&linenum,&ungetflag,buf)) != NULL)
   {
   if(line[0] == '\0') continue;     /* blank line */
   if(line[0] == '#') continue;      /* comment */
   nword = getwords(line,words);
   if(nword == 0) continue;          /* nothing on line */

   if(strcmp("DEF",words[0]) == 0)
      {
      do_prof(fpr,&medprofs[ip],&words[1],globnz,&h,fs,&linenum,&ungetflag,buf,medbnd);
      ip++;
      continue;
      }
   }

fclose(fpr);

if(fs) /* Adjust media stencil to add extra row at model top for free-surface */
   adj_medsten_fs(nx,ny,globnz,medsten,fs);

/*
   Set max. depths for basin profiles and check to ensure that input
   profiles are deep enough.
*/

set_max_depths(medsten,medprofs,ip,nx,ny);

*nprof = ip;
}

/*
   getmedsliceP3() loads the media parameters into an xz slice of the model.
   First, the 'HOST' (HST) media profiles are loaded into the array medf,
   then the appropriate basin profile values are loaded over these.  The media
   profiles stored in the array medptr in the structure mdp are as follows:

      mdp[iprof].medptr        -> lambda + 2*mu
      mdp[iprof].medptr +   nz -> lambda
      mdp[iprof].medptr + 2*nz -> 1/mu
      mdp[iprof].medptr + 3*nz -> rho
      mdp[iprof].medptr + 4*nz -> qp factor
      mdp[iprof].medptr + 5*nz -> qs factor

   These values are loaded into the array medf as follows:

      medf           -> lambda + 2*mu
      medf +   nx*nz -> lambda
      medf + 4*nx*nz -> 1/mu
      medf + 7*nx*nz -> rho
      medf + 9*nx*nz -> qp  (later turned into pk coefs.)
      medf +10*nx*nz -> qs  (later turned into sk coefs.)
      medf +11 nx*nz -> lambda (raw, never modified by avg. or atten. )
      medf +12*nx*nz -> mu     (raw, never modified by avg. or atten. )

   After getmedsliceP3() returns to main(), a call to avgmedslice() computes the
   smoothed values of rigidity (mu) and bouyancy (1/rho) which are then
   stored in the array medf and are ready to be used in the differencing
   operations.

*/

void getmedsliceP3(struct medstencil *ms,struct medprof *mdp,int nprof,float *medf,struct runparamsP3 *rpars,int iy,int dw,float *qbmax,float *qbmin)
{
float *mptr, *m2ptr;
float *mptr1, *mptr2, *mptr3, *mptr4;
int iprof, ix, izst, izend;

float *lam2mu, *mu, *lam;
float *qp, *qs, qp0, qs0, qfac;
float rfac, qbndf, qfac0;
int i, ip, powr, xpowr, ypowr, zpowr;
int iz, iend;

float one = 1.0;
float two = 2.0;

struct nodeinfo *ni;
int nx, ny, nz, nz1, globnz, fs;
float h;

ni = &(rpars->ni);
fs = rpars->freesurf;
h = rpars->h;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;
nz1 = ni->nz1;
globnz = ni->globnz;

	     /*  initialize host structure  */

iprof = 0;
while(strcmp("HST",mdp[iprof].prof_name) != 0)
   {
   iprof++;
   if(iprof == nprof)
      {
      fprintf(stderr,"\n***********  ERROR  ***********\n");
      fprintf(stderr,"Host media, 'HST', does not have an DEF ");
      fprintf(stderr,"entry in model file.\n");
      exit(-1);
      }
   }

for(ix=0;ix<nx;ix++)
   {
   mptr = medf + ix;
   m2ptr = (mdp[iprof].medptr) + nz1;
   store(nz,m2ptr,         1,mptr,         nx);  /* lambda + 2*mu */
   store(nz,m2ptr+  globnz,1,mptr+   nx*nz,nx);  /* lambda */
   store(nz,m2ptr+2*globnz,1,mptr+ 4*nx*nz,nx);  /* 1/mu */
   store(nz,m2ptr+3*globnz,1,mptr+ 7*nx*nz,nx);  /* rho */
   store(nz,m2ptr+4*globnz,1,mptr+ 9*nx*nz,nx);  /* qp */
   store(nz,m2ptr+5*globnz,1,mptr+10*nx*nz,nx);  /* qs */
   }

/* Apply damping function near boundaries */

qp = medf + 9*nx*nz;
qs = medf + 10*nx*nz;

rfac = one/(dw - NBND_PAD - 1);
qbndf = (*qbmin)/(*qbmax);
qfac0 = exp(rfac*log(qbndf));

ypowr = 0;
if(iy < dw && ni->minusId_y == -1)
   ypowr = dw - iy;
else if(iy >= ny-dw && ni->plusId_y == -1)
   ypowr = iy - (ny - 1 - dw);

for(iz=0;iz<nz;iz++)
   {
   zpowr = 0;
   if(iz < dw && fs == 0 && ni->minusId_z == -1)
      zpowr = dw - iz;
   else if(iz >= nz-dw && ni->plusId_z == -1)
      zpowr = iz - (nz - 1 - dw);

   for(ix=0;ix<nx;ix++)
      {
      ip = ix + iz*nx;

      xpowr = 0;
      if(ix < dw && ni->minusId_z == -1)
	 xpowr = dw - ix;
      else if(ix >= nx-dw && ni->plusId_z == -1)
	 xpowr = ix - (nx - 1 - dw);

      if(xpowr || ypowr || zpowr)
         {
         powr = xpowr + ypowr + zpowr - 1;
         qfac = one;
         while(powr--)
            qfac = qfac*qfac0;

	 qp0 = qp[ip];
	 if(qp0 > 2*(*qbmax))
	    qp0 = 2*(*qbmax);

	 qs0 = qs[ip];
	 if(qs0 > (*qbmax))
	    qs0 = (*qbmax);

	 qp[ip] = qfac*qp0;
	 qs[ip] = qfac*qs0;

	 if(qp[ip] < 5.0)
	    qp[ip] = 5.0;
	 if(qs[ip] < 5.0)
	    qs[ip] = 5.0;

         qp[ip] = -qp[ip];
         qs[ip] = -qs[ip];
         }
      }
   }

/*
   Store raw values of lambda and mu into medf array (11*nx*nz and
   12*nx*nz offsets).  These values will never be modified by averaging
   or attenuation (phase velocity adjustments).
*/

mptr1 = medf +    nx*nz;     /* lambda */
mptr2 = medf + 11*nx*nz;     /* lambda raw */
mptr3 = medf +  4*nx*nz;     /* 1/mu */
mptr4 = medf + 12*nx*nz;     /* mu raw */

for(iz=0;iz<nz;iz++)
   {
   for(ix=0;ix<nx;ix++)
      {
      ip = ix+iz*nx;
      mptr2[ip] = mptr1[ip];
      mptr4[ip] = one/mptr3[ip];
      }
   }
}
