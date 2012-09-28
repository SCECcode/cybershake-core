#include "mpi.h"
#include "include.h"
#include "function.h"
#include "func_mpi.h"

#define         MAXFILES       10000
#define         MAXLINE       256
#define         MAXMEM       1500   /* default max RAM in Mbytes */
#define MAX_NLEN 100

struct fdcf
   {
   float c4[4];
   float c2[2];
   };

struct gridparam
   {
   struct fdcf *cfx0;
   struct fdcf *cfx1;
   struct fdcf *cfy0;
   struct fdcf *cfy1;
   struct fdcf *cfz0;
   struct fdcf *cfz1;
   float *hx;
   float *hy;
   float *hz;
   float *xp;
   float *yp;
   float *zp;
   float hmin;
   float hmax;
   float dt;
   float dthx0;
   float dthxn;
   float dthy0;
   float dthyn;
   float dthz0;
   float dthzn;
   int order;
   int nx;
   int ny;
   int nz;
   };

void set_gridparams(char *,struct gridparam *, char *,int);
void setcoef(struct fdcf *,struct fdcf *,float *,int);
void gelim_double(double *,int,double *);
void getlens(FILE *,char *,float *,float *,float *,float *,int *);
void check_vmax(int,int,int,float *,float *,float *,struct gridparam *,int *,int *,int *,float *,float *,float *);
void check_vmin(int fdp,int fds,int fdd,float *vp,float *vs,float *den,struct gridparam *gp);
void set_minmax(float *,float *,float *,int);
void set_vpvs(float *,float *,int);

char filebuf[MAXLINE*MAXFILES];
char stringbuf[MAXLINE];

int main(int ac, char **av)
{
/* MPI-RWG: MPI stuff and distributed computation variables */
int myid, nproc, pnlen;
char procname[128];

struct gridparam gp;
char gridfile[256], gridout[256];
FILE *fpr, *fopfile();
float *fbuf, *vp, *vs, *dn;
int ixmax, iymax, izmax;
float vmax, hmin, dtmax;
int i, j, k, l, my_k, glob_k, nfiles;
int ip, kp, lp, is, id;
off_t off, head1_off, cur_off;
int fdr, fdp, fds, fdd;
char *infilebuf, *infile;
char filelist[512], logfile[128], logname[128];
char outfile[512], pmodfile[256], smodfile[256], dmodfile[256];
char logdir[128], outdir[256];
char str[512];

int nx, ny, nz, nx_in, ny_in, nz_in;
int iy, iy0, local_ny, glob_iy, iy_inc;

off_t bsize;
off_t maxmem = -1;

int maxfiles = MAXFILES;

int check_poisson_ratio = 1;

float pmin = -1.0;
float smin = -1.0;
float dmin = -1.0;

int smin_true, check_pmin;
int always_check_pmin = 1;

float pmax = 1.0e+10;
float smax = 1.0e+10;
float dmax = 1.0e+10;

float vpbackg = 1.0e+15;

float fdversion = 1.0;

/*  MPI-RWG: Start-up MPI */

mpi_init(&ac,&av,&nproc,&myid,procname,&pnlen);

sprintf(gridout,"gridout");
sprintf(outdir,".");
sprintf(logdir,".");
sprintf(logname,"merge_cvm");

setpar(ac, av);
mstpar("gridfile","s",gridfile);
getpar("gridout","s",gridout);
mstpar("filelist","s",filelist);

getpar("outdir","s",outdir);
mstpar("pmodfile","s",pmodfile);
mstpar("smodfile","s",smodfile);
mstpar("dmodfile","s",dmodfile);

getpar("logdir","s",logdir);
getpar("logname","s",logname);
getpar("maxfiles","d",&maxfiles);
getpar("maxmem","d",&maxmem);

getpar("check_poisson_ratio","d",&check_poisson_ratio);

getpar("pmin","f",&pmin);
getpar("smin","f",&smin);
getpar("dmin","f",&dmin);
getpar("pmax","f",&pmax);
getpar("smax","f",&smax);
getpar("dmax","f",&dmax);

getpar("vpbackg","f",&vpbackg);
getpar("always_check_pmin","d",&always_check_pmin);

getpar("fdversion","f",&fdversion);
endpar();

makedir(outdir);
makedir(logdir);

sprintf(logfile,"%s/%s-%.4d.log",logdir,logname,myid);
freopen(logfile,"w",stderr);

fprintf(stderr,"***** This is merge_sgt\n");
fprintf(stderr,"      Running on node= %s rank= %d\n\n",procname,myid);

fprintf(stderr,"**** Files to process are (filelist= %s):\n",filelist);
fflush(stderr);

if(myid != 0)
   gridout[0] = '\0';

set_gridparams(gridfile,&gp,gridout,0);  /* free surface = 0 */

nx = gp.nx;
ny = gp.ny;
nz = gp.nz;

head1_off = 0;
if(fdversion >= 2.0)
   head1_off = (off_t)(3*sizeof(int) + (nx+ny+nz)*sizeof(float));

infilebuf = (char *) check_malloc (MAXLINE*maxfiles*sizeof(char));

fpr = fopfile(filelist,"r");

k = 0;
my_k = 0;
glob_k = 0;
while(fgets(str,512,fpr) != NULL)
   {
   if((glob_k - k*nproc) == myid)
      {
      infile = infilebuf + my_k*MAXLINE;
      sscanf(str,"%s",infile);

      fprintf(stderr,"     %s\n",infile);
      fflush(stderr);

      my_k++;
      }

   glob_k++;

   if((glob_k - k*nproc) == nproc)
      k++;
   }
fclose(fpr);
nfiles = my_k;

fprintf(stderr,"\n");
fflush(stderr);

if(maxmem < 0)
   maxmem = MAXMEM;

fbuf = NULL;
for(k=0;k<nfiles;k++)
   {
   infile = infilebuf + k*MAXLINE;

   fprintf(stderr,"processing file= %s\n",infile);
   fflush(stderr);

   fdr = opfile_ro(infile);

   reed(fdr,&nx_in,sizeof(int));
   reed(fdr,&ny_in,sizeof(int));
   reed(fdr,&nz_in,sizeof(int));
   reed(fdr,&iy0,sizeof(int));
   reed(fdr,&local_ny,sizeof(int));

   if(nx_in != nx || ny_in != ny || nz_in != nz)
      {
      fprintf(stderr,"*** Model parameters not equal, exiting ...\n");
      fprintf(stderr,"                     nx_in= %6d nx= %6d\n",nx_in,nx);
      fprintf(stderr,"                     ny_in= %6d ny= %6d\n",ny_in,ny);
      fprintf(stderr,"                     nz_in= %6d nz= %6d\n",nz_in,nz);
      exit(-1);
      }

   iy_inc = local_ny;
   bsize = 3*iy_inc*nx*nz*sizeof(float);
   while(bsize > maxmem*1000000 && iy_inc > 1)
      {
      iy_inc--;
      bsize = 3*iy_inc*nx*nz*sizeof(float);
      }

   fbuf = (float *) check_realloc (fbuf,bsize);
   for(i=0;i<3*iy_inc*nx*nz;i++)
      fbuf[i] = -1.0;

   if(k == 0)   /* open and initialize output files */
      {
      if(myid == 0)  /* root node creates and initializes */
         {
	 sprintf(outfile,"%s/%s",outdir,pmodfile);
         fdp = croptrfile_sync(outfile);

	 sprintf(outfile,"%s/%s",outdir,smodfile);
         fds = croptrfile_sync(outfile);

	 sprintf(outfile,"%s/%s",outdir,dmodfile);
         fdd = croptrfile_sync(outfile);

         if(fdversion >= 2.0)
            {
            rite(fdp,&gp.nx,sizeof(int));
            rite(fdp,gp.xp,nx*sizeof(float));
            rite(fdp,&gp.nz,sizeof(int));
            rite(fdp,gp.zp,nz*sizeof(float));
            rite(fdp,&gp.ny,sizeof(int));
            rite(fdp,gp.yp,ny*sizeof(float));

            rite(fds,&gp.nx,sizeof(int));
            rite(fds,gp.xp,nx*sizeof(float));
            rite(fds,&gp.nz,sizeof(int));
            rite(fds,gp.zp,nz*sizeof(float));
            rite(fds,&gp.ny,sizeof(int));
            rite(fds,gp.yp,ny*sizeof(float));

            rite(fdd,&gp.nx,sizeof(int));
            rite(fdd,gp.xp,nx*sizeof(float));
            rite(fdd,&gp.nz,sizeof(int));
            rite(fdd,gp.zp,nz*sizeof(float));
            rite(fdd,&gp.ny,sizeof(int));
            rite(fdd,gp.yp,ny*sizeof(float));
            }

	 for(iy=0;iy<ny;iy=iy+3*iy_inc)
	    {
	    if((iy + 3*iy_inc) >= ny)
               bsize = (ny-iy)*nx*nz*sizeof(float);

            rite(fdp,fbuf,bsize);
            rite(fds,fbuf,bsize);
            rite(fdd,fbuf,bsize);
	    }

         close(fdp);
         close(fds);
         close(fdd);

         MPI_Barrier(MPI_COMM_WORLD);
         }
      else  /* other nodes wait, then just open */
         MPI_Barrier(MPI_COMM_WORLD);

      sprintf(outfile,"%s/%s",outdir,pmodfile);
      fdp = opfile(outfile);

      sprintf(outfile,"%s/%s",outdir,smodfile);
      fds = opfile(outfile);

      sprintf(outfile,"%s/%s",outdir,dmodfile);
      fdd = opfile(outfile);

      lseek(fdp,head1_off,SEEK_SET);
      lseek(fds,head1_off,SEEK_SET);
      lseek(fdd,head1_off,SEEK_SET);

      cur_off = head1_off;
      }

   bsize = nx*nz*sizeof(float);
   for(iy=0;iy<local_ny;iy=iy+iy_inc)
      {
      if(iy+iy_inc > local_ny)
         iy_inc = local_ny - iy;

      reed(fdr,fbuf,3*iy_inc*bsize);

      for(j=0;j<iy_inc;j++)
         {
	 glob_iy = iy + j + iy0;

	 vp = fbuf + 3*j*nx*nz;
	 vs = fbuf + (3*j + 1)*nx*nz;
	 dn = fbuf + (3*j + 2)*nx*nz;

/* check for negative poissons ratio, scale Vp until positive VERY KLUGY!!! */
	 if(check_poisson_ratio)
	    set_vpvs(vp,vs,nx*nz);

	 set_minmax(vp,&pmin,&pmax,nx*nz);
	 set_minmax(vs,&smin,&smax,nx*nz);
	 set_minmax(dn,&dmin,&dmax,nx*nz);

         off = head1_off + (off_t)(glob_iy*bsize) - cur_off;

         lseek(fdp,off,SEEK_CUR);
	 rite(fdp,vp,bsize);

         lseek(fds,off,SEEK_CUR);
	 rite(fds,vs,bsize);

         lseek(fdd,off,SEEK_CUR);
	 rite(fdd,dn,bsize);

         cur_off = head1_off + (glob_iy+1)*bsize;
	 }
      }

   close(fdr);
   }

close(fdp);
close(fds);
close(fdd);

fprintf(stderr,"\n");
fprintf(stderr,"**** File processing done.  Check for stability criteria.\n\n");
fflush(stderr);

if(myid == 0)
   {
   sprintf(outfile,"%s/%s",outdir,pmodfile);
   fdp = opfile(outfile);

   sprintf(outfile,"%s/%s",outdir,smodfile);
   fds = opfile(outfile);

   sprintf(outfile,"%s/%s",outdir,dmodfile);
   fdd = opfile(outfile);

   lseek(fdp,head1_off,SEEK_SET);
   lseek(fds,head1_off,SEEK_SET);
   lseek(fdd,head1_off,SEEK_SET);

   check_vmax(fdp,fds,fdd,vp,vs,dn,&gp,&ixmax,&iymax,&izmax,&vmax,&hmin,&dtmax);
   check_vmin(fdp,fds,fdd,vp,vs,dn,&gp);

   close(fdp);
   close(fds);
   close(fdd);

   MPI_Barrier(MPI_COMM_WORLD);
   }
else
   MPI_Barrier(MPI_COMM_WORLD);

MPI_Bcast(&ixmax,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&iymax,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&izmax,1,MPI_INT,0,MPI_COMM_WORLD);
MPI_Bcast(&vmax,1,MPI_FLOAT,0,MPI_COMM_WORLD);
MPI_Bcast(&hmin,1,MPI_FLOAT,0,MPI_COMM_WORLD);
MPI_Bcast(&dtmax,1,MPI_FLOAT,0,MPI_COMM_WORLD);

fprintf(stderr,"  At ix= %6d iy= %6d iz= %6d:\n",ixmax,iymax,izmax);
fprintf(stderr,"     vmax= %.6f km/s\n",vmax);
fprintf(stderr,"     hmin= %.6f km\n\n",hmin);
fprintf(stderr,"  -for 4th order, dt must be less than %.6f\n",6.0*dtmax/7.0);
fprintf(stderr,"  -for 2nd order, dt must be less than %.6f\n",dtmax);
fflush(stderr);

close(fdp);
close(fds);
close(fdd);
mpi_exit(0);
}

void set_gridparams(char *gfile,struct gridparam *gp, char *gout,int fs)
{
float x0[MAX_NLEN], x1[MAX_NLEN], dxlen[MAX_NLEN];
float y0[MAX_NLEN], y1[MAX_NLEN], dylen[MAX_NLEN];
float z0[MAX_NLEN], z1[MAX_NLEN], dzlen[MAX_NLEN];
float xlen, ylen, zlen;
float xx, yy, zz;
double dd;
int nxlen, nylen, nzlen;
int nx, ny, nz;
int i, j;

FILE *fpr, *fpw, *fopfile();
char str[512];

float hmin = 1.0e+15;
float hmax = -1.0e+15;
float szero = 0.0;

fpr = fopfile(gfile,"r");

while(fgets(str,512,fpr) != NULL)
   {
   if(strncmp(str,"xlen=",5) == 0)
      {
      sscanf(&str[5],"%f",&xlen);
      getlens(fpr,str,&xlen,x0,x1,dxlen,&nxlen);
      }
   else if(strncmp(str,"ylen=",5) == 0)
      {
      sscanf(&str[5],"%f",&ylen);
      getlens(fpr,str,&ylen,y0,y1,dylen,&nylen);
      }
   else if(strncmp(str,"zlen=",5) == 0)
      {
      sscanf(&str[5],"%f",&zlen);
      getlens(fpr,str,&zlen,z0,z1,dzlen,&nzlen);
      }
   }
fclose(fpr);

nx = 0;
for(i=0;i<nxlen;i++)
   nx = nx + (int)((x1[i] - x0[i])/dxlen[i] + 0.5);

gp->hx = (float *) check_malloc (nx*sizeof(float));
gp->xp = (float *) check_malloc (nx*sizeof(float));
gp->cfx0 = (struct fdcf *) check_malloc (nx*sizeof(struct fdcf));
gp->cfx1 = (struct fdcf *) check_malloc (nx*sizeof(struct fdcf));
gp->nx = nx;

j = 0;
dd = x0[0];
for(i=0;i<nx;i++)
   {
   gp->hx[i] = dxlen[j];
   dd = dd + gp->hx[i];
   gp->xp[i] = dd - 0.5*gp->hx[i];

   if(dxlen[j] < hmin)
      hmin = dxlen[j];
   if(dxlen[j] > hmax)
      hmax = dxlen[j];

   if(dd >= x1[j] && j < nxlen-1)
      j++;
   }

xx = gp->xp[0];     /* reset origin */

for(i=0;i<nx;i++)
   gp->xp[i] = gp->xp[i] - xx;

setcoef(gp->cfx0,gp->cfx1,gp->hx,gp->nx);

ny = 0;
for(i=0;i<nylen;i++)
   ny = ny + (int)((y1[i] - y0[i])/dylen[i] + 0.5);

gp->hy = (float *) check_malloc (ny*sizeof(float));
gp->yp = (float *) check_malloc (ny*sizeof(float));
gp->cfy0 = (struct fdcf *) check_malloc (ny*sizeof(struct fdcf));
gp->cfy1 = (struct fdcf *) check_malloc (ny*sizeof(struct fdcf));
gp->ny = ny;

j = 0;
dd = y0[0];
for(i=0;i<ny;i++)
   {
   gp->hy[i] = dylen[j];
   dd = dd + gp->hy[i];
   gp->yp[i] = dd - 0.5*gp->hy[i];

   if(dylen[j] < hmin)
      hmin = dylen[j];
   if(dylen[j] > hmax)
      hmax = dylen[j];

   if(dd >= y1[j] && j < nylen-1)
      j++;
   }

yy = gp->yp[0];     /* reset origin */

for(i=0;i<ny;i++)
   gp->yp[i] = gp->yp[i] - yy;

setcoef(gp->cfy0,gp->cfy1,gp->hy,gp->ny);

nz = 0;
for(i=0;i<nzlen;i++)
   nz = nz + (int)((z1[i] - z0[i])/dzlen[i] + 0.5);

if(fs)
   {
   z0[0] = -dzlen[0];
   nz++;
   }

gp->hz = (float *) check_malloc (nz*sizeof(float));
gp->zp = (float *) check_malloc (nz*sizeof(float));
gp->cfz0 = (struct fdcf *) check_malloc (nz*sizeof(struct fdcf));
gp->cfz1 = (struct fdcf *) check_malloc (nz*sizeof(struct fdcf));
gp->nz = nz;

j = 0;
dd = z0[0];
for(i=0;i<nz;i++)
   {
   gp->hz[i] = dzlen[j];
   dd = dd + gp->hz[i];
   gp->zp[i] = dd - 0.5*gp->hz[i];

   if(dzlen[j] < hmin)
      hmin = dzlen[j];
   if(dzlen[j] > hmax)
      hmax = dzlen[j];

   if(dd >= z1[j] && j < nzlen-1)
      j++;
   }

zz = gp->zp[0];     /* reset origin */
if(fs)
   zz = gp->zp[1];

for(i=0;i<nz;i++)
   gp->zp[i] = gp->zp[i] - zz;

setcoef(gp->cfz0,gp->cfz1,gp->hz,gp->nz);

gp->hmin = hmin;
gp->hmax = hmax;

if(gout[0] != '\0')
   {
   fpw = fopfile(gout,"w");

   fprintf(fpw,"xlen=%.5f\n",xlen);
   fprintf(fpw,"nx=%d\n",nx);
   for(i=0;i<nx;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->xp[i],gp->hx[i]);
      if(gp->hx[i] <= (float)(0.0))
         {
         fprintf(stderr,"**** Problem with zero hx, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fprintf(fpw,"ylen=%.5f\n",ylen);
   fprintf(fpw,"ny=%d\n",ny);
   for(i=0;i<ny;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->yp[i],gp->hy[i]);
      if(gp->hy[i] <= szero)
         {
         fprintf(stderr,"**** Problem with zero hy, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fprintf(fpw,"zlen=%.5f\n",zlen);
   fprintf(fpw,"nz=%d\n",nz);
   for(i=0;i<nz;i++)
      {
      fprintf(fpw,"%6d %13.5e %13.5e\n",i,gp->zp[i],gp->hz[i]);
      if(gp->hz[i] <= (float)(0.0))
         {
         fprintf(stderr,"**** Problem with zero hz, check '%s'\n",gout);
         fprintf(stderr,"     exiting...\n");
         exit(-1);
         }
      }

   fclose(fpw);
   }
}

void setcoef(struct fdcf *cg0,struct fdcf *cg1,float *hg,int ng)
{
float dd[4];
double am[16], bv[4];
int i, ig;
int sgn[4];

for(ig=0;ig<ng;ig++)
   {
   for(i=0;i<4;i++)
      {
      cg0[ig].c4[i] = 0.0;
      cg1[ig].c4[i] = 0.0;
      }

   for(i=0;i<2;i++)
      {
      cg0[ig].c2[i] = 0.0;
      cg1[ig].c2[i] = 0.0;
      }
   }

sgn[0] = 1;
sgn[1] = 1;
sgn[2] = -1;
sgn[3] = -1;

for(ig=1;ig<ng-1;ig++)
   {
   dd[0] = hg[ig+1] + 0.5*hg[ig];
   dd[1] = 0.5*hg[ig];
   dd[2] = 0.5*hg[ig];
   dd[3] = hg[ig-1] + 0.5*hg[ig];

   for(i=0;i<4;i++)
      {
      am[i]    = 1.0;
      am[i+4]  = sgn[i]*dd[i]*am[i];
      am[i+8]  = sgn[i]*dd[i]*am[i+4];
      am[i+12] = sgn[i]*dd[i]*am[i+8];
      }

   bv[0] = 0.0;
   bv[1] = 1.0;
   bv[2] = 0.0;
   bv[3] = 0.0;

   gelim_double(am,4,bv);

   cg0[ig].c4[0] = bv[0];
   cg0[ig].c4[1] = bv[1];
   cg0[ig].c4[2] = bv[2];
   cg0[ig].c4[3] = bv[3];
   }

for(ig=1;ig<ng-2;ig++)
   {
   dd[0] = hg[ig+1] + 0.5*hg[ig+2];
   dd[1] = 0.5*hg[ig+1];
   dd[2] = 0.5*hg[ig];
   dd[3] = hg[ig] + 0.5*hg[ig-1];

   for(i=0;i<4;i++)
      {
      am[i] = 1.0;
      am[i+4] = sgn[i]*dd[i]*am[i];
      am[i+8] = sgn[i]*dd[i]*am[i+4];
      am[i+12] = sgn[i]*dd[i]*am[i+8];
      }

   bv[0] = 0.0;
   bv[1] = 1.0;
   bv[2] = 0.0;
   bv[3] = 0.0;

   gelim_double(am,4,bv);

   cg1[ig].c4[0] = bv[0];
   cg1[ig].c4[1] = bv[1];
   cg1[ig].c4[2] = bv[2];
   cg1[ig].c4[3] = bv[3];
   }

for(ig=0;ig<ng-1;ig++)
   {
   cg0[ig].c2[0] = 1.0/hg[ig];
   cg0[ig].c2[1] = -cg0[ig].c2[0];

   cg1[ig].c2[0] = 2.0/(hg[ig+1] + hg[ig]);
   cg1[ig].c2[1] = -cg1[ig].c2[0];
   }
}

void getlens(FILE *fp,char *s,float *len,float *g0,float *g1,float *dg,int *n)
{
int i = 0;

fgets(s,512,fp);
sscanf(s,"%f %f %f",&g0[0],&g1[0],&dg[0]);
g0[0] = 0.0;  /* force to be origin */

while(g1[i] < (*len))
   {
   i++;
   fgets(s,512,fp);
   sscanf(s,"%f %f %f",&g0[i],&g1[i],&dg[i]);
   }

*n = i+1;
g1[i] = (*len);  /* force to be total length */
}

/*
Gauss elimination without pivoting (no row exchanges):
We solve Ax = b where A is n by n, and x and b have length n
by forming the decomposition A = LU via elimination.
Originally a contains the matrix A and b contains the
vector b.  At the end a contains the lower and upper
triangular matrices L and U and b contains the solution
vector x.  The diagonal of a contains the diagonal of U
(the pivots) since the diagonal elements of U are all 1's.
*/

void gelim_double(double *a,int n,double *b)
   {
        int i, j, jj;
        double *pa, *paj;
        double f, pivot;

         for (j=0; j<n-1; j++)   /* lu decomp of a (no row exchanges) */
            {
                pa= a + j*n;
                pivot = pa[j];
                for (i=j+1; i<n; i++)
                   {
                        pa = a + i*n;
                        paj= a + j*n;
                        f = pa[j]/pivot;
                        pa[j] = f;
                        for (jj=j+1; jj<n; jj++) pa[jj] -= f*paj[jj];
                   }
           }
         for (i=1; i<n; i++)        /* forward elimination on b */
            {
                pa = a + i*n;
                for (j=0; j<i; j++) b[i] -= pa[j]*b[j];
            }
         for (i=n-1; i>-1; i--)        /* back-substitution */
            {
                pa = a + i*n;
                for (j=n-1; j>i; j--) b[i] -= pa[j]*b[j];
                b[i] = b[i]/pa[i];
            }
    }

void check_vmax(int fdp,int fds,int fdd,float *vp,float *vs,float *den,struct gridparam *gp,int *ixx,int *iyx,int *izx,float *vx,float *hn,float *dtx)
{
int size_plane;
int j, ix, iy, iz, nx, ny, nz;
float f, h, amax, bmax, dmax;

nx = gp->nx;
ny = gp->ny;
nz = gp->nz;

size_plane = nx*nz*sizeof(float);

*dtx = 1.0e+15;
for(iy=0;iy<ny;iy++)
   {
   reed(fdp,(char *)vp,size_plane);
   reed(fds,(char *)vs,size_plane);
   reed(fdd,(char *)den,size_plane);

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
         j = ix + iz*nx;

         h = 1.0/sqrt(1.0/(gp->hx[ix]*gp->hx[ix]) + 1.0/(gp->hy[iy]*gp->hy[iy]) + 1.0/(gp->hz[iz]*gp->hz[iz]));

         f = h/vp[j];
         if(f < *dtx)
            {
            amax = vp[j];
            bmax = vs[j];
            dmax = den[j];

            *ixx = ix;
	    *iyx = iy;
	    *izx = iz;
            *vx = vp[j];
            *hn = sqrt(3.0)*h;
            *dtx = f;
            }
         }
      }
   }
lseek(fdp,0,0);
lseek(fds,0,0);
lseek(fdd,0,0);
}

void set_minmax(float *x,float *min,float *max,int n)
{
int i;

for(i=0;i<n;i++)
   {
   if(x[i] < *min)
      x[i] = *min;
   if(x[i] > *max)
      x[i] = *max;
   }
}

void set_vpvs(float *vp,float *vs,int n)
{
int i;
float fac = 1.45;   /* makes poisson's ratio about 0.05 */

for(i=0;i<n;i++)
   {
   if(vp[i]/vs[i] < fac)
      vs[i] = vp[i]/fac;
   }
}

void check_vmin(int fdp,int fds,int fdd,float *vp,float *vs,float *den,struct gridparam *gp)
{
int size_plane;
int j, ix, iy, iz, nx, ny, nz;
int ax, ay, az;
int bx, by, bz;
int dx, dy, dz;
int rx, ry, rz;
float f, h, amin, bmin, dmin, rmax;
float vpr, vsr;

nx = gp->nx;
ny = gp->ny;
nz = gp->nz;

size_plane = nx*nz*sizeof(float);

rmax = -1.0e+15;
amin = 1.0e+15;
bmin = 1.0e+15;
dmin = 1.0e+15;
for(iy=0;iy<ny;iy++)
   {
   reed(fdp,(char *)vp,size_plane);
   reed(fds,(char *)vs,size_plane);
   reed(fdd,(char *)den,size_plane);

   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
         j = ix + iz*nx;

         if(vp[j]/vs[j] > rmax)
	    {
	    rmax = vp[j]/vs[j];
	    vpr = vp[j];
	    vsr = vs[j];
            rx = ix;
            ry = iy;
            rz = iz;
	    }

         if(vp[j] < amin)
	    {
	    amin = vp[j];
            ax = ix;
            ay = iy;
            az = iz;
	    }

         if(vs[j] < bmin)
	    {
	    bmin = vs[j];
            bx = ix;
            by = iy;
            bz = iz;
	    }

         if(den[j] < dmin)
	    {
	    dmin = den[j];
            dx = ix;
            dy = iy;
            dz = iz;
	    }

         if(isnan(vp[j]))
	    fprintf(stderr,"vp= %.5f ix=%5d iy=%5d iz=%5d\n",vp[j],ix,iy,iz);

         if(isnan(vs[j]))
	    fprintf(stderr,"vs= %.5f ix=%5d iy=%5d iz=%5d\n",vs[j],ix,iy,iz);

         if(isnan(den[j]))
	    fprintf(stderr,"den= %.5f ix=%5d iy=%5d iz=%5d\n",den[j],ix,iy,iz);

	 if(vp[j]/vs[j] < sqrt(2.0))
            fprintf(stderr,"vp= %.5f vs= %.5f vp/vs= %.5f ix=%5d iy=%5d iz=%5d\n",vp[j],vs[j],vp[j]/vs[j],ix,iy,iz);

	 fflush(stderr);
         }
      }
   }

fprintf(stderr,"vp/vs= %.5f vp= %.5f vs= %.5f ix=%5d iy=%5d iz=%5d\n",rmax,vpr,vsr,rx,ry,rz);
fprintf(stderr,"vpmin= %.5f ix=%5d iy=%5d iz=%5d\n",amin,ax,ay,az);
fprintf(stderr,"vsmin= %.5f ix=%5d iy=%5d iz=%5d\n",bmin,bx,by,bz);
fprintf(stderr,"dnmin= %.5f ix=%5d iy=%5d iz=%5d\n",dmin,dx,dy,dz);
fflush(stderr);

lseek(fdp,0,0);
lseek(fds,0,0);
lseek(fdd,0,0);
}
