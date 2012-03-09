#include "include.h"

void gcproj(float *xf,float *yf,float *rlon,float *rlat,float *ref_rad,double *g0,double *b0,double *amat,double *ainv,int gflag)
{
double xp, yp, zp;
double xg, yg, zg;
double arg;
double cosG, sinG;
double cosB, sinB;

double rperd = RPERD;

if(gflag == 0)
   {
   arg = (*xf)/(*ref_rad) - (*b0);
   cosB = cos(arg);
   sinB = sin(arg);

   arg = (*yf)/(*ref_rad) - (*g0);
   cosG = cos(arg);
   sinG = sin(arg);

   arg = sqrt(1.0 + sinB*sinB*sinG*sinG);
   xp = sinG*cosB*arg;
   yp = sinB*cosG*arg;
   zp = sqrt(1.0 - xp*xp - yp*yp);

   xg = xp*amat[0] + yp*amat[1] + zp*amat[2];
   yg = xp*amat[3] + yp*amat[4] + zp*amat[5];
   zg = xp*amat[6] + yp*amat[7] + zp*amat[8];

   arg = sqrt(xg*xg + yg*yg)/zg;
   (*rlat) = 90.0 - atan(arg)/rperd;

   if(xg != (double)(0.0))
      {
      arg = yg/xg;
      (*rlon) = atan(arg)/rperd;
      }
   else
      (*rlon) = 0.0;

   /*
   RWG 05/08/08
   This is incorrect.  Replaced with following conditional on 'xg'.
   if(yg < (double)(0.0))
   */
   if(xg < (double)(0.0))
      (*rlon) = (*rlon) - 180.0;

   while((*rlon) < (double)(-180.0))
      (*rlon) = (*rlon) + 360.0;
   }
else
   {
   arg = (*rlon)*rperd;
   cosG = cos(arg);
   sinG = sin(arg);

   arg = (90.0 - (*rlat))*rperd;
   cosB = cos(arg);
   sinB = sin(arg);

   xg = sinB*cosG;
   yg = sinB*sinG;
   zg = cosB;

   xp = xg*ainv[0] + yg*ainv[1] + zg*ainv[2];
   yp = xg*ainv[3] + yg*ainv[4] + zg*ainv[5];
   zp = xg*ainv[6] + yg*ainv[7] + zg*ainv[8];

   sinG = xp/sqrt(1.0 - yp*yp);
   sinB = yp/sqrt(1.0 - xp*xp);

   *xf = (double)(*ref_rad)*(asin(sinB)+(*b0));
   *yf = (double)(*ref_rad)*(asin(sinG)+(*g0));
   }
}

void gen_matrices(double *amat,double *ainv,float *alpha,float *ref_lon,float *ref_lat)
{
double arg;
double cosA, sinA;
double cosT, sinT;
double cosP, sinP;
double det;

double rperd = RPERD;

arg = (double)(*alpha)*rperd;
cosA = cos(arg);
sinA = sin(arg);

arg = (double)(90.0-*ref_lat)*rperd;
cosT = cos(arg);
sinT = sin(arg);

arg = (double)(*ref_lon)*rperd;
cosP = cos(arg);
sinP = sin(arg);

amat[0] = cosA*cosT*cosP + sinA*sinP;
amat[1] = sinA*cosT*cosP - cosA*sinP;
amat[2] = sinT*cosP;
amat[3] = cosA*cosT*sinP - sinA*cosP;
amat[4] = sinA*cosT*sinP + cosA*cosP;
amat[5] = sinT*sinP;
amat[6] = -cosA*sinT;
amat[7] = -sinA*sinT;
amat[8] = cosT;

det = amat[0]*(amat[4]*amat[8] - amat[7]*amat[5])
    - amat[1]*(amat[3]*amat[8] - amat[6]*amat[5])
    + amat[2]*(amat[3]*amat[7] - amat[6]*amat[4]);

det = 1.0/det;
ainv[0] = det*amat[0];
ainv[1] = det*amat[3];
ainv[2] = det*amat[6];
ainv[3] = det*amat[1];
ainv[4] = det*amat[4];
ainv[5] = det*amat[7];
ainv[6] = det*amat[2];
ainv[7] = det*amat[5];
ainv[8] = det*amat[8];
}

void geocen(float *r,double x)
{
*r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
fprintf(stderr,"%20.10f %20.10f %20.10f\n",*r,x,atan((1.0 - (1.0/FLAT_CONST))*tan(x)));
}

void set_g2(float *g2,float *fc)
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

void latlon2km(float *arg,float *latkm,float *lonkm,float *rc,float *g2)
{
float cosA, sinA, g2s2, den;

float fone = 1.0;
float ftwo = 2.0;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((fone)/((fone) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((fone) + g2s2*((ftwo) + (*g2))))*den*den*den;
}

void init_modelc(struct runparams *rp)
{
FILE *fpr;
char string[512];
float lon, lat;
int ip, ix, iy, localny, iyleft, iyright;

double rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg, xlen, ylen, kmlon, kmlat;
float cosR, sinR, xx, yy, ar1, ar2, ar3, ar4; 

double g0, b0;
double amat[9], ainv[9];
int xy2ll = 0;
int ll2xy = 1;

xlen = (rp->nx)*(rp->h);
ylen = (rp->globny)*(rp->h);

iyleft = rp->ny1 + 2;
if(rp->ny1 == 0) iyleft = 0;

iyright = rp->ny2 - 3;
if(rp->ny2 == rp->globny) iyright = rp->ny2 - 1;

localny = iyright - iyleft + 1;

if(rp->geoproj == 0)
   {
   if(rp->xshift < 0.0 && rp->xshift > -xlen && rp->yshift < 0.0 && rp->yshift > -ylen)
      rp->center_origin = 1;
   else
      {
      rp->xshift = 0.0;
      rp->yshift = 0.0;
      rp->center_origin = 0;
      }

   rp->cosR = cos(rp->modelrot*RPERD);
   rp->sinR = sin(rp->modelrot*RPERD);

   radc = ERAD*RPERD;
   set_g2(&g2,&fc);

   latavg = rp->modellat;
   if(rp->center_origin == 0)  /* backward compatible */
      latavg = rp->modellat - 0.5*(xlen*rp->sinR + ylen*rp->cosR)/111.20;

   geocen(&latavg,(double)(latavg)*rperd);
   latlon2km(&latavg,&(rp->kmlat),&(rp->kmlon),&radc,&g2);

   ar1 = (rp->cosR)/(rp->kmlon);
   ar2 = (rp->sinR)/(rp->kmlon);
   ar3 = (rp->cosR)/(rp->kmlat);
   ar4 = (rp->sinR)/(rp->kmlat);

   for(iy=iyleft;iy<=iyright;iy++)
      {
      for(ix=0;ix<rp->nx;ix++)
         {
         xx = ix*rp->h;
         yy = iy*rp->h;
         
	 lon = rp->modellon + (xx + rp->xshift)*ar1 - (yy + rp->yshift)*ar2;
         lat = rp->modellat - (xx + rp->xshift)*ar4 - (yy + rp->yshift)*ar3;

         ip = ix + (iy - iyleft)*rp->nx;
         }
      }
   }
else if(rp->geoproj == 1)
   {
   rp->center_origin = 1;
   rp->xshift = -0.5*xlen;
   rp->yshift = -0.5*ylen;
   rp->erad = erad;

   gen_matrices(rp->amat,rp->ainv,&rp->modelrot,&rp->modellon,&rp->modellat);

   rp->g0 = (double)(0.5*ylen)/(double)(erad);
   rp->b0 = (double)(0.5*xlen)/(double)(erad);

   for(iy=iyleft;iy<=iyright;iy++)
      {
      for(ix=0;ix<rp->nx;ix++)
         {
         xx = ix*rp->h;
         yy = iy*rp->h;

         gcproj(&xx,&yy,&lon,&lat,&(rp->erad),&(rp->g0),&(rp->b0),rp->amat,rp->ainv,xy2ll);

         ip = ix + (iy - iyleft)*rp->nx;
         }
      }
   }
}

//Non-geoproj stuff
makedir(ipath)
char *ipath;
{
struct stat sbuf;
char stmp[256], str[128], path[1024];
int rtn, j;
mode_t mode = 00777;

strcpy(path,ipath);

j = 0;
while(path[j] != '\0')
   j++;

j--;
while(path[j] == '/')
   j--;
path[j+1] = '\0';
path[j+2] = '\0';

j = 0;
while(path[j] != '\0')
   {
   while(path[j] != '/' && path[j] != '\0')
      j++;

   if(j != 0)
      {
      strncpy(stmp,path,j);
      stmp[j] = '\0';

      rtn = stat(stmp,&sbuf); /* stat directory path to see if it already exists */

      if(rtn == -1 && errno == ENOENT) /* try to make the directory path */
         {
         rtn = mkdir(stmp,mode);

         if(rtn == -1)
            {
            if(errno != EEXIST)
               {
               sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
               perror(str);
               mpi_exit(-1);
               }
            }
         }

      else if(rtn == -1 && errno != ENOENT) /* some other problem */
         {
         sprintf(str,"problem with stat() on %s, exiting",stmp);
         perror(str);
         mpi_exit(-1);
         }
      }
   j++;
   }

rtn = mkdir(stmp,mode);
if(rtn == -1 && errno != EEXIST)
   {
   sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
   perror(str);
   mpi_exit(-1);
   }

return 0;
}

void mpi_exit(int val)
{
MPI_Finalize();
exit(val);
}

FILE *fopfile(name,mode)
char *name, *mode;
{
FILE *fp, *fopen();
int j;
char path[512];

if((fp = fopen(name,mode)) == NULL && errno == ENOENT)
   {
   j = strlen(name);
   while(name[j] != '/' && j > 1)
      j--;

   strncpy(path,name,j);
   path[j] = '\0';
   makedir(path);
   }

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

opfile_ro (name)
char *name;
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

void *check_malloc(size_t len)
{
void *ptr;

ptr = (void *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   mpi_exit(-1);
   }

return(ptr);
}

size_t reed(int fd, void *pntr, size_t length)
{
size_t temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%d attempted  %d read\n", length, temp);
   exit(-1);
   }
return(temp);
}

size_t rite(int fd, void *pntr, size_t length)
{
size_t temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%u attempted  %u written\n", length, temp);
   exit(-1);
   }
return(temp);
}

int cropfile_rw(char *name)
{
int fd;
if ((fd = open(name,CROP_FLAGS,0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return(fd);
}

