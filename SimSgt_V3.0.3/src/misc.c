/*
   misc.c comtains the following functions:
      
      copy()
      copy_sgn()
      copyusage()
      diftime()
      init_field()
      makedir()
      mirror()
      print_report()
      zero()
*/

#include "include.h"

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

	 /*
fprintf(stderr,"**** errno=%d path=%s\n",errno,stmp);
fprintf(stderr,"\tEACCES=%d\n \tEEXIST=%d\n \tEFAULT=%d\n \tEIO=%d\n \tELOOP=%d\n \tEMLINK=%d\n \tEMULTIHOP=%d\n \tENAMETOOLONG=%d\n \tENOENT=%d\n \tENOLINK=%d\n \tENOSPC=%d\n \tENOTDIR=%d\n \tEROFS=%d\n\n",EACCES,
EEXIST,
EFAULT,
EIO,
ELOOP,
EMLINK,
EMULTIHOP,
ENAMETOOLONG,
ENOENT,
ENOLINK,
ENOSPC,
ENOTDIR,
EROFS);
fflush(stderr);
	 */

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

/*
   Double check to make sure directory exists.  This is a brute-force
   method, but I ran inot problems with automounted directories using the
   error-checking above.  RWG 9/20/99
*/

rtn = mkdir(stmp,mode);
if(rtn == -1 && errno != EEXIST)
   {
   sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
   perror(str);
   mpi_exit(-1);
   }
}

void floor_c(float *c,float *r, int n)
{
int i;

for(i=0;i<n;i++)
   c[i] = c[i] + r[i];
}

void init_field_val(float v, float *x,int n)
{
int i;
double frand();

for(i=0;i<n;i++)
   x[i] = v*frand();
}

init_field(x,n)
register float *x;
int n; 
{
int i;
double frand();

for(i=0;i<n;i++)
   x[i] = 1.0e-15*frand();

   /*
   x[i] = 0.0;
   */
}

zero(x,n)
register float *x;
int n;
{
float zero = 0.0;
int i;

for(i=0;i<n;i++)
   x[i] = zero;
}

copy_sgn(q1,q2,n,sgn)
register float *q1, *q2;
int n, sgn;
{
int i;

for(i=0;i<n;i++)
   q2[i] = sgn*q1[i];
}

copy(q1,q2,n,stride)
register float *q1, *q2;
int n, stride;
{
int i;

for(i=0;i<n;i++)
   q2[i] = q1[i*stride];
}

add_random_field(pv,rf,nx,nz)
float *pv, *rf;
int nx, nz;
{
float *vx, *vy, *vz;
float *txx, *tyy, *tzz, *txy, *txz, *tyz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
int i;

vx = pv;
vy = pv + nx*nz;
vz = pv + 2*nx*nz;
txx = pv + 3*nx*nz;
tyy = pv + 4*nx*nz;
tzz = pv + 5*nx*nz;
txy = pv + 6*nx*nz;
txz = pv + 7*nx*nz;
tyz = pv + 8*nx*nz;
cxx = pv + 9*nx*nz;
cyy = pv + 10*nx*nz;
czz = pv + 11*nx*nz;
cxy = pv + 12*nx*nz;
cxz = pv + 13*nx*nz;
cyz = pv + 14*nx*nz;

for(i=0;i<nx*nz;i++)
   {
   vx[0] = vx[0] + rf[0];
   vy[0] = vy[0] + rf[0];
   vz[0] = vz[0] + rf[0];
   txx[0] = txx[0] + rf[0];
   tyy[0] = tyy[0] + rf[0];
   tzz[0] = tzz[0] + rf[0];
   txy[0] = txy[0] + rf[0];
   txz[0] = txz[0] + rf[0];
   tyz[0] = tyz[0] + rf[0];
   cxx[0] = cxx[0] + rf[0];
   cyy[0] = cyy[0] + rf[0];
   czz[0] = czz[0] + rf[0];
   cxy[0] = cxy[0] + rf[0];
   cxz[0] = cxz[0] + rf[0];
   cyz[0] = cyz[0] + rf[0];

   vx++; vy++; vz++; rf++;
   txx++; tyy++; tzz++; txy++; txz++; tyz++;
   cxx++; cyy++; czz++; cxy++; cxz++; cyz++;
   }
}

void *check_realloc(void *ptr,size_t len)
{

ptr = (void *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   mpi_exit(-1);
   }

return(ptr);
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

static	long	frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand()
   {
	frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
	return( (double)(frandx)/1073741824.0 -1.0 );
   }

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(seed)
int *seed;
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

dtimes(tu1,rt1,dtu,drt,ts2,len)
struct tms tu1, *dtu;
size_t len;
clock_t rt1, *drt;
float *ts2;
{
struct tms tu2;
clock_t rt2;

rt2 = times(&tu2);

dtu->tms_utime = dtu->tms_utime + tu2.tms_utime - tu1.tms_utime;
dtu->tms_stime = dtu->tms_stime + tu2.tms_stime - tu1.tms_stime;
dtu->tms_cutime = dtu->tms_utime + tu2.tms_cutime - tu1.tms_cutime;
dtu->tms_cstime = dtu->tms_cstime + tu2.tms_cstime - tu1.tms_cstime;

*drt = *drt + rt2 - rt1;
*ts2 = *ts2 + (float)(len)*1.0e-06;
}

init_report(rt0,rt1,u2,u1,dtu,drt,ts2,ts1,report)
struct tms *u1, *u2, *dtu;
float *ts1, *ts2;
int report;
clock_t rt0, *rt1, *drt;
{
fprintf(stderr,"\nStatus given at initialization and every %d time steps\n",report);
fprintf(stderr,"All times in seconds\n\n");
fprintf(stderr,"                  Data Transfer      Transfer & Computation      Cumulative\n");
fprintf(stderr,"time step       CPU    Mbyte  %%Real        CPU  %%Real            CPU  %%Real\n");

dtu->tms_utime = 0;
dtu->tms_stime = 0;
dtu->tms_cutime = 0;
dtu->tms_cstime = 0;
*drt = 0;
*ts2 = 0.0;
*ts1 = 0.0;

print_report(rt0,rt1,u2,u1,dtu,drt,ts2,ts1,0);
}

print_report(rt0,rt1,u2,u1,dtu,drt,ts2,ts1,it)
struct tms *u1, *u2, *dtu;
float *ts1, *ts2;
int it;
clock_t rt0, *rt1, *drt;
{
clock_t rt2;
float uos, uou, uor;
float tos, tou, tor;

tou = (float)(dtu->tms_utime)/(float)(sysconf(_SC_CLK_TCK));
tos = (float)(dtu->tms_stime)/(float)(sysconf(_SC_CLK_TCK));
tor = (double)(*drt)/(double)(sysconf(_SC_CLK_TCK)) + 1.0e-05;

rt2 = times(u2);
uou = (float)(u2->tms_utime - u1->tms_utime)/(float)(sysconf(_SC_CLK_TCK));
uos = (float)(u2->tms_stime - u1->tms_stime)/(float)(sysconf(_SC_CLK_TCK));
uor = (double)(rt2 - *rt1)/(double)(sysconf(_SC_CLK_TCK));

if(it < 0)
   fprintf(stderr,"         ");
else
   fprintf(stderr,"%9d",it);

fprintf(stderr,"%10.2f",tos + tou);
fprintf(stderr,"%9.2f",(*ts2)-(*ts1));
fprintf(stderr,"%7.2f",(tos + tou)/(tor));
fprintf(stderr,"%11.2f",uos + uou);
fprintf(stderr,"%7.2f",(uos + uou)/(uor));
fprintf(stderr,"%14.0f.",(float)(u2->tms_utime + u2->tms_stime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%7.2f\n",(double)(u2->tms_utime + u2->tms_stime)/(double)(rt2-rt0));
fflush(stderr);

u1->tms_utime = u2->tms_utime;
u1->tms_stime = u2->tms_stime;
u1->tms_cutime = u2->tms_cutime;
u1->tms_cstime = u2->tms_cstime;
*rt1 = rt2;

dtu->tms_utime = 0;
dtu->tms_stime = 0;
dtu->tms_cutime = 0;
dtu->tms_cstime = 0;
*drt = 0;
*ts1 = *ts2;
}

Xprint_report(rt1,u2,u1,it)
struct tms *u1, *u2;
int it, *rt1;
{
int rt2;
float uos, uou, uor;

rt2 = times(u2);
uou = (float)(u2->tms_utime - u1->tms_utime)/(float)(sysconf(_SC_CLK_TCK));
uos = (float)(u2->tms_stime - u1->tms_stime)/(float)(sysconf(_SC_CLK_TCK));
uor = (float)(rt2 - *rt1)/(float)(sysconf(_SC_CLK_TCK));

if(it < 0)
   fprintf(stderr,"            ");
else
   fprintf(stderr,"  %7d   ",it);

fprintf(stderr,"%7.3f  ",uos);
fprintf(stderr,"%7.3f   ",uou);
fprintf(stderr,"%7.0f.    ",(float)(u2->tms_stime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%7.0f.    ",(float)(u2->tms_utime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%7.0f.    ",(float)(u2->tms_utime + u2->tms_stime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%6.2f\n",(uos + uou)/(uor));

u1->tms_utime = u2->tms_utime;
u1->tms_stime = u2->tms_stime;
u1->tms_cutime = u2->tms_cutime;
u1->tms_cstime = u2->tms_cstime;
*rt1 = rt2;
}

total_report(rt1,u2,tsize)
struct tms *u2;
float tsize;
clock_t rt1;
{
clock_t rt2;

rt2 = times(u2);

fprintf(stderr,"\n*** usage totals:     Mbytes Transfered    System CPU      User CPU   %%Real\n");
fprintf(stderr,"                      ");
fprintf(stderr,"%17.2f",tsize);
fprintf(stderr,"%13.0f.",(float)(u2->tms_stime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%13.0f.",(float)(u2->tms_utime)/(float)(sysconf(_SC_CLK_TCK)));
fprintf(stderr,"%8.2f\n",(double)(u2->tms_utime + u2->tms_stime)/(double)(rt2-rt1));
fflush(stderr);
}

void copy_swap(char *obuf,char *ibuf,int n)
{

while(n--)
   {
   obuf[3] = ibuf[0];
   obuf[2] = ibuf[1];
   obuf[1] = ibuf[2];
   obuf[0] = ibuf[3];

   obuf = obuf + 4;
   ibuf = ibuf + 4;
   }
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

void swap_in_place_d(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[7];
   cbuf[7] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[6];
   cbuf[6] = cv;

   cv = cbuf[2];
   cbuf[2] = cbuf[5];
   cbuf[5] = cv;

   cv = cbuf[3];
   cbuf[3] = cbuf[4];
   cbuf[4] = cv;

   cbuf = cbuf + 8;
   }
}

void set_tmpdir(char *name,int rank,char *list,char *dir)
{
FILE *fpr, *fopfile();
int iflag, id;
char pn[256], td[512], str[1024];

fpr = fopfile(list,"r");

iflag = 0;
while(fgets(str,1024,fpr) != NULL)
   {
   sscanf(str,"%s %d %s",pn,&id,td);
   if(strcmp(pn,name) == 0 && id == rank)
      {
      iflag = 1;
      sprintf(dir,"%s",td);
      }
   }
if(iflag == 0)
   {
   fprintf(stderr,"**** ERROR: Node '%s %d' not found in tmpdirlist= %s\n",name,rank,list);
   mpi_exit(-1);
   }
}

void check_procname(char *name,int *len)
{
int l;

if(name[*len-1] == '.')
   name[*len-1] = '\0';

*len = *len - 1;
}

void get_nodetype_neighbor(struct nodeinfo *ni)
{
int ip, ipx, ipy, ipz;
int xv, yv, zv;

ip = 0;
for(ipy=0;ipy<ni->nproc_y;ipy++)
   {
/*
   if(ipy == 0)
      yv = NODE_MIN_BOUNDARY;
   else if(ipy == (ni->nproc_y - 1))
      yv = NODE_MAX_BOUNDARY;
   else
      yv = NODE_INTERIOR;
*/

   for(ipz=0;ipz<ni->nproc_z;ipz++)
      {
/*
      if(ipz == 0)
         zv = NODE_MIN_BOUNDARY;
      else if(ipz == (ni->nproc_z - 1))
         zv = NODE_MAX_BOUNDARY;
      else
         zv = NODE_INTERIOR;
*/

      for(ipx=0;ipx<ni->nproc_x;ipx++)
         {
/*
         if(ipx == 0)
            xv = NODE_MIN_BOUNDARY;
         else if(ipx == (ni->nproc_x - 1))
            xv = NODE_MAX_BOUNDARY;
         else
            xv = NODE_INTERIOR;
*/

	 if(ni->segmentId == ip)
	    {
	    ni->procId_x = ipx;
	    ni->procId_y = ipy;
	    ni->procId_z = ipz;

	    ni->minusId_x = ip-1;
	    ni->plusId_x  = ip+1;
	    ni->minusId_z = ip-(ni->nproc_x);
	    ni->plusId_z  = ip+(ni->nproc_x);
	    ni->minusId_y = ip-(ni->nproc_x)*(ni->nproc_z);
	    ni->plusId_y  = ip+(ni->nproc_x)*(ni->nproc_z);
	    
            if(ipx == 0)
	       ni->minusId_x = -1;
            if(ipx == (ni->nproc_x - 1))
	       ni->plusId_x = -1;
	    
            if(ipy == 0)
	       ni->minusId_y = -1;
            if(ipy == (ni->nproc_y - 1))
	       ni->plusId_y = -1;
	    
            if(ipz == 0)
	       ni->minusId_z = -1;
            if(ipz == (ni->nproc_z - 1))
	       ni->plusId_z = -1;

/*
	    if(xv == NODE_MIN_BOUNDARY)
	       ni->minusId_x = -1;
	    if(xv == NODE_MAX_BOUNDARY)
	       ni->plusId_x = -1;

	    if(yv == NODE_MIN_BOUNDARY)
	       ni->minusId_y = -1;
	    if(yv == NODE_MAX_BOUNDARY)
	       ni->plusId_y = -1;

	    if(zv == NODE_MIN_BOUNDARY)
	       ni->minusId_z = -1;
	    if(zv == NODE_MAX_BOUNDARY)
	       ni->plusId_z = -1;
*/
	    }

	 ip++;
         }
      }
   }
}

/*
   Maximize volume to surface area by making the volumes as close to cubes as possible,
   that is, attempt to make loc_nx=loc_ny=loc_nz=D (constant).

   This means

      D*nproc_x=nx
      D*nproc_y=ny
      D*nproc_z=nz

   or

      (D*nproc_x)*(D*nproc_y)*(D*nproc_z) = D^3(nx*ny*nz)

   and using nproc=nproc_x*nproc_y*nproc_z, we get

      D = ((nx*ny*nz)/nproc)^(1/3)

   which then gives

      nproc_x = nx*(nproc/(nx*ny*nz))^(1/3)
      nproc_y = ny*(nproc/(nx*ny*nz))^(1/3)
      nproc_z = nz*(nproc/(nx*ny*nz))^(1/3)

   To make sure the integers align and we get the right total nprocs, solve this using
   (start with nproc_z since nz is probably the smallest -> strongest constraint on D).

      nproc_z = nz*(nproc/(nx*ny*nz))^(1/3)
      nproc_x = nx*(nproc/(nx*ny*nproc_z)^(1/2)
      nproc_y = nproc/(nproc_x*nproc_z)

*/

int get_nproc(struct nodeinfo *ni)
{
int np;
int ipt, ip2, ip3, np2;
float fnp, fnx, fny, fnz, inv_fmp;

inv_fmp = 1.0/(float)(ni->min_nproc);
fnp = (float)(ni->nproc);
fnx = (float)(ni->globnx);
fny = (float)(ni->globny);
fnz = (float)(ni->globnz);

if((ni->nproc_z) < 0)
   {
   ni->nproc_z = (int)(inv_fmp*fnz*exp(log(fnp/(fnx*fny*fnz))/3.0) + 0.5);
   if((ni->nproc_z) < 1)
      ni->nproc_z = 1;
   ni->nproc_z = (ni->min_nproc)*(ni->nproc_z);
   }

if((ni->nproc_x) < 0)
   {
   ni->nproc_x = (int)(inv_fmp*fnx*exp(log(fnp/(fnx*fny*(ni->nproc_z)))/2.0) + 0.5);
   if((ni->nproc_x) < 1)
      ni->nproc_x = 1;
   ni->nproc_x = (ni->min_nproc)*(ni->nproc_x);
   }

ni->nproc_y = (int)(inv_fmp*fnp/((float)(ni->nproc_x)*(float)(ni->nproc_z)) + 0.5);
if((ni->nproc_y) < 1)
   ni->nproc_y = 1;
ni->nproc_y = (ni->min_nproc)*(ni->nproc_y);

np = (ni->nproc_x)*(ni->nproc_y)*(ni->nproc_z);

if(np != ni->nproc)
   {
   ip3 = (int)(exp(log(fnp)/3.0) + 0.5);

   ipt = 1;
   while(2*ipt <= ip3 && (ni->nproc)%ipt == 0 && (ni->nproc/ipt)%2 == 0)
      ipt = 2*ipt;

   ni->nproc_z = ipt;

   np2 = ni->nproc/ni->nproc_z;
   ip2 = (int)(exp(log(1.0*np2)/2.0) + 0.5);

   ipt = 1;
   while(2*ipt <= ip2 && (np2)%ipt == 0 && (np2/ipt)%2 == 0)
      ipt = 2*ipt;

   ni->nproc_x = ipt;
   ni->nproc_y = np2/ni->nproc_x;

   np = (ni->nproc_x)*(ni->nproc_y)*(ni->nproc_z);
   }

return(np);
}

void get_n1n2(int flag,struct nodeinfo *ni)
{
int i;
float fslice, fn1, fn2;

/* find nx1,nx2,nx
        ny1,ny2,ny
        nz1,nz2,nz

   flag=0 -> balance only on interior nodes
   flag=1 -> balance across all nodes
*/

if(flag == 0)
   {
   fslice = (float)(ni->globnx)/(float)(ni->nproc_x) + 3.0;
   fn1 = -2.0;        /* start at -2 so 1st n2 is i=fslice+1  */
   }
else
   {
   fslice = (float)(ni->globnx + (ni->nproc_x-1.0)*4.0)/(float)(ni->nproc_x) - 1.0;
   fn1 = 0.0;
   }

for(i=0;i<=ni->procId_x;i++)
   {
   fn2 = fn1 + fslice;

   ni->nx1 = (int)(fn1 + 0.5);
   ni->nx2 = (int)(fn2 + 0.5);

   fn1 = fn2 - 3.0;
   }

ni->ixminus = ni->nx1 + 2;
if(ni->minusId_x < 0)
   {
   ni->nx1 = 0;
   ni->ixminus = 0;
   }

ni->ixplus = ni->nx2 - 2;
if(ni->plusId_x < 0)
   {
   ni->nx2 = ni->globnx - 1;
   ni->ixplus = ni->globnx - 1;
   }

ni->nx2 = ni->nx2 + 1;
ni->loc_nx = ni->nx2 - ni->nx1;

if(flag == 0)
   {
   fslice = (float)(ni->globny)/(float)(ni->nproc_y) + 3.0;
   fn1 = -2.0;        /* start at -2 so 1st n2 is i=fslice+1  */
   }
else
   {
   fslice = (float)(ni->globny + (ni->nproc_y-1.0)*4.0)/(float)(ni->nproc_y) - 1.0;
   fn1 = 0.0;
   }

for(i=0;i<=ni->procId_y;i++)
   {
   fn2 = fn1 + fslice;

   ni->ny1 = (int)(fn1 + 0.5);
   ni->ny2 = (int)(fn2 + 0.5);

   fn1 = fn2 - 3.0;
   }

ni->iyminus = ni->ny1 + 2;
if(ni->minusId_y < 0)
   {
   ni->ny1 = 0;
   ni->iyminus = 0;
   }

ni->iyplus = ni->ny2 - 2;
if(ni->plusId_y < 0)
   {
   ni->ny2 = ni->globny - 1;
   ni->iyplus = ni->globny - 1;
   }

ni->ny2 = ni->ny2 + 1;
ni->loc_ny = ni->ny2 - ni->ny1;

if(flag == 0)
   {
   fslice = (float)(ni->globnz)/(float)(ni->nproc_z) + 3.0;
   fn1 = -2.0;        /* start at -2 so 1st n2 is i=fslice+1  */
   }
else
   {
   fslice = (float)(ni->globnz + (ni->nproc_z-1.0)*4.0)/(float)(ni->nproc_z) - 1.0;
   fn1 = 0.0;
   }

for(i=0;i<=ni->procId_z;i++)
   {
   fn2 = fn1 + fslice;

   ni->nz1 = (int)(fn1 + 0.5);
   ni->nz2 = (int)(fn2 + 0.5);

   fn1 = fn2 - 3.0;
   }

ni->izminus = ni->nz1 + 2;
if(ni->minusId_z < 0)
   {
   ni->nz1 = 0;
   ni->izminus = 0;
   }

ni->izplus = ni->nz2 - 2;
if(ni->plusId_z < 0)
   {
   ni->nz2 = ni->globnz - 1;
   ni->izplus = ni->globnz - 1;
   }

ni->nz2 = ni->nz2 + 1;
ni->loc_nz = ni->nz2 - ni->nz1;
}

void get_n1n2_indv(int flag,int gn,int np,int id,int minusId,int plusId,int *n1,int *n2,int *nn)
{
int i;
float fslice, fn1, fn2;

if(flag == 0)
   {
   fslice = (float)(gn)/(float)(np) + 3.0;
   fn1 = -2.0;        /* start at -2 so 1st n2 is i=fslice+1  */
   }
else
   {
   fslice = (float)(gn + (np-1.0)*4.0)/(float)(np) - 1.0;
   fn1 = 0.0;
   }

for(i=0;i<=id;i++)
   {
   fn2 = fn1 + fslice;

   *n1 = (int)(fn1 + 0.5);
   *n2 = (int)(fn2 + 0.5);

   fn1 = fn2 - 3.0;
   }

if(minusId < 0)
   *n1 = 0;

if(plusId < 0)
   *n2 = gn - 1;

*n2 = *n2 + 1;      /* SN: in 0 convention */
*nn = *n2 - *n1;    /* nn is now local number of planes */
}

void get_ny1ny2(int flag,int gny,int np,int id,int typ,int *ny1,int *ny2,int *ny)
{
int i;
float fslice, fny1, fny2;

if(flag == 0)
   {
   fslice = (float)(gny)/(float)(np) + 3.0;
   fny1 = -2.0;        /* start at -2 so 1st ny2 is iy=fslice+1  */
   }
else
   {
   fslice = (float)(gny + (np-1.0)*4.0)/(float)(np) - 1.0;
   fny1 = 0.0;
   }

for(i=0;i<=id;i++)
   {
   fny2 = fny1 + fslice;

   *ny1 = (int)(fny1 + 0.5);
   *ny2 = (int)(fny2 + 0.5);

   fny1 = fny2 - 3.0;
   }

if(typ == PARL_HEAD)
   *ny1 = 0;

if(typ == PARL_TAIL)
   *ny2 = gny - 1;

*ny2 = *ny2 + 1;      /* SN: in 0 convention */
*ny = *ny2 - *ny1;    /* ny is now local number of planes */
}

double geocenX(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2X(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2kmX(arg,latkm,lonkm,rc,g2)
float *arg, *latkm, *lonkm, *rc, *g2;
{
float cosA, sinA, g2s2, den;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((1.0)/((1.0) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((1.0) + g2s2*((2.0) + (*g2))))*den*den*den;
}

void resample(float *s,int nt,float *dt,int isamp,int ntpad,int ntrsmp,float *newdt,float *p,int ord,float *perc, float *tp)
{
float df, f, f0, fl, fl2, fac;
int i, j;

float one = 1.0;

int minus = -1;
int plus = 1;

taper_norm(s,dt,nt,tp);
zero(s+nt,(ntpad)-(nt));

for(i=ntpad-1;i>=0;i--)
   {
   s[2*i] = s[i];
   s[2*i + 1] = 0.0;
   }

fourg_(s,&ntpad,&minus,p);

if(isamp > 0)
   zero(s+ntpad,2*ntrsmp-ntpad);
else if(isamp < 0)
   {
   if(ord)  /* lowpass at 100*(*perc) % of new Nyquist */
      {
      f0 = (*perc)/(2.0*(*newdt));
      df = 1.0/(ntrsmp*(*newdt));
      for(i=1;i<ntrsmp/2;i++)
         {
         f = i*df;

         fl = f/f0;
         fl2 = fl*fl;
         fl = fl2;
         for(j=1;j<ord;j++)
            fl = fl*fl2;

         fac = one/(one + fl);

         s[2*i] = fac*s[2*i];
         s[2*i + 1] = fac*s[2*i + 1];
         }
      }

   s[ntrsmp] = s[ntrsmp+1] = 0.0; /* zero nyquist */
   }

for(i=1;i<ntrsmp/2;i++)  /* replicate with complex-conjugate */
   {
   s[2*(ntrsmp-i)] = s[2*i];
   s[2*(ntrsmp-i) + 1] = -s[2*i + 1];
   }

fourg_(s,&ntrsmp,&plus,p);

for(i=0;i<ntrsmp;i++)
   s[i] = s[2*i];

norm(s,newdt,ntrsmp);
}

void taper_norm(float *g,float *dt,int nt,float *tp)
{
float fac, df, arg;
int i;
int ntap;

ntap = nt*(*tp);

for(i=0;i<nt-ntap;i++)
   g[i] = g[i]*(*dt);

df = 3.14159/(float)(ntap);
for(i=nt-ntap;i<nt;i++)
   {
   arg = (i-(nt-(ntap+1)))*df;
   fac = (*dt)*0.5*(1.0 + cos(arg));
   g[i] = g[i]*fac;
   }
}

void norm(float *g,float *dt,int nt)
{
float fac;

fac = 1.0/((*dt)*nt);
while(nt--)
   {
   g[0] = g[0]*fac;
   g++;
   }
}

double nt_tol(float fnt,int gnt)
{
double diff;

diff = ((double)(fnt) - (double)(gnt));
if(diff < 0.0)
   diff = -diff;

/*
fprintf(stderr,"diff= %15.10e\n",diff);
*/

return(diff);
}

void check_nan(float *p,float *m,int nx,int nz,int iy,int it,int j)
{
int ip, ix, iz, k;

for(k=0;k<9;k++)
   {
   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
         ip = ix + iz*nx + k*nx*nz;

         if(isnan(p[ip]))
            {
            fprintf(stderr,"%13.5e ix= %d iy= %d iz= %d it= %d field= %d call= %d\n",p[ip],ix,iy,iz,it,k,j);
            fprintf(stderr,"l2m= %13.5e\n",m[ix+iz*nx]);
            fprintf(stderr,"lam= %13.5e\n",m[ix+iz*nx + 1*nx*nz]);
            fprintf(stderr,"mxy= %13.5e\n",m[ix+iz*nx + 2*nx*nz]);
            fprintf(stderr,"mxz= %13.5e\n",m[ix+iz*nx + 3*nx*nz]);
            fprintf(stderr,"myz= %13.5e\n",m[ix+iz*nx + 4*nx*nz]);
            fprintf(stderr,"bxx= %13.5e\n",m[ix+iz*nx + 5*nx*nz]);
            fprintf(stderr,"byy= %13.5e\n",m[ix+iz*nx + 6*nx*nz]);
            fprintf(stderr,"bzz= %13.5e\n",m[ix+iz*nx + 7*nx*nz]);
            fprintf(stderr,"akx= %13.5e\n",m[ix+iz*nx + 8*nx*nz]);
            fprintf(stderr,"pkx= %13.5e\n",m[ix+iz*nx + 9*nx*nz]);
            fprintf(stderr,"skx= %13.5e\n",m[ix+iz*nx +10*nx*nz]);
            fprintf(stderr,"lrw= %13.5e\n",m[ix+iz*nx +11*nx*nz]);
            fprintf(stderr,"mrw= %13.5e\n",m[ix+iz*nx +12*nx*nz]);

            fflush(stderr);
            mpi_exit(-1);
            }
         }
      }
   }
}

void copy_nodeinfo(struct nodeinfo *ni2,struct nodeinfo *ni1)
{
ni2->nproc = ni1->nproc;
ni2->nproc_x = ni1->nproc_x;
ni2->nproc_y = ni1->nproc_y;
ni2->nproc_z = ni1->nproc_z;
ni2->min_nproc = ni1->min_nproc;
ni2->segmentId = ni1->segmentId;
ni2->minusId_x = ni1->minusId_x;
ni2->plusId_x = ni1->plusId_x;
ni2->minusId_y = ni1->minusId_y;
ni2->plusId_y = ni1->plusId_y;
ni2->minusId_z = ni1->minusId_z;
ni2->plusId_z = ni1->plusId_z;
ni2->procId_x = ni1->procId_x;
ni2->procId_y = ni1->procId_y;
ni2->procId_z = ni1->procId_z;
ni2->globnx = ni1->globnx;
ni2->nx1 = ni1->nx1;
ni2->nx2 = ni1->nx2;
ni2->ixminus = ni1->ixminus;
ni2->ixplus = ni1->ixplus;
ni2->loc_nx = ni1->loc_nx;
ni2->globny = ni1->globny;
ni2->ny1 = ni1->ny1;
ni2->ny2 = ni1->ny2;
ni2->iyminus = ni1->iyminus;
ni2->iyplus = ni1->iyplus;
ni2->loc_ny = ni1->loc_ny;
ni2->globnz = ni1->globnz;
ni2->nz1 = ni1->nz1;
ni2->nz2 = ni1->nz2;
ni2->izminus = ni1->izminus;
ni2->izplus = ni1->izplus;
ni2->loc_nz = ni1->loc_nz;
strcpy(ni2->procname,ni1->procname);
}
