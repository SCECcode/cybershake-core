#include "include.h"
#include "structure.h"
#include        "function.h"

main(int ac,char **av)
{
FILE *fopfile(), *fp;
float h, *rv, *dv;
float invh, invh2, y2, r;
int *rinc, *dinc, nrad, ndep;
int ix, iy, iz, i, id, ir;
int nx, ny, nz, np, npz;
int xsrc, ysrc;
int *ixv, *iyv, *izv;

int ixmin = -999;
int ixmax = -999;
int iymin = -999;
int iymax = -999;
int izstart = 2;
int izmax = -999;

char radiusfile[512], outfile[512], str[512], nedfile[128];
float xp, yp, north, east, dep, modelrot;
float cosR, sinR;
float rperd = 0.017453293;

char incords[512], outcords[512];
float *lon, *lat;

int cord_flag = 0;

incords[0] = '\0';

setpar(ac, av);
mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("h","f",&h);
mstpar("xsrc","d",&xsrc);
mstpar("ysrc","d",&ysrc);
mstpar("radiusfile","s",radiusfile);
mstpar("outfile","s",outfile);

mstpar("nedfile","s",nedfile);
mstpar("modelrot","f",&modelrot);

getpar("ixmin","d",&ixmin);
getpar("ixmax","d",&ixmax);
getpar("iymin","d",&iymin);
getpar("iymax","d",&iymax);
getpar("izstart","d",&izstart);
getpar("izmax","d",&izmax);

getpar("incords","s",incords);
if(incords[0] != '\0')
   {
   cord_flag = 1;
   mstpar("outcords","s",outcords);
   }

endpar();

if(ixmin < 0)
   ixmin = 0;
if(ixmax < 0 || ixmax > nx)
   ixmax = nx;

if(iymin < 0)
   iymin = 0;
if(iymax < 0 || iymax > ny)
   iymax = ny;

if(izmax < 0 || izmax > nz-izstart)
   izmax = nz - izstart;

invh = 1.0/h;
invh2 = 1.0/(h*h);

fp = fopfile(radiusfile,"r");

fscanf(fp,"%d",&nrad);

rv = (float *) check_malloc (nrad*sizeof(float));
rinc = (int *) check_malloc (nrad*sizeof(int));

for(i=0;i<nrad;i++)
   fscanf(fp,"%f",&rv[i]);
for(i=0;i<nrad;i++)
   fscanf(fp,"%d",&rinc[i]);

fscanf(fp,"%d",&ndep);

dv = (float *) check_malloc (ndep*sizeof(float));
dinc = (int *) check_malloc (ndep*sizeof(int));

for(i=0;i<ndep;i++)
   fscanf(fp,"%f",&dv[i]);
for(i=0;i<ndep;i++)
   fscanf(fp,"%d",&dinc[i]);

fclose(fp);

/* renormalize values to save computations */
for(i=0;i<nrad;i++)
   rv[i] = rv[i]*rv[i]*invh2;
for(i=0;i<ndep;i++)
   dv[i] = dv[i]*invh;

if(cord_flag)
   {
   lon = (float *) check_malloc (nx*ny*sizeof(float));
   lat = (float *) check_malloc (nx*ny*sizeof(float));
   fp = fopfile(incords,"r");
   
   for(i=0;i<nx*ny;i++)
      fscanf(fp,"%f %f %*d %*d",&lon[i],&lat[i]);

   fclose(fp);
   }

ixv = (int *) check_malloc ((ixmax-ixmin)*(iymax-iymin)*sizeof(int));
iyv = (int *) check_malloc ((ixmax-ixmin)*(iymax-iymin)*sizeof(int));

np = 0;
for(iy=iymin;iy<iymax;iy++)
   {
   y2 = (ysrc-iy);
   y2 = y2*y2;

   for(ix=ixmin;ix<ixmax;ix++)
      {
      r = (xsrc-ix)*(xsrc-ix) + y2;
      ir = 0;
      while(rv[ir] < r && ir != nrad-1)
         ir++;

      if(ix%rinc[ir] == 0 && iy%rinc[ir] == 0)
         {
	 ixv[np] = ix;
	 iyv[np] = iy;
	 np++;
	 }
      }
   }

ixv = (int *) check_realloc (ixv,np*sizeof(int));
iyv = (int *) check_realloc (iyv,np*sizeof(int));

izv = (int *) check_malloc (izmax*sizeof(int));

npz = 0;
for(iz=izstart;iz<izmax+izstart;iz++)
   {
   id = 0;
   while(dv[id] < (float)(iz-1.0) && id != ndep-1)
      id++;

   if((iz-izstart)%dinc[id] == 0)
      {
      izv[npz] = iz;
      npz++;
      }
   }

izv = (int *) check_realloc (izv,npz*sizeof(int));

fprintf(stderr,"np= %d\n",np*npz);

fp = fopfile(nedfile,"w");

cosR = cos((90.0+modelrot)*rperd);
sinR = sin((90.0+modelrot)*rperd);

fprintf(fp,"%d\n",np);
for(ir=0;ir<np;ir++)
   {
   xp = (xsrc - ixv[ir])*h;
   yp = (ysrc - iyv[ir])*h;

   north = xp*cosR - yp*sinR;
   east  = xp*sinR + yp*cosR;

   fprintf(fp,"%13.8f %13.8f %13.8f\n",north,east,sqrt(xp*xp + yp*yp));
   }

fprintf(fp,"%d\n",npz);
for(id=0;id<npz;id++)
   fprintf(fp,"%13.8f\n",(izv[id]-1.0)*h);

fclose(fp);

fp = fopfile(outfile,"w");

fprintf(fp,"%d\n",np*npz);
for(id=0;id<npz;id++)
   {
   for(ir=0;ir<np;ir++)
      {
      fprintf(fp,"%6d %6d %6d\n",ixv[ir],iyv[ir],izv[id]);
      }
   }

fclose(fp);

if(cord_flag)
   {
   fp = fopfile(outcords,"w");

   for(ir=0;ir<np;ir++)
      {
      i = ixv[ir]+iyv[ir]*nx;
      fprintf(fp,"%10.4f %10.4f\n",lon[i],lat[i]);
      }

   fclose(fp);
   }
}
