#include "include.h"
#include "structure.h"
#include "function.h"
#include <stdint.h>

#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         BLOCK_SIZE      1000000
#define         BNDPAD          3

void get_filepar(char *str,char *file,int *nh,int *latfirst);

main(int ac,char **av)
{
FILE *fopfile(), *fp, *fpr;
float h, *rv, *dv;
float y2, r;
int *rinc, *dinc, nrad, ndep;
int ix, iy, iz, i, id, ir, ifound;
int nx, ny, nz, np, npz;
int xsrc, ysrc;
int *ixv, *iyv, *izv, *izlevel;
double invh, invh2;

int bcnt = 1;
int sflag = 1;
long long *indx, ll_int;

int nhead = 0;
int latfirst = 0;

int ixmin = -999;
int ixmax = -999;
int iymin = -999;
int iymax = -999;
int izstart = 2;
int izmax = -999;

char radiusfile[512], outfile[512], str[1024], nedfile[128];
float xp, yp, north, east, dep, modellon, modellat, modelrot;
float cosR, sinR;

char faultlist[512], infile[1024];
int ixp, iyp, izp;
float flon, flat, fdep, xx, yy, zz;

int geoproj = 1;
int xy2ll = 0;
int ll2xy = 1;
float kmlon, kmlat;

double g0, b0;
double amat[9], ainv[9];

double rperd = RPERD;
float ref_rad = ERAD;

radiusfile[0] = '\0';
faultlist[0] = '\0';
nedfile[0] = '\0';

setpar(ac, av);

mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);
mstpar("h","f",&h);
mstpar("xsrc","d",&xsrc);
mstpar("ysrc","d",&ysrc);

mstpar("modellon","f",&modellon);
mstpar("modellat","f",&modellat);
mstpar("modelrot","f",&modelrot);

mstpar("outfile","s",outfile);

getpar("radiusfile","s",radiusfile);
getpar("faultlist","s",faultlist);
getpar("nedfile","s",nedfile);

getpar("geoproj","d",&geoproj);

getpar("ixmin","d",&ixmin);
getpar("ixmax","d",&ixmax);
getpar("iymin","d",&iymin);
getpar("iymax","d",&iymax);
getpar("izstart","d",&izstart);
getpar("izmax","d",&izmax);

endpar();

gen_matrices(amat,ainv,&modelrot,&modellon,&modellat);

g0 = (double)(0.5*ny*h)/(double)(ref_rad);
b0 = (double)(0.5*nx*h)/(double)(ref_rad);

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

invh = 1.0/(double)(h);
invh2 = invh*invh;

ixv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));
iyv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));
izv = (int *) check_malloc (bcnt*BLOCK_SIZE*sizeof(int));

if(radiusfile[0] != '\0')
   {
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

   izlevel = (int *) check_malloc (nz*sizeof(int));
   npz = 0;
   for(iz=izstart;iz<izmax+izstart;iz++)
      {
      id = 0;
      while(dv[id] < (float)(iz-1.0) && id != ndep-1)
         id++;

      if((iz-izstart)%dinc[id] == 0)
         {
         izlevel[npz] = iz;
         npz++;
         }
      }

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
	    for(iz=0;iz<npz;iz++)
	       {
	       np++;

               if(np > bcnt*BLOCK_SIZE)
                  {
                  bcnt++;

                  ixv = (int *) check_realloc (ixv,bcnt*BLOCK_SIZE*sizeof(int));
                  iyv = (int *) check_realloc (iyv,bcnt*BLOCK_SIZE*sizeof(int));
                  izv = (int *) check_realloc (izv,bcnt*BLOCK_SIZE*sizeof(int));
                  }

	       ixv[np-1] = ix;
	       iyv[np-1] = iy;
	       izv[np-1] = izlevel[iz];
	       }
	    }
         }
      }
   }

/*  now read in fault coords */

if(faultlist[0] != '\0')
   {
   fp = fopfile(faultlist,"r");

   while(fgets(str,1024,fp) != NULL)
      {
      get_filepar(str,infile,&nhead,&latfirst);

      fprintf(stderr,"%s nheader=%d latfirst=%d\n",infile,nhead,latfirst);
      fpr = fopfile(infile,"r");

      i = 0;
      while(i < nhead)
         {
         fgets(str,1024,fpr);
	 i++;
	 }

/*
      i = 0;
      while(str[i] != '\n' && i < 1023)
         i++;
      str[i] = '\0';
      strcpy(infile,str);

      fprintf(stderr,"%s\n",infile);
      fpr = fopfile(infile,"r");

      fgets(str,1024,fpr);
      fgets(str,1024,fpr);
      fgets(str,1024,fpr);
      fgets(str,1024,fpr);
      fgets(str,1024,fpr);
      fgets(str,1024,fpr);
*/

      while(fgets(str,1024,fpr) != NULL)
         {
	 if(latfirst)
            sscanf(str,"%f %f %f",&flat,&flon,&fdep);
	 else
            sscanf(str,"%f %f %f",&flon,&flat,&fdep);

         gcproj(&xx,&yy,&flon,&flat,&ref_rad,&g0,&b0,amat,ainv,ll2xy);

	 /*  go to closest point  */
         ixp = (int)((double)(xx)*invh + 0.5);
         iyp = (int)((double)(yy)*invh + 0.5);
         izp = (int)((double)(fdep)*invh + 1.5);

	 if(ixp >= BNDPAD && ixp < nx-BNDPAD && iyp >= BNDPAD && iyp < ny-BNDPAD && izp >= 1 && izp < nz-BNDPAD)
	    {
	    ifound = 0;
            for(ir=0;ir<np;ir++)
	       {
	       if(ixv[ir] == ixp && iyv[ir] == iyp && izv[ir] == izp)
	          {
	          ifound = 1;
	          break;
	          }
	       }

	    if(ifound == 0)
	       {
	       np++;

               if(np > bcnt*BLOCK_SIZE)
                  {
                  bcnt++;

                  ixv = (int *) check_realloc (ixv,bcnt*BLOCK_SIZE*sizeof(int));
                  iyv = (int *) check_realloc (iyv,bcnt*BLOCK_SIZE*sizeof(int));
                  izv = (int *) check_realloc (izv,bcnt*BLOCK_SIZE*sizeof(int));
                  }

	       ixv[np-1] = ixp;
	       iyv[np-1] = iyp;
	       izv[np-1] = izp;
	       }
	    }
         }

      fclose(fpr);
      }

   fclose(fp);
   }

fprintf(stderr,"np= %d\n",np);

ixv = (int *) check_realloc (ixv,np*sizeof(int));
iyv = (int *) check_realloc (iyv,np*sizeof(int));
izv = (int *) check_realloc (izv,np*sizeof(int));
indx = (long long *) check_malloc (np*sizeof(long long));

for(ir=0;ir<np;ir++)
   indx[ir] = (long long)(ixv[ir])*(long long)(100000000) + (long long)(iyv[ir])*(long long)(10000) + (long long)(izv[ir]);

/* (dumb) sort by increasing indx
*/
sflag = 1;
while(sflag)
   {
   sflag = 0;
   for(ir=0;ir<np-1;ir++)
      {
      if(indx[ir] > indx[ir+1])
         {
	 ix        = ixv[ir+1]; iy        = iyv[ir+1]; iz        = izv[ir+1];
	 ixv[ir+1] = ixv[ir];   iyv[ir+1] = iyv[ir];   izv[ir+1] = izv[ir];
	 ixv[ir]   = ix;        iyv[ir]   = iy;        izv[ir]   = iz;

	 ll_int = indx[ir+1]; indx[ir+1] = indx[ir]; indx[ir] = ll_int;
	 sflag = 1;
	 }
      }
   }

if(nedfile[0] != '\0')
   {
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
   }

if(geoproj == 0)
   sprintf(str,"localized flat earth on spherical ellipsoid (kmlon=%f, kmlat=%f)",kmlon,kmlat);
else if(geoproj == 1)
   sprintf(str,"great circle arcs");

fp = fopfile(outfile,"w");

fprintf(fp,"# geoproj= %d (%s)\n",geoproj,str);
fprintf(fp,"# modellon= %12.5f modellat= %12.5f modelrot= %6.2f\n",modellon,modellat,modelrot);
fprintf(fp,"# xlen= %12.5f ylen= %12.5f\n",nx*h,ny*h);
fprintf(fp,"#\n");

fprintf(fp,"%d\n",np);
for(ir=0;ir<np;ir++)
   {
   xx = ixv[ir]*h;
   yy = iyv[ir]*h;
   zz = (izv[ir]-1)*h;   /* subtract 1 for free-surface */
   gcproj(&xx,&yy,&flon,&flat,&ref_rad,&g0,&b0,amat,ainv,xy2ll);
   fprintf(fp,"%6d %6d %6d %.12lld %11.5f %11.5f %12.5f\n",ixv[ir],iyv[ir],izv[ir],indx[ir],flon,flat,zz);
   }

fclose(fp);
}

void get_filepar(char *str,char *file,int *nh,int *latfirst)
{
int i;
char *sptr;

int lend = 0;

sscanf(str,"%s",file);
*nh = 0;
*latfirst = 0;

sptr = str;

/* 1st skip past 'file' */
while(sptr[0] != ' ' && sptr[0] != '\t')
   {
   if(sptr[0] == '\n' || sptr[0] == '\0')
      return;
   else
      sptr++;
   }

/* now, skip past white space */
while(sptr[0] == ' ' || sptr[0] == '\t')
   {
   if(sptr[0] == '\n' || sptr[0] == '\0')
      return;
   else
      sptr++;
   }

/* now, get 'nheaders' and/or 'latfirst', if they are given */
while(lend == 0)
   {
   i = 0;
   while(sptr[i] != '=')
      {
      if(sptr[i] == '\n' || sptr[i] == '\0')
         return;
      else
         i++;
      }

   if(strncmp(sptr,"nheader",7) == 0)
      sscanf(sptr+i+1,"%d",nh);

   if(strncmp(sptr,"latfirst",8) == 0)
      sscanf(sptr+i+1,"%d",latfirst);

   while(sptr[0] != ' ' && sptr[0] != '\t')
      {
      if(sptr[0] == '\n' || sptr[0] == '\0')
         {
         lend = 1;
         break;
         }
      else
         sptr++;
      }

   while(sptr[0] == ' ' || sptr[0] == '\t')
      {
      if(sptr[0] == '\n' || sptr[0] == '\0')
         {
         lend = 1;
         break;
         }
      else
         sptr++;
      }
   }
}
