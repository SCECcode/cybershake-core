#include <sys/file.h>
#include <stdio.h>
#include <math.h>

#define		MAX_STAT 10000

char *check_malloc();

struct statinfo
   {
   char name[8];
   float lat;
   float lon;
   int x;
   int y;
   int z;
   }

main(ac,av)
int ac;
char **av;
{
FILE *fopfile(), *fp;
struct statinfo si[MAX_STAT];
float h, mlat, mlon;
float kperd_e, mlat1, mlon1, xlen, ylen, xs, ys;
int i, ns, nx, ny, test;
float cosR, sinR, xr, yr;

char infile[512], outfile[512], str[512];

float dperr = 57.2957795;
float rperd = 0.017453293;
float kperd_n = 111.19;
float rotate = 0.0;

setpar(ac, av);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("h","f",&h);
getpar("rotate","f",&rotate);
endpar();

xlen = nx*h;
ylen = ny*h;

mlat1 = mlat - ylen/kperd_n;
kperd_e = kperd_n*cos(0.5*(mlat1+mlat)*rperd);
mlon1 = mlon + xlen/kperd_e;

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

fp = fopfile(infile,"r");

ns = 0;
while(fgets(str,512,fp) != NULL)
   {
   sscanf(str,"%f %f %s",&si[ns].lon,&si[ns].lat,si[ns].name);

   test = 1;
   for(i=0;i<ns;i++)
      {
      if(strcmp(si[i].name,si[ns].name) == 0)
	 {
	 test = 0;
	 fprintf(stderr,"Stat '%s' duplicated in input file, only first occurrence will be used\n",si[ns].name);
	 }
      }

   if(test)
      {
      xs = (si[ns].lon - mlon)*kperd_e;
      ys = (mlat - si[ns].lat)*kperd_n;

      xr = xs*cosR + ys*sinR;
      yr = -xs*sinR + ys*cosR;

      si[ns].x = (int)(xr/h + 0.5);
      si[ns].y = (int)(yr/h + 0.5);
      si[ns].z = 1;

      if(si[ns].x >= 0 && si[ns].x < nx && si[ns].y >= 0 && si[ns].y <= ny)
         ns++;
      }
   }
fclose(fp);

fp = fopfile(outfile,"w");

fprintf(fp,"%d\n",ns);
for(i=0;i<ns;i++)
   fprintf(fp,"%5d %5d %5d %s\n",si[i].x,si[i].y,si[i].z,si[i].name);
fclose(fp);
}
