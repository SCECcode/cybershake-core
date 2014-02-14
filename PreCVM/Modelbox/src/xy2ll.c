/*

   Assumes that local (x,y) coordinate system has y positive along strike
   and x positive along 90 deg clockwise rotation from y positive direction.

*/

#include <sys/file.h>
#include <stdio.h>
#include <math.h>

#define		MAX_STAT 10000
#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

void *check_malloc(int);
FILE *fopfile(char*, char*);

main(ac,av)
int ac;
char **av;
{
FILE *fopfile(), *fpr, *fpw;
float lat, lon, strike;
float kperd_n, kperd_e, ar1, ar2, ar3, ar4;
float xp, yp, zp, plon, plat;
float cosT, sinT;
int nr;

char infile[512], outfile[512], str[512];
char name[16];

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

setpar(ac, av);
mstpar("lat","f",&lat);
mstpar("lon","f",&lon);
mstpar("strike","f",&strike);
mstpar("infile","s",infile);
mstpar("outfile","s",outfile);
endpar();

cosT = cos(strike*rperd);
sinT = sin(strike*rperd);

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = geocen(lat*rperd);
latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

ar1 = sinT/kperd_e;
ar2 = cosT/kperd_e;
ar3 = sinT/kperd_n;
ar4 = cosT/kperd_n;

name[0] = '\0';

fpw = fopfile(outfile,"w");
fpr = fopfile(infile,"r");

while(fgets(str,512,fpr) != NULL)
   {
   nr = sscanf(str,"%s %f %f",name,&xp,&yp);

   plon = lon + yp*ar1 + xp*ar2;
   plat = lat + yp*ar4 - xp*ar3;

   fprintf(fpw,"%10.4f %10.4f %s %8.2f %8.2f\n",plon,plat,name,xp,yp);
   }

fclose(fpr);
fclose(fpw);
}

double geocen(x)
double x;
{
double r;
r = atan((1.0 - (1.0/FLAT_CONST))*tan(x));
return(r);
}

set_g2(g2,fc)
float *g2, *fc;
{
float f;

f = (1.0)/(*fc);
*g2 = ((2.0)*f - f*f)/(((1.0) - f)*((1.0) - f));
}

latlon2km(arg,latkm,lonkm,rc,g2)
float *arg, *latkm, *lonkm, *rc, *g2;
{
float cosA, sinA, g2s2, den;

cosA = cos((*arg));
sinA = sin((*arg));
g2s2 = (*g2)*sinA*sinA;

den = sqrt((FONE)/((FONE) + g2s2));
*lonkm = (*rc)*cosA*den;
*latkm = (*rc)*(sqrt((FONE) + g2s2*((FTWO) + (*g2))))*den*den*den;
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

void *check_malloc(int len)
{
char *ptr;

ptr = (char *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}
