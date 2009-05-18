#include <features.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/times.h>

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <sys/file.h>
#include <sys/resource.h>
#include <sys/signal.h>
#include <sys/syscall.h>

#define         FLAT_CONST      298.256
#define         ERAD            6378.139
#define         RPERD           0.017453292
#define         FONE            (float)(1.0)
#define         FTWO            (float)(2.0)

main(ac,av)
int ac;
char **av;
{
float mlat, mlon, lon, lat;
float kperd_n, kperd_e, xs, ys;
float cosR, sinR, xr, yr;

float rotate = 0.0;
float h = 1.0;

float rperd = RPERD;
float erad = ERAD;
float fc = FLAT_CONST;
float g2, radc, latavg;
double geocen();

double xr0, yr0, dlon, dlat, dxr, dyr;
int iway = 0;
int utm_zone = 11;

printf("Enter model origin (long. first, then lat.) and model rotation:\n");
scanf("%f %f %f",&mlon,&mlat,&rotate);

cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

radc = ERAD*RPERD;
set_g2(&g2,&fc);

latavg = geocen(mlat*rperd);
latlon2km(&latavg,&kperd_n,&kperd_e,&radc,&g2);

printf("ke=%12.4f kn=%12.4f latavg=%10.4f\n\n",kperd_e,kperd_n,latavg/rperd);

dlon = mlon;
dlat = mlat;
geoutm_(&dlon,&dlat,&xr0,&yr0,&utm_zone,&iway);

printf("Enter longitude and latitude:\n\n");
while(scanf("%f %f",&lon,&lat)==2)
   {
   xs = (lon - mlon)*kperd_e;
   ys = (mlat - lat)*kperd_n;

   xr = xs*cosR + ys*sinR;
   yr = -xs*sinR + ys*cosR;

   printf("%10.5f %10.5f %10.5f %10.5f\t",lon,lat,xr,yr);

   dlon = lon;
   dlat = lat;
   geoutm_(&dlon,&dlat,&dxr,&dyr,&utm_zone,&iway);

   printf("UTM: %10.5f %10.5f\n",0.001*(dxr-xr0),0.001*(dyr-yr0));
   }
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
