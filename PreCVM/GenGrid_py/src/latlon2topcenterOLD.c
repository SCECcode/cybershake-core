#include <sys/file.h>
#include <stdio.h>
#include <math.h>

main(ac,av)
int ac;
char **av;
{
float elat, elon, mlat, mlon;
float kperd_e, xtop, ytop;
float x0, y0, cosR, sinR;

float dperr = 57.2957795;
float rperd = 0.017453293;
float kperd_n = 111.19;
float rotate = 0.0;

setpar(ac, av);
mstpar("elat","f",&elat);
mstpar("elon","f",&elon);
mstpar("mlat","f",&mlat);
mstpar("mlon","f",&mlon);
getpar("rotate","f",&rotate);
endpar();

kperd_e = kperd_n*cos(0.5*(elat+mlat)*rperd);
cosR = cos(rotate*rperd);
sinR = sin(rotate*rperd);

x0 = (elon - mlon)*kperd_e;
y0 = (mlat - elat)*kperd_n;
xtop = x0*cosR + y0*sinR;
ytop = -x0*sinR + y0*cosR;

printf("%8.1f %8.1f\n",xtop,ytop);
}
