#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void *check_malloc(size_t len)
{
void *ptr;

ptr = (void *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

double gaus_rand(float *sigma,float *mean,long *seed)
{
double r = 0.0;
double six = 6.0;
double one = 1.0;
double half = 0.5;
int i;

for(i=0;i<12;i++)
   r = r + half*(one + sfrand(seed));

return((double)((r - six)*(*sigma) + *mean));
}

void zapit(float *s, int n)
{
  while(n--)
    {
      s[0] = 0.0;
      s++;
    }
}

void set_ne(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*sn = ((*slat) - (*elat))*kperd_n;
*se = ((*slon) - (*elon))*kperd_e;
}

void set_ll(float *elon,float *elat,float *slon,float *slat,float *sn,float *se)
{
float kperd_n, kperd_e;
double e2, den, g2, lat0;
float cosA, sinA;

double rperd = 0.017453293;
double radius = 6378.139;
double f = 298.256;

f = 1.0/f;
e2 = 2.0*f - f*f;
g2 = e2/((1.0 - f)*(1.0 - f));

lat0 = atan((1.0 - f)*tan((*elat)*rperd));

cosA = cos(lat0);
sinA = sin(lat0);

den = sqrt(1.0/(1.0 + g2*sinA*sinA));
kperd_e = rperd*radius*cosA*den;
kperd_n = rperd*radius*(sqrt(1.0 + g2*sinA*sinA*(2.0 + g2)))*den*den*den;

*slat = (*sn)/kperd_n + *elat;
*slon = (*se)/kperd_e + *elon;
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
