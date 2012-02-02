#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

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

void _zapit(float *s, int n)
{
  while(n--)
    {
      s[0] = 0.0;
      s++;
    }
}

