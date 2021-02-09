#include <stdio.h>
#include <stdlib.h>

#include "structure.h"
#include "functions.h"

void *check_malloc(size_t len)
{
void *ptr;

//fprintf(stderr, "Allocating %ld bytes.\n", len);
ptr = (void *) malloc (len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   perror("");
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
   perror("");
   exit(-1);
   }

return(ptr);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

double gaussian_rand(float *sigma,float *mean,long *seed)
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

zapit(s,n)
float *s;
int n;
{
while(n--)
   {
   s[0] = 0.0;
   s++;
   }
}

