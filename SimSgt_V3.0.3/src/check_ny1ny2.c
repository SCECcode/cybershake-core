#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define PARL_HEAD 1
#define PARL_BODY 2
#define PARL_TAIL 3

void get_ny1ny2(int,int,int,int,int,int *,int *,int *);

int main(int ac,char **av)
{
int globny, ny, ny1, ny2, nproc, i, nodeType;
int min_ny, max_ny;

int nx = -1;
int nz = -1;

/*
fprintf(stdout,"CLK_TCK= %d  CLOCKS_PER_SEC= %d\n",sysconf(_SC_CLK_TCK),CLOCKS_PER_SEC);
exit(-1);
*/

setpar(ac,av);
mstpar("nproc","d",&nproc);
mstpar("globny","d",&globny);
getpar("nx","d",&nx);
getpar("nz","d",&nz);
endpar();

min_ny = globny;
max_ny = 0;

for(i=0;i<nproc;i++)
   {
   nodeType = PARL_BODY;
   if(i == 0)
      nodeType = PARL_HEAD;
   if(i+1 == nproc)
      nodeType = PARL_TAIL;

   get_ny1ny2(1,globny,nproc,i,nodeType,&ny1,&ny2,&ny);

   /*
   fprintf(stderr,"NEW proc= %5d localny= %5d: ny1= %5d ny2= %5d globny= %5d\n",i,ny,ny1,ny2,globny);

   get_ny1ny2(0,globny,nproc,i,nodeType,&ny1,&ny2,&ny);

   fprintf(stderr,"OLD proc= %5d localny= %5d: ny1= %5d ny2= %5d globny= %5d\n",i,ny,ny1,ny2,globny);
   */

   if(ny < min_ny)
      min_ny = ny;
   if(ny > max_ny)
      max_ny = ny;
   }

fprintf(stderr,"***** min_ny= %5d",min_ny);
if(nx > 0 && nz > 0)
   fprintf(stderr," (approximate RAM= %.0f Mbyte)",112.0*min_ny*nx*nz*1.0e-06);
fprintf(stderr,"\n");

fprintf(stderr,"***** max_ny= %5d",max_ny);
if(nx > 0 && nz > 0)
   fprintf(stderr," (approximate RAM= %.0f Mbyte)",112.0*max_ny*nx*nz*1.0e-06);
fprintf(stderr,"\n");

}

void get_ny1ny2(int flag,int gny,int np,int id,int typ,int *ny1,int *ny2,int *ny)
{
int i;
float fslice, fny1, fny2;

if(flag == 0)
   {
   fslice = (float)(gny)/(float)(np) + 3.0;
   fny1 = -2.0;
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

*ny2 = *ny2 + 1;
*ny = *ny2 - *ny1;
}
