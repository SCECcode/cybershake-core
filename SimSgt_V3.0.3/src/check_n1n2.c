#include "include.h"

int main(int ac,char **av)
{
struct nodeinfo ninfo;
float ftp, fnx, fny, fnz, inv_fmp;
int nx, globnx, nx1, nx2, min_nx, max_nx;
int ny, globny, ny1, ny2, min_ny, max_ny;
int nz, globnz, nz1, nz2, min_nz, max_nz;

int i, target_nproc;
size_t mlen;
int nbuf;
int bflag = 1;

ninfo.nproc_x = -1;
ninfo.nproc_y = -1;
ninfo.nproc_z = -1;
ninfo.min_nproc = 1;

setpar(ac,av);

mstpar("target_nproc","d",&target_nproc);
mstpar("ny","d",&ninfo.globny);
mstpar("nx","d",&ninfo.globnx);
mstpar("nz","d",&ninfo.globnz);

getpar("min_nproc","d",&ninfo.min_nproc);
getpar("nproc_x","d",&ninfo.nproc_x);
getpar("nproc_y","d",&ninfo.nproc_y);
getpar("nproc_z","d",&ninfo.nproc_z);

endpar();

if(ninfo.min_nproc < 1)
   ninfo.min_nproc = 1;
inv_fmp = 1.0/(float)(ninfo.min_nproc);

min_nx = ninfo.globnx;
max_nx = 0;
min_ny = ninfo.globny;
max_ny = 0;
min_nz = ninfo.globnz;
max_nz = 0;

/*
   blfag=0 -> balance only on interior nodes
   blfag=1 -> balance across all nodes

   generally bflag=1 is most efficient
*/

bflag = 1;

ninfo.nproc = target_nproc;
ninfo.nproc = get_nproc(&ninfo);
get_nodetype_neighbor(&ninfo);
get_n1n2(bflag,&ninfo);

fprintf(stderr,"Target nproc= %d\n",target_nproc);
fprintf(stderr,"Actual nproc= %d (nproc_x= %d nproc_y= %d nproc_z= %d)\n",ninfo.nproc,ninfo.nproc_x,ninfo.nproc_y,ninfo.nproc_z);

for(i=0;i<ninfo.nproc_x;i++)
   {
   if(ninfo.loc_nx < min_nx)
      min_nx = ninfo.loc_nx;
   if(ninfo.loc_nx > max_nx)
      max_nx = ninfo.loc_nx;
   }

for(i=0;i<ninfo.nproc_y;i++)
   {
   if(ninfo.loc_ny < min_ny)
      min_ny = ninfo.loc_ny;
   if(ninfo.loc_ny > max_ny)
      max_ny = ninfo.loc_ny;
   }

for(i=0;i<ninfo.nproc_z;i++)
   {
   if(ninfo.loc_nz < min_nz)
      min_nz = ninfo.loc_nz;
   if(ninfo.loc_nz > max_nz)
      max_nz = ninfo.loc_nz;
   }

nbuf = 3*min_nx*min_nz;
if(4*(N_WAVE_VARS-3)*min_nx*min_ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*min_nx*min_ny;
if(2*(N_WAVE_VARS-3)*min_nz*min_ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*min_nz*min_ny;
mlen = sizeof(float)*((N_WAVE_VARS + N_MED_VARS)*min_nx*min_ny*min_nz + nbuf + 6*(min_nx+min_nz));

fprintf(stderr,"*** min_nx= %5d min_ny= %5d min_nz= %5d",min_nx,min_ny,min_nz);
fprintf(stderr," (approximate RAM= %.0f Mbyte)\n",mlen*1.0e-06);

nbuf = 3*max_nx*max_nz;
if(4*(N_WAVE_VARS-3)*max_nx*max_ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*max_nx*max_ny;
if(2*(N_WAVE_VARS-3)*max_nz*max_ny > nbuf)
   nbuf = 4*(N_WAVE_VARS-3)*max_nz*max_ny;
mlen = sizeof(float)*((N_WAVE_VARS + N_MED_VARS)*max_nx*max_ny*max_nz + nbuf + 6*(max_nx+max_nz));

fprintf(stderr,"*** max_nx= %5d max_ny= %5d max_nz= %5d",max_nx,max_ny,max_nz);
fprintf(stderr," (approximate RAM= %.0f Mbyte)\n",mlen*1.0e-06);

/*
for(i=0;i<ninfo.nproc_x;i++)
   {
   ninfo.procId_x = i;
   ninfo.minusId_x = 1;
   ninfo.plusId_x = 1;
   if(i == 0)
      ninfo.minusId_x = -1;
   if(i+1 == ninfo.nproc_x)
      ninfo.plusId_x = -1;

   get_n1n2_indv(1,ninfo.globnx,ninfo.nproc_x,ninfo.procId_x,ninfo.minusId_x,ninfo.plusId_x,&nx1,&nx2,&nx);

   if(nx < min_nx)
      min_nx = nx;
   if(nx > max_nx)
      max_nx = nx;
   }

for(i=0;i<ninfo.nproc_y;i++)
   {
   ninfo.procId_y = i;
   ninfo.minusId_y = 1;
   ninfo.plusId_y = 1;
   if(i == 0)
      ninfo.minusId_y = -1;
   if(i+1 == ninfo.nproc_y)
      ninfo.plusId_y = -1;

   get_n1n2_indv(1,ninfo.globny,ninfo.nproc_y,ninfo.procId_y,ninfo.minusId_y,ninfo.plusId_y,&ny1,&ny2,&ny);

   if(ny < min_ny)
      min_ny = ny;
   if(ny > max_ny)
      max_ny = ny;
   }

for(i=0;i<ninfo.nproc_z;i++)
   {
   ninfo.procId_z = i;
   ninfo.minusId_z = 1;
   ninfo.plusId_z = 1;
   if(i == 0)
      ninfo.minusId_z = -1;
   if(i+1 == ninfo.nproc_z)
      ninfo.plusId_z = -1;

   get_n1n2_indv(1,ninfo.globnz,ninfo.nproc_z,ninfo.procId_z,ninfo.minusId_z,ninfo.plusId_z,&nz1,&nz2,&nz);

   if(nz < min_nz)
      min_nz = nz;
   if(nz > max_nz)
      max_nz = nz;
   }

fprintf(stderr,"*** min_nx= %5d min_ny= %5d min_nz= %5d",min_nx,min_ny,min_nz);
fprintf(stderr," (approximate RAM= %.0f Mbyte)\n",112.0*min_nx*min_ny*min_nz*1.0e-06);

fprintf(stderr,"*** max_nx= %5d max_ny= %5d max_nz= %5d",max_nx,max_ny,max_nz);
fprintf(stderr," (approximate RAM= %.0f Mbyte)\n",112.0*max_nx*max_ny*max_nz*1.0e-06);
*/
}
