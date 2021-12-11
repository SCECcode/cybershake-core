#include "include.h"

int main(int ac,char **av)
{
struct nodeinfo ninfo;
int nx, globnx, nx1, nx2, min_nx, max_nx;
int ny, globny, ny1, ny2, min_ny, max_ny;
int nz, globnz, nz1, nz2, min_nz, max_nz;

int npret, i, pnlen;

char logfile[256], name[128], logdir[128];

ninfo.nproc_x = -1;
ninfo.nproc_y = -1;
ninfo.nproc_z = -1;
ninfo.min_nproc = 1;

sprintf(name,"get_nodeinfo");
sprintf(logdir,".");

setpar(ac,av);

mstpar("nproc","d",&ninfo.nproc);
mstpar("ny","d",&ny);
mstpar("nx","d",&nx);
mstpar("nz","d",&nz);

getpar("min_nproc","d",&ninfo.min_nproc);
getpar("nproc_x","d",&ninfo.nproc_x);
getpar("nproc_y","d",&ninfo.nproc_y);
getpar("nproc_z","d",&ninfo.nproc_z);

getpar("name","s",name);
getpar("logdir","s",logdir);

endpar();

ninfo.globnx = nx;
ninfo.globny = ny;
ninfo.globnz = nz;

makedir(logdir);

if(ninfo.min_nproc < 1)
   ninfo.min_nproc = 1;

mpi_init(&ac,&av,&ninfo.nproc,&ninfo.segmentId,ninfo.procname,&pnlen);

sprintf(logfile,"%s/%s-%.5d.rlog",logdir,name,ninfo.segmentId);
freopen(logfile,"w",stderr);

npret = get_nproc(&ninfo);
get_nodetype_neighbor(&ninfo);
get_n1n2(1,&ninfo);

fprintf(stderr,"nproc= %d (nproc_x= %d nproc_y= %d nproc_z= %d)\n",ninfo.nproc,ninfo.nproc_x,ninfo.nproc_y,ninfo.nproc_z);
fprintf(stderr,"      Running on node= %s\n",ninfo.procname);
fprintf(stderr,"               nodeId= %5d\n",ninfo.segmentId);
fprintf(stderr,"            minusId_x= %5d     plusId_x= %5d\n",ninfo.minusId_x,ninfo.plusId_x);
fprintf(stderr,"            minusId_y= %5d     plusId_y= %5d\n",ninfo.minusId_y,ninfo.plusId_y);
fprintf(stderr,"            minusId_z= %5d     plusId_z= %5d\n\n",ninfo.minusId_z,ninfo.plusId_z);
fflush(stderr);

fprintf(stderr,"**** Distribution of model subsets:\n");
fprintf(stderr,"     nx1= %5d ixm= %5d ixp= %5d nx2= %5d loc_nx= %5d globnx= %5d\n",ninfo.nx1,ninfo.ixminus,ninfo.ixplus,ninfo.nx2,ninfo.loc_nx,ninfo.globnx);
fprintf(stderr,"     ny1= %5d iym= %5d iyp= %5d ny2= %5d loc_ny= %5d globny= %5d\n",ninfo.ny1,ninfo.iyminus,ninfo.iyplus,ninfo.ny2,ninfo.loc_ny,ninfo.globny);
fprintf(stderr,"     nz1= %5d izm= %5d izp= %5d nz2= %5d loc_nz= %5d globnz= %5d\n",ninfo.nz1,ninfo.izminus,ninfo.izplus,ninfo.nz2,ninfo.loc_nz,ninfo.globnz);
fflush(stderr);

if((ninfo.minusId_x == -1 && ninfo.nx1 != ninfo.ixminus) || (ninfo.minusId_x >= 0 && ninfo.nx1 != ninfo.ixminus - 2))
   fprintf(stderr,"OOPS ixminus problem nx1= %5d ixm= %5d\n",ninfo.nx1,ninfo.ixminus);

if((ninfo.plusId_x == -1 && ninfo.nx2 != ninfo.ixplus+1) || (ninfo.plusId_x >= 0 && ninfo.nx2 != ninfo.ixplus + 3))
   fprintf(stderr,"OOPS ixplus problem ixp= %5d nx2= %5d\n",ninfo.ixplus,ninfo.nx2);

if((ninfo.minusId_y == -1 && ninfo.ny1 != ninfo.iyminus) || (ninfo.minusId_y >= 0 && ninfo.ny1 != ninfo.iyminus - 2))
   fprintf(stderr,"OOPS iyminus problem ny1= %5d iym= %5d\n",ninfo.ny1,ninfo.iyminus);

if((ninfo.plusId_y == -1 && ninfo.ny2 != ninfo.iyplus+1) || (ninfo.plusId_y >= 0 && ninfo.ny2 != ninfo.iyplus + 3))
   fprintf(stderr,"OOPS iyplus problem iyp= %5d ny2= %5d\n",ninfo.iyplus,ninfo.ny2);

if((ninfo.minusId_z == -1 && ninfo.nz1 != ninfo.izminus) || (ninfo.minusId_z >= 0 && ninfo.nz1 != ninfo.izminus - 2))
   fprintf(stderr,"OOPS izminus problem nz1= %5d izm= %5d\n",ninfo.nz1,ninfo.izminus);

if((ninfo.plusId_z == -1 && ninfo.nz2 != ninfo.izplus+1) || (ninfo.plusId_z >= 0 && ninfo.nz2 != ninfo.izplus + 3))
   fprintf(stderr,"OOPS izplus problem izp= %5d nz2= %5d\n",ninfo.izplus,ninfo.nz2);

mpi_final("PROGRAM get_nodeinfo IS FINISHED");
exit(0);
}
