#include "include.h"
#include "function.h"

int main(int ac, char **av)
{
float *vp, *vs, *dn;
int ip, ix, iy, iz, nx, ny, nz;
int fdp, fds, fdd;
char pmodfile[256], smodfile[256], dmodfile[256];

int ixp = -1;
int iyp = -1;
int izp = -1;

setpar(ac, av);
mstpar("pmodfile","s",pmodfile);
mstpar("smodfile","s",smodfile);
mstpar("dmodfile","s",dmodfile);

mstpar("nx","d",&nx);
mstpar("ny","d",&ny);
mstpar("nz","d",&nz);

getpar("ixp","d",&ixp);
getpar("iyp","d",&iyp);
getpar("izp","d",&izp);
endpar();

vp = (float *) check_malloc (nx*nz*sizeof(float));
vs = (float *) check_malloc (nx*nz*sizeof(float));
dn = (float *) check_malloc (nx*nz*sizeof(float));

fdp = opfile_ro(pmodfile);
fds = opfile_ro(smodfile);
fdd = opfile_ro(dmodfile);
fprintf(stderr,"sizeof(size_t)= %d\n",sizeof(size_t));
fprintf(stderr,"sizeof(off_t)= %d\n",sizeof(off_t));
fprintf(stderr,"current offsets: fdp= %Ld fds= %Ld fdd= %Ld\n",(off_t)(lseek(fdp,(off_t)0,SEEK_CUR)),(off_t)(lseek(fds,(off_t)0,SEEK_CUR)),(off_t)(lseek(fdd,(off_t)0,SEEK_CUR)));

for(iy=0;iy<ny;iy++)
   {
   fprintf(stderr,"iy= %5d\n",iy);
   reed(fdp,vp,nx*nz*sizeof(float));
   reed(fds,vs,nx*nz*sizeof(float));
   reed(fdd,dn,nx*nz*sizeof(float));
   for(iz=0;iz<nz;iz++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 ip = ix + iz*nx;
	 if(ix==ixp && iy==iyp && iz==izp)
	    {
	    fprintf(stderr,"current offsets: fdp= %Ld fds= %Ld fdd= %Ld\n",(off_t)(lseek(fdp,(off_t)0,SEEK_CUR)),(off_t)(lseek(fds,(off_t)0,SEEK_CUR)),(off_t)(lseek(fdd,(off_t)0,SEEK_CUR)));
	    fprintf(stderr,"lam= %13.5e mu= %13.5e vp= %13.5e vs= %13.5e dn= %13.5e ix= %5d iy= %5d iz= %5d\n",(vp[ip]*vp[ip]-2.0*vs[ip]*vs[ip])*dn[ip],vs[ip]*vs[ip]*dn[ip],vp[ip],vs[ip],dn[ip],ix,iy,iz);
	    }

	 if(vp[ip]/vs[ip] <= sqrt(2.0))
	    {
	    fprintf(stderr,"ERROR vp= %13.5e vs= %13.5e dn= %13.5e ix= %5d iy= %5d iz= %5d\n",vp[ip],vs[ip],dn[ip],ix,iy,iz);
	    }
         }
      }
   }
close(fdp);
close(fds);
close(fdd);
}
