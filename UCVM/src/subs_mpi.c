#include "mpi.h"
#include "include.h"
#include "func_mpi.h"

void mpi_sndrcv2(void *snd1,void *rcv1,void *snd2,void *rcv2,int len,int id,int sf,int rf)
{
MPI_Status s;

MPI_Sendrecv(snd1,len,MPI_CHAR,id,sf,rcv1,len,MPI_CHAR,id,rf,MPI_COMM_WORLD,&s);
MPI_Sendrecv(snd2,len,MPI_CHAR,id,sf+1,rcv2,len,MPI_CHAR,id,rf+1,MPI_COMM_WORLD,&s);
}

void mpi_global_val(void *snd,void *rcv,char *name,int nu,MPI_Datatype typ,MPI_Op op)
{
if(typ == MPI_INT)
   fprintf(stderr,"\t resolving '%s': \tsend= %d \t",name,*((int *)snd));
if(typ == MPI_FLOAT)
   fprintf(stderr,"\t resolving '%s': \tsend= %g \t",name,*((float *)snd));
if(typ == MPI_DOUBLE)
   fprintf(stderr,"\t resolving '%s': \tsend= %lg \t",name,*((double *)snd));

fflush(stderr);

MPI_Reduce(snd,rcv,nu,typ,op,0,MPI_COMM_WORLD);
MPI_Bcast(rcv,nu,typ,0,MPI_COMM_WORLD);

if(typ == MPI_INT)
   fprintf(stderr,"recv= %d\n",*((int *)rcv));
if(typ == MPI_FLOAT)
   fprintf(stderr,"recv= %g\n",*((float *)rcv));
if(typ == MPI_DOUBLE)
   fprintf(stderr,"recv= %lg\n",*((double *)rcv));

fflush(stderr);
}

void mpi_exit(int val)
{
MPI_Finalize();
exit(val);
}

void mpi_barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
 return;

}

void mpi_init(int *ac,char ***av,int *np,int *id,char *pname,int *len)
{
MPI_Init(ac,av);
MPI_Comm_size(MPI_COMM_WORLD,np);
MPI_Comm_rank(MPI_COMM_WORLD,id);

MPI_Get_processor_name(pname,len);
check_procname(pname,len);
}

void mpi_final(char *s)
{
fprintf(stderr,"%s\n",s);
MPI_Finalize();
}

void check_procname(char *name,int *len)
{
int l;

if(name[*len-1] == '.')
   name[*len-1] = '\0';

*len = *len - 1;
}

void get_iylimits(int np,int id,int *iy0,int *iy1,int ny)
{
float fslice, fy0, fy1;
int i;

fslice = (float)(ny)/(float)(np);
fy0 = 0.0;

for(i=0;i<=id;i++)
   {
   fy1 = fy0 + fslice;

   *iy0 = (int)(fy0 + 0.5);
   *iy1 = (int)(fy1 + 0.5);

   fy0 = fy1;
   }
}
