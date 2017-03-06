#include "include.h"

void mpi_sndrcv2(void *snd1,void *rcv1,void *snd2,void *rcv2,int len,int id,int sf,int rf)
{
MPI_Status s;
char *sp1, *sp2, *rp1, *rp2;
int rlen;
int limit_len = 10000000; /* RWG 11/16/2009: 10 Mb seems to be OK? */

rlen = len;
sp1 = snd1;
rp1 = rcv1;
sp2 = snd2;
rp2 = rcv2;
while(rlen >= limit_len)
   {
   MPI_Sendrecv(sp1,limit_len,MPI_CHAR,id,sf,rp1,limit_len,MPI_CHAR,id,rf,MPI_COMM_WORLD,&s);
   MPI_Sendrecv(sp2,limit_len,MPI_CHAR,id,sf+1,rp2,limit_len,MPI_CHAR,id,rf+1,MPI_COMM_WORLD,&s);

   rlen = rlen - limit_len;
   sp1 = sp1 + limit_len;
   rp1 = rp1 + limit_len;
   sp2 = sp2 + limit_len;
   rp2 = rp2 + limit_len;
   }

if(rlen > 0)
   {
   MPI_Sendrecv(sp1,rlen,MPI_CHAR,id,sf,rp1,rlen,MPI_CHAR,id,rf,MPI_COMM_WORLD,&s);
   MPI_Sendrecv(sp2,rlen,MPI_CHAR,id,sf+1,rp2,rlen,MPI_CHAR,id,rf+1,MPI_COMM_WORLD,&s);
   }
}

void mpi_global_val(void *snd,void *rcv,char *name,int nu,MPI_Datatype typ,MPI_Op op)
{

/*
switch(typ)
   {
   case MPI_INT:
      fprintf(stderr,"\t resolving '%s': \tsend= %d \t",name,*((int *)snd));
      break;
   case MPI_FLOAT:
      fprintf(stderr,"\t resolving '%s': \tsend= %g \t",name,*((float *)snd));
      break;
   case MPI_DOUBLE:
      fprintf(stderr,"\t resolving '%s': \tsend= %lg \t",name,*((double *)snd));
      break;
   }
*/

if(typ == MPI_INT)
   fprintf(stderr,"\t resolving '%s': \tsend= %d \t",name,*((int *)snd));
if(typ == MPI_FLOAT)
   fprintf(stderr,"\t resolving '%s': \tsend= %g \t",name,*((float *)snd));
if(typ == MPI_DOUBLE)
   fprintf(stderr,"\t resolving '%s': \tsend= %lg \t",name,*((double *)snd));

fflush(stderr);

MPI_Reduce(snd,rcv,nu,typ,op,0,MPI_COMM_WORLD);
MPI_Bcast(rcv,nu,typ,0,MPI_COMM_WORLD);

/*
switch(typ)
   {
   case MPI_INT:
      fprintf(stderr,"recv= %d\n",*((int *)rcv));
      break;
   case MPI_FLOAT:
      fprintf(stderr,"recv= %g\n",*((float *)rcv));
      break;
   case MPI_DOUBLE:
      fprintf(stderr,"recv= %lg\n",*((double *)rcv));
      break;
   }
*/

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

size_t exchange_pvP3(float *pvf,float *tbuf,struct nodeinfo *ni,int it,int iflag)
{
size_t tot_blen;
float *pvptr, *tbptr, *snd1, *snd2, *rcv1, *rcv2;;
int ix, iy, iz, iv, nx, ny, nz, nv, ib, nb, noff;
int blen, id, sndtag, rcvtag;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

if(iflag == VEL)
   {
   noff = 0;
   nv = 3;
   }
else if(iflag == TAU)
   {
   noff = 3*nx*nz;
   nv = N_WAVE_VARS - 3;
   }

sndtag = it + ni->segmentId;
tot_blen = 0;

if(ni->minusId_x >= 0)
   {
   id = ni->minusId_x;
   rcvtag = it + id;

   tbptr = tbuf;
   for(ix=2;ix<4;ix++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       tbptr[ib] = pvptr[iz*nx];

	       ib++;
	       }
            }
         }

      nb = ib;
      tbptr = tbptr + nb;
      }

   blen = nb*sizeof(float);
   snd1 = tbuf + 0*nb;
   rcv1 = tbuf + 2*nb;
   snd2 = tbuf + 1*nb;
   rcv2 = tbuf + 3*nb;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tbptr = tbuf + 2*nb;
   for(ix=0;ix<2;ix++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       pvptr[iz*nx] = tbptr[ib];

	       ib++;
	       }
            }
         }

      tbptr = tbptr + nb;
      }

   tot_blen = tot_blen + 4*blen;
   }

if(ni->plusId_x >= 0)
   {
   id = ni->plusId_x;
   rcvtag = it + id;

   tbptr = tbuf;
   for(ix=nx-4;ix<nx-2;ix++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       tbptr[ib] = pvptr[iz*nx];

	       ib++;
	       }
            }
         }

      nb = ib;
      tbptr = tbptr + nb;
      }

   blen = nb*sizeof(float);
   snd1 = tbuf + 0*nb;
   rcv1 = tbuf + 2*nb;
   snd2 = tbuf + 1*nb;
   rcv2 = tbuf + 3*nb;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tbptr = tbuf + 2*nb;
   for(ix=nx-2;ix<nx;ix++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       pvptr[iz*nx] = tbptr[ib];

	       ib++;
	       }
            }
         }

      tbptr = tbptr + nb;
      }

   tot_blen = tot_blen + 4*blen;
   }

if(ni->minusId_y >= 0)
   {
   id = ni->minusId_y;
   rcvtag = it + id;

   blen = nv*nx*nz*sizeof(float);
   snd1 = pvf + 2*N_WAVE_VARS*nx*nz + noff;
   rcv1 = pvf + 0*N_WAVE_VARS*nx*nz + noff;
   snd2 = pvf + 3*N_WAVE_VARS*nx*nz + noff;
   rcv2 = pvf + 1*N_WAVE_VARS*nx*nz + noff;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tot_blen = tot_blen + 4*blen;
   }

if(ni->plusId_y >= 0)
   {
   id = ni->plusId_y;
   rcvtag = it + id;

   blen = nv*nx*nz*sizeof(float);
   snd1 = pvf + (ny-4)*N_WAVE_VARS*nx*nz + noff;
   rcv1 = pvf + (ny-2)*N_WAVE_VARS*nx*nz + noff;
   snd2 = pvf + (ny-3)*N_WAVE_VARS*nx*nz + noff;
   rcv2 = pvf + (ny-1)*N_WAVE_VARS*nx*nz + noff;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tot_blen = tot_blen + 4*blen;
   }

if(ni->minusId_z >= 0)
   {
   id = ni->minusId_z;
   rcvtag = it + id;

   tbptr = tbuf;
   for(iz=2;iz<4;iz++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       tbptr[ib] = pvptr[ix + iz*nx + iv*nx*nz];

	       ib++;
	       }
            }
         }

      nb = ib;
      tbptr = tbptr + nb;
      }

   blen = nb*sizeof(float);
   snd1 = tbuf + 0*nb;
   rcv1 = tbuf + 2*nb;
   snd2 = tbuf + 1*nb;
   rcv2 = tbuf + 3*nb;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tbptr = tbuf + 2*nb;
   for(iz=0;iz<2;iz++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       pvptr[ix + iz*nx + iv*nx*nz] = tbptr[ib];

	       ib++;
	       }
            }
         }

      tbptr = tbptr + nb;
      }

   tot_blen = tot_blen + 4*blen;
   }

if(ni->plusId_z >= 0)
   {
   id = ni->plusId_z;
   rcvtag = it + id;

   tbptr = tbuf;
   for(iz=nz-4;iz<nz-2;iz++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       tbptr[ib] = pvptr[ix + iz*nx + iv*nx*nz];

	       ib++;
	       }
            }
         }

      nb = ib;
      tbptr = tbptr + nb;
      }

   blen = nb*sizeof(float);
   snd1 = tbuf + 0*nb;
   rcv1 = tbuf + 2*nb;
   snd2 = tbuf + 1*nb;
   rcv2 = tbuf + 3*nb;

   mpi_sndrcv2(snd1,rcv1,snd2,rcv2,blen,id,sndtag,rcvtag);

   tbptr = tbuf + 2*nb;
   for(iz=nz-2;iz<nz;iz++)
      {
      ib = 0;
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       pvptr[ix + iz*nx + iv*nx*nz] = tbptr[ib];

	       ib++;
	       }
            }
         }

      tbptr = tbptr + nb;
      }

   tot_blen = tot_blen + 4*blen;
   }

return(tot_blen);
}

size_t exchange1_pvP3(float *pvf,float *tbuf,struct nodeinfo *ni,int it,int iflag)
{
MPI_Status s;
size_t tot_blen;
float *pvptr, *snd1, *rcv1;
int ix, iy, iz, iv, nx, ny, nz, nv, ib, nb, noff;
int blen, id, sndtag, rcvtag;

nx = ni->loc_nx;
ny = ni->loc_ny;
nz = ni->loc_nz;

if(iflag == VEL)
   {
   noff = 0;
   nv = 3;
   }
else if(iflag == TAU)
   {
   noff = 3*nx*nz;
   nv = N_WAVE_VARS - 3;
   }

sndtag = it + ni->segmentId;
tot_blen = 0;

if(ni->minusId_x >= 0)
   {
   id = ni->minusId_x;
   rcvtag = it + id;

   nb = 2*nv*ny*nz;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(ix=2;ix<4;ix++)
      {
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       snd1[ib] = pvptr[iz*nx];

	       ib++;
	       }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(ix=0;ix<2;ix++)
      {
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       pvptr[iz*nx] = rcv1[ib];

	       ib++;
	       }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

if(ni->plusId_x >= 0)
   {
   id = ni->plusId_x;
   rcvtag = it + id;

   nb = 2*nv*ny*nz;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(ix=nx-4;ix<nx-2;ix++)
      {
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       snd1[ib] = pvptr[iz*nx];

	       ib++;
	       }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(ix=nx-2;ix<nx;ix++)
      {
      for(iy=0;iy<ny;iy++)
         {
         for(iv=0;iv<nv;iv++)
            {
            pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix + iv*nx*nz;

            for(iz=0;iz<nz;iz++)
               {
	       pvptr[iz*nx] = rcv1[ib];

	       ib++;
	       }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

if(ni->minusId_y >= 0)
   {
   id = ni->minusId_y;
   rcvtag = it + id;

   nb = 2*nv*nx*nz;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(iy=2;iy<4;iy++)
      {
      pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

      for(iv=0;iv<nv;iv++)
         {
         for(iz=0;iz<nz;iz++)
            {
            for(ix=0;ix<nx;ix++)
               {
               snd1[ib] = pvptr[ix+iz*nx + iv*nx*nz];

               ib++;
               }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(iy=0;iy<2;iy++)
      {
      pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

      for(iv=0;iv<nv;iv++)
         {
         for(iz=0;iz<nz;iz++)
            {
            for(ix=0;ix<nx;ix++)
               {
               pvptr[ix + iz*nx + iv*nx*nz] = rcv1[ib];

               ib++;
               }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

if(ni->plusId_y >= 0)
   {
   id = ni->plusId_y;
   rcvtag = it + id;

   nb = 2*nv*nx*nz;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(iy=ny-4;iy<ny-2;iy++)
      {
      pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

      for(iv=0;iv<nv;iv++)
         {
         for(iz=0;iz<nz;iz++)
            {
            for(ix=0;ix<nx;ix++)
               {
               snd1[ib] = pvptr[ix+iz*nx + iv*nx*nz];

               ib++;
               }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(iy=ny-2;iy<ny;iy++)
      {
      pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

      for(iv=0;iv<nv;iv++)
         {
         for(iz=0;iz<nz;iz++)
            {
            for(ix=0;ix<nx;ix++)
               {
               pvptr[ix + iz*nx + iv*nx*nz] = rcv1[ib];

               ib++;
               }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

if(ni->minusId_z >= 0)
   {
   id = ni->minusId_z;
   rcvtag = it + id;

   nb = 2*nv*nx*ny;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(iz=2;iz<4;iz++)
      {
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       snd1[ib] = pvptr[ix + iz*nx + iv*nx*nz];

	       ib++;
	       }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(iz=0;iz<2;iz++)
      {
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       pvptr[ix + iz*nx + iv*nx*nz] = rcv1[ib];

	       ib++;
	       }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

if(ni->plusId_z >= 0)
   {
   id = ni->plusId_z;
   rcvtag = it + id;

   nb = 2*nv*nx*ny;
   blen = nb*sizeof(float);
   snd1 = tbuf;
   rcv1 = tbuf + nb;

   ib = 0;
   for(iz=nz-4;iz<nz-2;iz++)
      {
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       snd1[ib] = pvptr[ix + iz*nx + iv*nx*nz];

	       ib++;
	       }
            }
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

   ib = 0;
   for(iz=nz-2;iz<nz;iz++)
      {
      for(iy=0;iy<ny;iy++)
         {
         pvptr = pvf + iy*N_WAVE_VARS*nx*nz + noff;

         for(iv=0;iv<nv;iv++)
            {
            for(ix=0;ix<nx;ix++)
               {
	       pvptr[ix + iz*nx + iv*nx*nz] = rcv1[ib];

	       ib++;
	       }
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

return(tot_blen);
}
