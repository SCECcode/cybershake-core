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

void mpi_global_vel_minmax(struct vel_minmax *vpars,struct nodeinfo *ni)
{
MPI_Status s;
struct vel_minmax vtmp, vorig;
int ib;
int srtag, blen;

vorig.vpmax = vpars->vpmax;
vorig.ix_vpmax = vpars->ix_vpmax;
vorig.iy_vpmax = vpars->iy_vpmax;
vorig.iz_vpmax = vpars->iz_vpmax;

vorig.vpmin = vpars->vpmin;
vorig.ix_vpmin = vpars->ix_vpmin;
vorig.iy_vpmin = vpars->iy_vpmin;
vorig.iz_vpmin = vpars->iz_vpmin;

vorig.vsmax = vpars->vsmax;
vorig.ix_vsmax = vpars->ix_vsmax;
vorig.iy_vsmax = vpars->iy_vsmax;
vorig.iz_vsmax = vpars->iz_vsmax;

vorig.vsmin = vpars->vsmin;
vorig.ix_vsmin = vpars->ix_vsmin;
vorig.iy_vsmin = vpars->iy_vsmin;
vorig.iz_vsmin = vpars->iz_vsmin;

vorig.vpvsmax = vpars->vpvsmax;
vorig.ix_vpvsmax = vpars->ix_vpvsmax;
vorig.iy_vpvsmax = vpars->iy_vpvsmax;
vorig.iz_vpvsmax = vpars->iz_vpvsmax;

vorig.vpvsmin = vpars->vpvsmin;
vorig.ix_vpvsmin = vpars->ix_vpvsmin;
vorig.iy_vpvsmin = vpars->iy_vpvsmin;
vorig.iz_vpvsmin = vpars->iz_vpvsmin;

blen = sizeof(struct vel_minmax);
srtag = ni->segmentId;

if(ni->segmentId != 0)
   {
   MPI_Send(vpars,blen,MPI_CHAR,0,srtag,MPI_COMM_WORLD);
   MPI_Recv(vpars,blen,MPI_CHAR,0,srtag,MPI_COMM_WORLD,&s);
   }
else
   {
   for(ib=1;ib<ni->nproc;ib++)
      {
      MPI_Recv(&vtmp,blen,MPI_CHAR,ib,ib,MPI_COMM_WORLD,&s);

      if(vtmp.vpmax > vpars->vpmax)
         {
	 vpars->vpmax = vtmp.vpmax;
	 vpars->ix_vpmax = vtmp.ix_vpmax;
	 vpars->iy_vpmax = vtmp.iy_vpmax;
	 vpars->iz_vpmax = vtmp.iz_vpmax;
	 }

      if(vtmp.vpmin < vpars->vpmin)
         {
	 vpars->vpmin = vtmp.vpmin;
	 vpars->ix_vpmin = vtmp.ix_vpmin;
	 vpars->iy_vpmin = vtmp.iy_vpmin;
	 vpars->iz_vpmin = vtmp.iz_vpmin;
	 }

      if(vtmp.vsmax > vpars->vsmax)
         {
	 vpars->vsmax = vtmp.vsmax;
	 vpars->ix_vsmax = vtmp.ix_vsmax;
	 vpars->iy_vsmax = vtmp.iy_vsmax;
	 vpars->iz_vsmax = vtmp.iz_vsmax;
	 }

      if(vtmp.vsmin < vpars->vsmin)
         {
	 vpars->vsmin = vtmp.vsmin;
	 vpars->ix_vsmin = vtmp.ix_vsmin;
	 vpars->iy_vsmin = vtmp.iy_vsmin;
	 vpars->iz_vsmin = vtmp.iz_vsmin;
	 }

      if(vtmp.vpvsmax > vpars->vpvsmax)
         {
	 vpars->vpvsmax = vtmp.vpvsmax;
	 vpars->ix_vpvsmax = vtmp.ix_vpvsmax;
	 vpars->iy_vpvsmax = vtmp.iy_vpvsmax;
	 vpars->iz_vpvsmax = vtmp.iz_vpvsmax;
	 }

      if(vtmp.vpvsmin < vpars->vpvsmin)
         {
	 vpars->vpvsmin = vtmp.vpvsmin;
	 vpars->ix_vpvsmin = vtmp.ix_vpvsmin;
	 vpars->iy_vpvsmin = vtmp.iy_vpvsmin;
	 vpars->iz_vpvsmin = vtmp.iz_vpvsmin;
	 }
      }

   for(ib=1;ib<ni->nproc;ib++)
      MPI_Send(vpars,blen,MPI_CHAR,ib,ib,MPI_COMM_WORLD);
   }

fprintf(stderr,"\t resolving Vpmax: \tsend= %g \trecv= %g",vorig.vpmax,vpars->vpmax);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vpmax,vpars->iy_vpmax,vpars->iz_vpmax);

fprintf(stderr,"\t resolving Vpmin: \tsend= %g \trecv= %g",vorig.vpmin,vpars->vpmin);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vpmin,vpars->iy_vpmin,vpars->iz_vpmin);

fprintf(stderr,"\t resolving Vsmax: \tsend= %g \trecv= %g",vorig.vsmax,vpars->vsmax);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vsmax,vpars->iy_vsmax,vpars->iz_vsmax);

fprintf(stderr,"\t resolving Vsmin: \tsend= %g \trecv= %g",vorig.vsmin,vpars->vsmin);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vsmin,vpars->iy_vsmin,vpars->iz_vsmin);

fprintf(stderr,"\t resolving VpVsmax: \tsend= %g \trecv= %g",vorig.vpvsmax,vpars->vpvsmax);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vpvsmax,vpars->iy_vpvsmax,vpars->iz_vpvsmax);

fprintf(stderr,"\t resolving VpVsmin: \tsend= %g \trecv= %g",vorig.vpvsmin,vpars->vpvsmin);
fprintf(stderr,"\t(ix= %5d iy= %5d iz= %5d)\n",vpars->ix_vpvsmin,vpars->iy_vpvsmin,vpars->iz_vpvsmin);

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

size_t exchange1_wptrs_pvP3(float *pvf,float *tbuf,struct nodeinfo *ni,int it,int iflag,struct tms *dtu0,int pflag)
{
MPI_Status s;
size_t tot_blen;
float *pvptr, *snd1, *rcv1;
int ix, iy, iz, iv, nx, ny, nz, nv, ib, nb, noff;
int blen, id, sndtag, rcvtag;

float *pvptr1, *tptr1, *pvptr2, *tptr2;
int blen2, ix1, ix2, iy1, iy2, iz1, iz2;

struct tms use0;
clock_t rt0;

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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

   ix1 = 2;
   ix2 = 3;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

         for(iz=0;iz<nz;iz++)
            {
            tptr1[0] = pvptr1[0];
            tptr2[0] = pvptr2[0];

            tptr1++;
            tptr2++;
            pvptr1 += nx;
            pvptr2 += nx;
            }
         }
      }

   rt0 = times(&use0);
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

   ix1 = 0;
   ix2 = 1;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

         for(iz=0;iz<nz;iz++)
            {
            pvptr1[0] = tptr1[0];
            pvptr2[0] = tptr2[0];

            tptr1++;
            tptr2++;
            pvptr1 += nx;
            pvptr2 += nx;
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

   ix1 = nx-4;
   ix2 = nx-3;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

         for(iz=0;iz<nz;iz++)
            {
            tptr1[0] = pvptr1[0];
            tptr2[0] = pvptr2[0];

            tptr1++;
            tptr2++;
            pvptr1 += nx;
            pvptr2 += nx;
            }
         }
      }

   rt0 = times(&use0);
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

   ix1 = nx-2;
   ix2 = nx-1;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

         for(iz=0;iz<nz;iz++)
            {
            pvptr1[0] = tptr1[0];
            pvptr2[0] = tptr2[0];

            tptr1++;
            tptr2++;
            pvptr1 += nx;
            pvptr2 += nx;
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

   blen2 = nv*nx*nz*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

/*
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

         iy1 = 2;
         iy2 = 3;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   rt0 = times(&use0);
   for(iy=2;iy<4;iy++)
      {
      snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff;
      rcv1 = pvf + (iy-2)*N_WAVE_VARS*nx*nz + noff;

      MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
      }
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   /*
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;
   */

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

         iy1 = 0;
         iy2 = 1;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

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

   blen2 = nv*nx*nz*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

/*
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

         iy1 = ny-4;
         iy2 = ny-3;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   rt0 = times(&use0);
   for(iy=ny-4;iy<ny-2;iy++)
      {
      snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff;
      rcv1 = pvf + (iy+2)*N_WAVE_VARS*nx*nz + noff;

      MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
      }
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   /*
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;
   */

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

         iy1 = ny-2;
         iy2 = ny-1;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

   iz1 = 2;
   iz2 = 3;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         for(ix=0;ix<nx;ix++)
            {
            tptr1[0] = pvptr1[0];
            tptr2[0] = pvptr2[0];

            tptr1++;
            tptr2++;
            pvptr1++;
            pvptr2++;
            }
         }
      }

   rt0 = times(&use0);
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

   iz1 = 0;
   iz2 = 1;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         for(ix=0;ix<nx;ix++)
            {
            pvptr1[0] = tptr1[0];
            pvptr2[0] = tptr2[0];

            tptr1++;
            tptr2++;
            pvptr1++;
            pvptr2++;
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

   iz1 = nz-4;
   iz2 = nz-3;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         for(ix=0;ix<nx;ix++)
            {
            tptr1[0] = pvptr1[0];
            tptr2[0] = pvptr2[0];

            tptr1++;
            tptr2++;
            pvptr1++;
            pvptr2++;
            }
         }
      }

   rt0 = times(&use0);
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   if(pflag)
      dtimes_lite(use0,rt0,dtu0);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

   iz1 = nz-2;
   iz2 = nz-1;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         for(ix=0;ix<nx;ix++)
            {
            pvptr1[0] = tptr1[0];
            pvptr2[0] = tptr2[0];

            tptr1++;
            tptr2++;
            pvptr1++;
            pvptr2++;
            }
         }
      }

   tot_blen = tot_blen + 2*blen;
   }

return(tot_blen);
}

size_t exchange1_nocpy_pvP3(float *pvf,float *tbuf,struct nodeinfo *ni,int it,int iflag)
{
MPI_Status s;
size_t tot_blen;
float *pvptr, *snd1, *rcv1;
int ix, iy, iz, iv, nx, ny, nz, nv, ib, nb, noff;
int blen, id, sndtag, rcvtag;

float *pvptr1, *tptr1, *pvptr2, *tptr2;
int blen2, ix1, ix2, iy1, iy2, iz1, iz2;

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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

/*
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
*/

         ix1 = 2;
         ix2 = 3;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

               for(iz=0;iz<nz;iz++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1 += nx;
                  pvptr2 += nx;
                  }
               }
            }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

/*
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
*/

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

         ix1 = 0;
         ix2 = 1;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

               for(iz=0;iz<nz;iz++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1 += nx;
                  pvptr2 += nx;
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

/*
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
*/

         ix1 = nx-4;
         ix2 = nx-3;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

               for(iz=0;iz<nz;iz++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1 += nx;
                  pvptr2 += nx;
                  }
               }
            }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
   */

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;

/*
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
*/

   tptr1 = tbuf;
   tptr2 = tbuf + nv*ny*nz;

         ix1 = nx-2;
         ix2 = nx-1;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix1 + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + ix2 + iv*nx*nz;

               for(iz=0;iz<nz;iz++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1 += nx;
                  pvptr2 += nx;
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

   blen2 = nv*nx*nz*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

/*
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

         iy1 = 2;
         iy2 = 3;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   for(iy=2;iy<4;iy++)
      {
      snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff;
      rcv1 = pvf + (iy-2)*N_WAVE_VARS*nx*nz + noff;

      MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
      }

/*
   iy1 = 2;
   pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff;
   iy2 = 0;
   pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff;

   MPI_Sendrecv(pvptr1,blen2,MPI_CHAR,id,sndtag,pvptr2,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);

   iy1 = 3;
   pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff;
   iy2 = 1;
   pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff;

   MPI_Sendrecv(pvptr1,blen2,MPI_CHAR,id,sndtag,pvptr2,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
*/

   /*
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   */
   rcv1 = tbuf;

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

         iy1 = 0;
         iy2 = 1;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

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

   blen2 = nv*nx*nz*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

/*
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

         iy1 = ny-4;
         iy2 = ny-3;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   for(iy=ny-4;iy<ny-2;iy++)
      {
      snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff;
      rcv1 = pvf + (iy+2)*N_WAVE_VARS*nx*nz + noff;

      MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
      }

/*
   iy1 = ny-4;
   pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff;
   iy2 = ny-2;
   pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff;

   MPI_Sendrecv(pvptr1,blen2,MPI_CHAR,id,sndtag,pvptr2,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);

   iy1 = ny-3;
   pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff;
   iy2 = ny-1;
   pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff;

   MPI_Sendrecv(pvptr1,blen2,MPI_CHAR,id,sndtag,pvptr2,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
*/

   /*
   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;
   */

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*nz;

         iy1 = ny-2;
         iy2 = ny-1;
         for(iv=0;iv<nv;iv++)
            {
            for(iz=0;iz<nz;iz++)
               {
               pvptr1 = pvf + iy1*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;
               pvptr2 = pvf + iy2*N_WAVE_VARS*nx*nz + noff + iz*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

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

   blen2 = 2*nx*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

/*
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

         iz1 = 2;
         iz2 = 3;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   iz1 = 2;
   iz2 = 0;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         rcv1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;
   */

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

         iz1 = 0;
         iz2 = 1;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

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

   blen2 = 2*nx*sizeof(float);

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

/*
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

         iz1 = nz-4;
         iz2 = nz-3;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  tptr1[0] = pvptr1[0];
                  tptr2[0] = pvptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   iz1 = nz-4;
   iz2 = nz-2;
   for(iy=0;iy<ny;iy++)
      {
      for(iv=0;iv<nv;iv++)
         {
         snd1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
         rcv1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

         MPI_Sendrecv(snd1,blen2,MPI_CHAR,id,sndtag,rcv1,blen2,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);
         }
      }

   /*
   MPI_Sendrecv(snd1,blen,MPI_CHAR,id,sndtag,rcv1,blen,MPI_CHAR,id,rcvtag,MPI_COMM_WORLD,&s);

   MPI_Sendrecv_replace(tbuf,blen,MPI_CHAR,id,sndtag,id,rcvtag,MPI_COMM_WORLD,&s);
   rcv1 = tbuf;
   */

/*
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

   tptr1 = tbuf;
   tptr2 = tbuf + nv*nx*ny;

         iz1 = nz-2;
         iz2 = nz-1;
         for(iy=0;iy<ny;iy++)
            {
            for(iv=0;iv<nv;iv++)
               {
               pvptr1 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz1*nx + iv*nx*nz;
               pvptr2 = pvf + iy*N_WAVE_VARS*nx*nz + noff + iz2*nx + iv*nx*nz;

               for(ix=0;ix<nx;ix++)
                  {
                  pvptr1[0] = tptr1[0];
                  pvptr2[0] = tptr2[0];

                  tptr1++;
                  tptr2++;
                  pvptr1++;
                  pvptr2++;
                  }
               }
            }
*/

   tot_blen = tot_blen + 2*blen;
   }

return(tot_blen);
}

void mpi_global_prof1d(struct media_input *mi,struct nodeinfo *ni)
{
MPI_Status s;
struct profile1d *prof1d, *ptmp;
float afac, fac2;
int iz, ib;
int srtag, blen;

prof1d = mi->prof1d;

afac = 1.0/((1.0*ni->nproc_x)*(1.0*ni->nproc_y));
blen = ni->globnz*sizeof(struct profile1d);
srtag = ni->segmentId;

if(ni->segmentId != 0)
   {
   MPI_Send(prof1d,blen,MPI_CHAR,0,srtag,MPI_COMM_WORLD);
   MPI_Recv(prof1d,blen,MPI_CHAR,0,srtag,MPI_COMM_WORLD,&s);
   }
else
   {
   ptmp = (struct profile1d *) check_malloc (blen);

   for(iz=0;iz<ni->globnz;iz++)
      {
      if(prof1d[iz].vp > 0.0)
         {
         prof1d[iz].vp = afac*prof1d[iz].vp;
         prof1d[iz].vs = afac*prof1d[iz].vs;
         prof1d[iz].dn = afac*prof1d[iz].dn;
         prof1d[iz].qp = afac*prof1d[iz].qp;
         prof1d[iz].qs = afac*prof1d[iz].qs;
         }
      else
         {
         prof1d[iz].vp = 0.0;
         prof1d[iz].vs = 0.0;
         prof1d[iz].dn = 0.0;
         prof1d[iz].qp = 0.0;
         prof1d[iz].qs = 0.0;
         }
      }

   for(ib=1;ib<ni->nproc;ib++)
      {
      MPI_Recv(ptmp,blen,MPI_CHAR,ib,ib,MPI_COMM_WORLD,&s);

      for(iz=0;iz<ni->globnz;iz++)
	 {
	 if(ptmp[iz].vp > 0.0)
	    {
	    prof1d[iz].vp = prof1d[iz].vp + afac*ptmp[iz].vp;
	    prof1d[iz].vs = prof1d[iz].vs + afac*ptmp[iz].vs;
	    prof1d[iz].dn = prof1d[iz].dn + afac*ptmp[iz].dn;
	    prof1d[iz].qp = prof1d[iz].qp + afac*ptmp[iz].qp;
	    prof1d[iz].qs = prof1d[iz].qs + afac*ptmp[iz].qs;
	    }
	 }
      }

   for(iz=0;iz<ni->globnz;iz++)
      {
      if(prof1d[iz].vp/prof1d[iz].vs > mi->vpvs_max_boundary)    /* check vp/vs ratio */
         prof1d[iz].vp = prof1d[iz].vs*(mi->vpvs_max_boundary);

      if(prof1d[iz].vp < mi->vp_min_boundary)   /* check vp_min, preserve vp/vs ratio */
         {
         fac2 = prof1d[iz].vs/prof1d[iz].vp;

         prof1d[iz].vp = mi->vp_min_boundary;
         prof1d[iz].vs = fac2*prof1d[iz].vs;
         }
      }

   for(iz=0;iz<NBND_PAD_TOPCOPY;iz++)
      {
      prof1d[iz].vp = prof1d[NBND_PAD_TOPCOPY].vp;
      prof1d[iz].vs = prof1d[NBND_PAD_TOPCOPY].vs;
      prof1d[iz].dn = prof1d[NBND_PAD_TOPCOPY].dn;
      prof1d[iz].qp = prof1d[NBND_PAD_TOPCOPY].qp;
      prof1d[iz].qs = prof1d[NBND_PAD_TOPCOPY].qs;
      }

/* XXXXYYYY
   for(iz=0;iz<ni->globnz;iz++)
      {
      prof1d[iz].vp = 5.2000000;
      prof1d[iz].vs = 3.0000000;
      prof1d[iz].dn = 2.5000000;
      prof1d[iz].qp = 20000.0;
      prof1d[iz].qs = 10000.0;
      }
 */

   for(ib=1;ib<ni->nproc;ib++)
      MPI_Send(prof1d,blen,MPI_CHAR,ib,ib,MPI_COMM_WORLD);

   free(ptmp);
   }
}
