#include <stdio.h>
#include "pmcl3d.h"

int inisource(int      rank,    int     IFAULT, int     NSRC,   int     READ_STEP, int     NST,     int     *SRCPROC, int    NZ, 
              MPI_Comm MCW,     int     nxt,    int     nyt,    int     nzt,       int     *coords, int     maxdim,   int    *NPSRC,
              PosInf   *ptpsrc, Grid1D  *ptaxx, Grid1D  *ptayy, Grid1D  *ptazz,    Grid1D  *ptaxz,  Grid1D  *ptayz,   Grid1D *ptaxy, char *INSRC)
{
   int i, j, k, npsrc, srcproc, master=0;
   int nbx, nex, nby, ney, nbz, nez;
   PosInf tpsrc=NULL, tpsrcp =NULL;
   Grid1D taxx =NULL, tayy   =NULL, tazz =NULL, taxz =NULL, tayz =NULL, taxy =NULL;
   Grid1D taxxp=NULL, tayyp  =NULL, tazzp=NULL, taxzp=NULL, tayzp=NULL, taxyp=NULL;
   if(NSRC<1) return 0;

   npsrc   = 0;
   srcproc = -1;
// Indexing is based on 1: [1, nxt], etc. Include 1st layer ghost cells
   nbx     = nxt*coords[0] + 1 - 2*loop;
   nex     = nbx + nxt + 4*loop - 1;
   nby     = nyt*coords[1] + 1 - 2*loop;
   ney     = nby + nyt + 4*loop - 1;
   nbz     = 1;
   nez     = nzt;
   // IFAULT=1 has bug! READ_STEP does not work, it tries to read NST all at once - Efe
   if(IFAULT<=1 || IFAULT==5)
   {
	  printf("Reading source.\n");
	  fflush(stdout);
      tpsrc = Alloc1P(NSRC*maxdim);
      taxx  = Alloc1D(NSRC*READ_STEP);
      tayy  = Alloc1D(NSRC*READ_STEP);
      tazz  = Alloc1D(NSRC*READ_STEP);
      taxz  = Alloc1D(NSRC*READ_STEP);
      tayz  = Alloc1D(NSRC*READ_STEP);
      taxy  = Alloc1D(NSRC*READ_STEP);

      if(rank==master)
      {
      	 FILE   *file;
         int    tmpsrc[3];
         Grid1D tmpta;
         if(IFAULT == 1 || IFAULT==5 ){
          file = fopen(INSRC,"rb");
          tmpta = Alloc1D(NST*6);
         }
         else if(IFAULT == 0) file = fopen(INSRC,"r");
         if(!file)
         {
            printf("can't open file %s", INSRC);
	    return 0;
         }
         if(IFAULT == 1 || IFAULT==5 ){
          for(i=0;i<NSRC;i++)
          { 
            if(fread(tmpsrc,sizeof(int),3,file) && fread(tmpta,sizeof(float),NST*6,file))
            {
               tpsrc[i*maxdim]   = tmpsrc[0];
               tpsrc[i*maxdim+1] = tmpsrc[1];
               tpsrc[i*maxdim+2] = NZ+1-tmpsrc[2];
               for(j=0;j<READ_STEP;j++)
               {
                  taxx[i*READ_STEP+j] = tmpta[j*6];
                  tayy[i*READ_STEP+j] = tmpta[j*6+1];
                  tazz[i*READ_STEP+j] = tmpta[j*6+2];
                  taxz[i*READ_STEP+j] = tmpta[j*6+3];
                  tayz[i*READ_STEP+j] = tmpta[j*6+4];
                  taxy[i*READ_STEP+j] = tmpta[j*6+5];
               }
            }
          }
          Delloc1D(tmpta);
         }
         else if(IFAULT == 0)
          for(i=0;i<NSRC;i++)
          {
            fscanf(file, " %d %d %d ",&tmpsrc[0], &tmpsrc[1], &tmpsrc[2]); 
            tpsrc[i*maxdim]   = tmpsrc[0];
            tpsrc[i*maxdim+1] = tmpsrc[1];
            tpsrc[i*maxdim+2] = NZ+1-tmpsrc[2];
            //printf("SOURCE: %d,%d,%d\n",tpsrc[0],tpsrc[1],tpsrc[2]);
            for(j=0;j<READ_STEP;j++){
              fscanf(file, " %f, %f, %f, %f, %f, %f ",
              //fscanf(file, " %f %f %f %f %f %f ",
                &taxx[i*READ_STEP+j], &tayy[i*READ_STEP+j],
                &tazz[i*READ_STEP+j], &taxz[i*READ_STEP+j],
                &tayz[i*READ_STEP+j], &taxy[i*READ_STEP+j]);
              //printf("SOURCE VAL %d: %f,%f\n",j,taxx[j],tayy[j]);
            }
          }
         fclose(file);
      }
      MPI_Bcast(tpsrc, NSRC*maxdim,    MPI_INT,  master, MCW);
      MPI_Bcast(taxx,  NSRC*READ_STEP, MPI_FLOAT, master, MCW); 
      MPI_Bcast(tayy,  NSRC*READ_STEP, MPI_FLOAT, master, MCW);
      MPI_Bcast(tazz,  NSRC*READ_STEP, MPI_FLOAT, master, MCW);
      MPI_Bcast(taxz,  NSRC*READ_STEP, MPI_FLOAT, master, MCW);
      MPI_Bcast(tayz,  NSRC*READ_STEP, MPI_FLOAT, master, MCW);
      MPI_Bcast(taxy,  NSRC*READ_STEP, MPI_FLOAT, master, MCW);
      for(i=0;i<NSRC;i++)
      {
          if( tpsrc[i*maxdim]   >= nbx && tpsrc[i*maxdim]   <= nex && tpsrc[i*maxdim+1] >= nby 
           && tpsrc[i*maxdim+1] <= ney && tpsrc[i*maxdim+2] >= nbz && tpsrc[i*maxdim+2] <= nez)
          {
              srcproc = rank;
              npsrc ++;
          }
      }
      if(npsrc > 0)
      {
          tpsrcp = Alloc1P(npsrc*maxdim);
          taxxp  = Alloc1D(npsrc*READ_STEP);
          tayyp  = Alloc1D(npsrc*READ_STEP);
          tazzp  = Alloc1D(npsrc*READ_STEP);
          taxzp  = Alloc1D(npsrc*READ_STEP);
          tayzp  = Alloc1D(npsrc*READ_STEP);
          taxyp  = Alloc1D(npsrc*READ_STEP);
          k      = 0;
          for(i=0;i<NSRC;i++)
          {
              if( tpsrc[i*maxdim]   >= nbx && tpsrc[i*maxdim]   <= nex && tpsrc[i*maxdim+1] >= nby
               && tpsrc[i*maxdim+1] <= ney && tpsrc[i*maxdim+2] >= nbz && tpsrc[i*maxdim+2] <= nez)
               {
                 tpsrcp[k*maxdim]   = tpsrc[i*maxdim]   - nbx - 1;
                 tpsrcp[k*maxdim+1] = tpsrc[i*maxdim+1] - nby - 1;
                 tpsrcp[k*maxdim+2] = tpsrc[i*maxdim+2] - nbz + 1;
				printf("Source at coordinates tpsrcp[k*maxdim]=%d, tpsrcp[k*maxdim+1]=%d, tpsrcp[k*maxdim+2]=%d\n", tpsrcp[k*maxdim],  tpsrcp[k*maxdim+1], tpsrcp[k*maxdim+2]);
                 for(j=0;j<READ_STEP;j++)
                 {
                    taxxp[k*READ_STEP+j] = taxx[i*READ_STEP+j];
                    tayyp[k*READ_STEP+j] = tayy[i*READ_STEP+j];
                    tazzp[k*READ_STEP+j] = tazz[i*READ_STEP+j];
                    taxzp[k*READ_STEP+j] = taxz[i*READ_STEP+j];
                    tayzp[k*READ_STEP+j] = tayz[i*READ_STEP+j];
                    taxyp[k*READ_STEP+j] = taxy[i*READ_STEP+j];
                  }
                  k++;
               }
          }
      }
      Delloc1D(taxx);
      Delloc1D(tayy);
      Delloc1D(tazz);
      Delloc1D(taxz);
      Delloc1D(tayz);
      Delloc1D(taxy); 
      Delloc1P(tpsrc);       

      *SRCPROC = srcproc;
      *NPSRC   = npsrc;
      *ptpsrc  = tpsrcp;
      *ptaxx   = taxxp;
      *ptayy   = tayyp;
      *ptazz   = tazzp;
      *ptaxz   = taxzp;
      *ptayz   = tayzp;
      *ptaxy   = taxyp;
   }
   return 0;
}

void addsrc(int i,      float DH,   float DT,   int NST,    int npsrc,  int READ_STEP, int dim, int igreen, int nzt, PosInf psrc,
            Grid3D d1,  Grid3D u1,  Grid3D v1,  Grid3D w1,
            Grid1D axx, Grid1D ayy, Grid1D azz, Grid1D axz, Grid1D ayz, Grid1D axy,
            Grid3D xx,  Grid3D yy,  Grid3D zz,  Grid3D xy,  Grid3D yz,  Grid3D xz)
{
  float vtst, vtst1, tmp;
  int idx, idy, idz, j;
  vtst = DT/(DH*DH*DH);
  vtst1 = 1.0/(DH*DH);

  i   = i - 1;
  for(j=0;j<npsrc;j++) 
  {
     idx = psrc[j*dim]   + 1 + 4*loop;
     idy = psrc[j*dim+1] + 1 + 4*loop;
     idz = psrc[j*dim+2] + awp_align - 1;
	 printf("Adding source at array indices (%d, %d, %d)\n", idx, idy, idz);
     if(igreen == -1)
     {
        xx[idx][idy][idz] = xx[idx][idy][idz] - vtst*axx[j*READ_STEP+i];
        yy[idx][idy][idz] = yy[idx][idy][idz] - vtst*ayy[j*READ_STEP+i];
        zz[idx][idy][idz] = zz[idx][idy][idz] - vtst*azz[j*READ_STEP+i];
        xz[idx][idy][idz] = xz[idx][idy][idz] - vtst*axz[j*READ_STEP+i];
        yz[idx][idy][idz] = yz[idx][idy][idz] - vtst*ayz[j*READ_STEP+i];
        xy[idx][idy][idz] = xy[idx][idy][idz] - vtst*axy[j*READ_STEP+i];
     }
     else if(igreen == 1)
     {
        u1[idx][idy][idz] = u1[idx][idy][idz] + vtst*axx[j*READ_STEP+i]/d1[idx][idy][idz];
     }
     else if(igreen == 2)
     {
        v1[idx][idy][idz] = v1[idx][idy][idz] + vtst*ayy[j*READ_STEP+i]/d1[idx][idy][idz]; 
     }
     else if(igreen == 3)
     {
        w1[idx][idy][idz] = w1[idx][idy][idz] + vtst*azz[j*READ_STEP+i]/d1[idx][idy][idz];
     }
     else if(igreen == -2)
     {
        u1[idx][idy][idz] = u1[idx][idy][idz] + vtst*axx[j*READ_STEP+i]/d1[idx][idy][idz];
        v1[idx][idy][idz] = v1[idx][idy][idz] + vtst*ayy[j*READ_STEP+i]/d1[idx][idy][idz];
        w1[idx][idy][idz] = w1[idx][idy][idz] + vtst*azz[j*READ_STEP+i]/d1[idx][idy][idz];
     }
     else if(igreen == 4){
        tmp = vtst1*axx[j*READ_STEP+i];
        xz[idx][idy][nzt+awp_align-1] = tmp;
        tmp = tmp*2;   
        xz[idx][idy][nzt+awp_align] = tmp - xz[idx][idy][nzt+awp_align-2];   
        xz[idx][idy][nzt+awp_align+1] = tmp - xz[idx][idy][nzt+awp_align-3];   
     }
     else if(igreen == 5){
        tmp = vtst1*ayy[j*READ_STEP+i];
        yz[idx][idy][nzt+awp_align-1] = tmp;
        tmp = tmp*2;   
        yz[idx][idy][nzt+awp_align] = tmp - yz[idx][idy][nzt+awp_align-2];   
        yz[idx][idy][nzt+awp_align+1] = tmp - yz[idx][idy][nzt+awp_align-3];   
     }
     else if(igreen == 6){
        tmp = 2.0*vtst1*azz[j*READ_STEP+i];
        zz[idx][idy][nzt+awp_align] = tmp - zz[idx][idy][nzt+awp_align-1];   
        zz[idx][idy][nzt+awp_align+1] = tmp - zz[idx][idy][nzt+awp_align-2];   
     }
/*
     printf("xx=%1.6g\n",xx[idx][idy][idz]);
     printf("yy=%1.6g\n",yy[idx][idy][idz]);
     printf("zz=%1.6g\n",zz[idx][idy][idz]);
     printf("xz=%1.6g\n",xz[idx][idy][idz]);
     printf("yz=%1.6g\n",yz[idx][idy][idz]);
     printf("xy=%1.6g\n",xy[idx][idy][idz]);
*/
  }
  return;
}
