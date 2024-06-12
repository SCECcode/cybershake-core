/*
********************************************************************************
* Grid3D.c                                                                     *
* programming in C language                                                    *
* 3D data structure                                                            * 
********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pmcl3d.h"

Grid3D Alloc3D(int nx, int ny, int nz)
{
   int i, j, k;
   Grid3D U = (Grid3D)malloc(sizeof(float**)*nx + sizeof(float *)*nx*ny +sizeof(float)*nx*ny*nz);

   if (!U){
       printf("Cannot allocate 3D float array\n");
       exit(-1);
   }
   for(i=0;i<nx;i++){
       U[i] = ((float**) U) + nx + i*ny;
    }

   float *Ustart = (float *) (U[nx-1] + ny);
   for(i=0;i<nx;i++)
       for(j=0;j<ny;j++)
           U[i][j] = Ustart + i*ny*nz + j*nz;

   for(i=0;i<nx;i++)
       for(j=0;j<ny;j++)
           for(k=0;k<nz;k++)
              U[i][j][k] = 0.0f;

   return U;
}


Grid1D Alloc1D(int nx)
{
   int i;
   Grid1D U = (Grid1D)malloc(sizeof(float)*nx);

   if (!U){
       printf("Cannot allocate 2D float array\n");
       exit(-1);
   }

   for(i=0;i<nx;i++)
       U[i] = 0.0f;

   return U;
}


PosInf Alloc1P(int nx)
{
   int i;
   PosInf U = (PosInf)malloc(sizeof(int)*nx);

   if (!U){
       printf("Cannot allocate 2D integer array\n");
       exit(-1);
   }

   for(i=0;i<nx;i++)
       U[i] = 0;

   return U;
}

void Delloc3D(Grid3D U)
{
   if (U) 
   {
      free(U);
      U = NULL;
   }

   return;
}

void Delloc1D(Grid1D U)
{
   if (U)
   {
      free(U);
      U = NULL;
   }

   return;
}

void Delloc1P(PosInf U)
{
   if (U)
   {
      free(U);
      U = NULL;
   }

   return;
}
