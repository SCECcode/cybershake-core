#include <stdio.h>
#include "pmcl3d.h"

int inisgt(int rank, char *INSGT, int *sgt_numsta_all, int *sgt_numsta, int *sgt_max_numsta,
    int **sgt_sta, 
    int **sgt_sta_indices, int *sgtproc, 
    MPI_Comm MCW, int NZ, int nxt, int nyt, int nzt, int *coord){

  int master=0, i;
  PosInf sgt_sta_all;
  PosInf t_sgt_sta_indices;
  PosInf t_sgt_sta;
  int nbgx, nbgy, nbgz, nedx, nedy, nedz;

  *sgt_numsta = 0;
  nbgx = coord[0]*nxt + 1;
  nbgy = coord[1]*nyt + 1;
  nbgz = 1;
  nedx = nbgx + nxt - 1;
  nedy = nbgy + nyt - 1;
  nedz = nbgz + nzt - 1;
  if(rank%10==0) printf("SGT: inside rank=%d, nbgx,edx:(%d,%d), (%d,%d), (%d,%d)\n",rank,nbgx,nedx,nbgy,nedy,nbgz,nedz);

  if(rank==master){
    FILE *f;
    f = fopen(INSGT, "r");
    if(f){
      fscanf(f," %d",sgt_numsta_all); 
      sgt_sta_all = Alloc1P((*sgt_numsta_all)*3);
      for(i=0;i<*sgt_numsta_all;i++){
        fscanf(f," %d %d %d",&sgt_sta_all[i*3],&sgt_sta_all[i*3+1],&sgt_sta_all[i*3+2]);
      }
    }
    else{
      printf("ERROR! Cannot read INSGT file: %s\n",INSGT);
      return -1;
    }
  }
  MPI_Bcast(sgt_numsta_all, 1,                  MPI_INT,  master, MCW);
  if(rank!=master){ 
    sgt_sta_all = Alloc1P((*sgt_numsta_all)*3);
  }
  MPI_Bcast(sgt_sta_all,    (*sgt_numsta_all)*3,   MPI_INT,  master, MCW);
  t_sgt_sta_indices = Alloc1P(*sgt_numsta_all);
  *sgtproc = -1;
  for(i=0;i<*sgt_numsta_all;i++){
    sgt_sta_all[i*3+2] = NZ+1-sgt_sta_all[i*3+2];
    if(sgt_sta_all[i*3] >= nbgx && sgt_sta_all[i*3] <= nedx 
          && sgt_sta_all[i*3+1] >= nbgy && sgt_sta_all[i*3+1] <= nedy
          && sgt_sta_all[i*3+2] >= nbgz && sgt_sta_all[i*3+2] <= nedz){
      t_sgt_sta_indices[(*sgt_numsta)] = i;
      (*sgt_numsta)++;
      *sgtproc = rank;
    }
  }
  if(*sgtproc == rank) printf("SGT: rank=%d, numsta=%d\n",rank,*sgt_numsta);
  MPI_Allreduce(sgt_numsta, sgt_max_numsta, 1, MPI_INT, MPI_MAX, MCW);
  int ind;
  t_sgt_sta = Alloc1P((*sgt_numsta)*3);
  *sgt_sta_indices = Alloc1P((*sgt_numsta));
  for(i=0;i<*sgt_numsta;i++){
    (*sgt_sta_indices)[i] = t_sgt_sta_indices[i];
    ind = (*sgt_sta_indices)[i];
    t_sgt_sta[i*3]   = sgt_sta_all[ind*3] - nbgx +1;
    t_sgt_sta[i*3+1] = sgt_sta_all[ind*3+1] - nbgy +1;
    t_sgt_sta[i*3+2] = sgt_sta_all[ind*3+2] - nbgz +1;
  }
  *sgt_sta = t_sgt_sta;
  Delloc1P(sgt_sta_all);
  Delloc1P(t_sgt_sta_indices);

return 0;
}
