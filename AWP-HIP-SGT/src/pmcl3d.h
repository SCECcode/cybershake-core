/*  
********************************************************************************
* pmcl3d.h                                                                     *
* programming in C language                                                    *
* all pmcl3d data types are defined here                                       * 
********************************************************************************
*/

#include <mpi.h>
#include "cuda_to_hip.h"
#include "pmcl3d_cons.h"
#include "kernel.h"

#ifdef __RESTRICT 
#define RESTRICT restrict 
#else
#define RESTRICT 
#endif

#ifndef _PMCL3D_H
#define _PMCL3D_H

typedef float *RESTRICT *RESTRICT *RESTRICT Grid3D;
typedef float *RESTRICT Grid1D;
typedef int   *RESTRICT PosInf;

int inisgt(int rank, char *INSGT, int *sgt_numsta_all, int *sgt_numsta, int *sgt_max_numsta,
    int **sgt_sta, 
    int **sgt_sta_indices, int *sgtproc, 
    MPI_Comm MCW, int NZ, int nxt, int nyt, int nzt, int *coord);

void command(int argc, char **argv,
             float *TMAX, float *DH, float *DT, float *ARBC, float *PHT,
             int *NPC, int *ND, int *NSRC, int *NST, int *NVAR,
             int *NVE, int *MEDIASTART, int *IFAULT, 
             int *READ_STEP, int *READ_STEP_GPU,
             int *NTISKP, int *WRITE_STEP,
             int *NX, int *NY, int *NZ, int *PX, int *PY, 
             int *NBGX, int *NEDX, int *NSKPX,
             int *NBGY, int *NEDY, int *NSKPY,
             int *NBGZ, int *NEDZ, int *NSKPZ,
             float *FL, float *FH, float *FP, int *IDYNA, int *SoCalQ,
             char  *INSRC,  char *INVEL, char *OUT,
             char  *CHKFILE,
             int *IGREEN, int *NTISKP_SGT, char  *INSGT);

int inisource(int      rank,    int     IFAULT, int     NSRC,   int     READ_STEP, 
              int     NST,     int     *SRCPROC, int    NZ,
              MPI_Comm MCW,     int     nxt,    int     nyt,    int     nzt,       
              int     *coords, int     maxdim,   int    *NPSRC,
              PosInf   *ptpsrc, Grid1D  *ptaxx, Grid1D  *ptayy, Grid1D  *ptazz,    
              Grid1D  *ptaxz,  Grid1D  *ptayz,   Grid1D *ptaxy, char *INSRC);

void addsrc(int i,      float DH,   float DT,   int NST,    int npsrc,  int READ_STEP, int dim, int igreen, int nzt, PosInf psrc,
            Grid3D d1,  Grid3D u1,  Grid3D v1,  Grid3D w1,
            Grid1D axx, Grid1D ayy, Grid1D azz, Grid1D axz, Grid1D ayz, Grid1D axy,
            Grid3D xx,  Grid3D yy,  Grid3D zz,  Grid3D xy,  Grid3D yz,  Grid3D xz);

void inimesh(int MEDIASTART, Grid3D d1, Grid3D mu, Grid3D lam, Grid3D qp, Grid3D qs, float *taumax, float *taumin,
             int nvar, float FP,  float FL, float FH, int nxt, int nyt, int nzt, int PX, int PY, int NX, int NY,
             int NZ, int *coords, MPI_Comm MCW, int IDYNA, int NVE, int SoCalQ, char *INVEL,
             float *vse, float *vpe, float *dde);

int writeCHK(char *chkfile, int ntiskp, float dt, float dh, 
      int nxt, int nyt, int nzt,
      int nt, float arbc, int npc, int nve,
      float fl, float fh, float fp, 
      float *vse, float *vpe, float *dde);

void mediaswap(Grid3D d1, Grid3D mu,     Grid3D lam,    Grid3D qp,     Grid3D qs,
               int rank,  int x_rank_L,  int x_rank_R,  int y_rank_F,  int y_rank_B,
               int nxt,   int nyt,       int nzt,       MPI_Comm MCW);

void tausub( Grid3D tau, float taumin,float taumax);

void inicrj(float ARBC, int *coords, int nxt, int nyt, int nzt, int NX, int NY, int ND, Grid1D dcrjx, Grid1D dcrjy, Grid1D dcrjz);

void init_texture(int nxt,  int nyt,  int nzt,  Grid3D tau1,  Grid3D tau2,  Grid3D vx1,  Grid3D vx2,
                  int xls,  int xre,  int yls,  int yre);

void PostRecvMsg_X(float* RL_M, float* RR_M, MPI_Comm MCW, MPI_Request* request, int* count, int msg_size, int rank_L, int rank_R);

void PostRecvMsg_Y(float* RF_M, float* RB_M, MPI_Comm MCW, MPI_Request* request, int* count, int msg_size, int rank_F, int rank_B);

void Cpy2Device_source(int npsrc, int READ_STEP,
      int index_offset,
      Grid1D taxx, Grid1D tayy, Grid1D tazz,
      Grid1D taxz, Grid1D tayz, Grid1D taxy,
      float *d_taxx, float *d_tayy, float *d_tazz,
      float *d_taxz, float *d_tayz, float *d_taxy);

void Cpy2Host_VX(float* u1, float* v1, float* w1, float* h_m, int nxt, int nyt, int nzt, cudaStream_t St, int rank, int flag);

void Cpy2Host_VY(float* s_u1, float* s_v1, float* s_w1, float* h_m, int nxt, int nzt, cudaStream_t St, int rank);

void Cpy2Device_VX(float* u1, float* v1, float* w1,        float* L_m,       float* R_m, int nxt,    
                   int nyt,   int nzt,   cudaStream_t St1, cudaStream_t St2, int rank_L, int rank_R);

void Cpy2Device_VY(float* u1,   float *v1,  float *w1,  float* f_u1, float* f_v1, float* f_w1, float* b_u1,      float* b_v1,
                   float* b_w1, float* F_m, float* B_m, int nxt,     int nyt,     int nzt,     cudaStream_t St1, cudaStream_t St2,
                   int rank_F,  int rank_B, int rank);

void PostSendMsg_X(float* SL_M, float* SR_M, MPI_Comm MCW, MPI_Request* request, int* count, int msg_size,
                   int rank_L,  int rank_R,  int rank,     int flag);

void PostSendMsg_Y(float* SF_M, float* SB_M, MPI_Comm MCW, MPI_Request* request, int* count, int msg_size,
                   int rank_F,  int rank_B,  int rank,     int flag);

Grid3D Alloc3D(int nx, int ny, int nz);
Grid1D Alloc1D(int nx);    
PosInf Alloc1P(int nx);

void Delloc3D(Grid3D U);
void Delloc1D(Grid1D U);
void Delloc1P(PosInf U);

#endif
