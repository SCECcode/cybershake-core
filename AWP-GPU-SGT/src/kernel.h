#ifndef _KERNEL_H
#define _KERNEL_H

__global__ void dvelcx(float* u1,    float* v1,    float* w1,    float* xx,  float* yy, float* zz, float* xy, float* xz, float* yz,
                       float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, int s_i,   int e_i);

__global__ void dvelcy(float* u1,    float* v1,    float* w1,    float* xx,  float* yy,   float* zz,   float* xy, float* xz, float* yz,
                       float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, float* s_u1, float* s_v1, float* s_w1, int s_j, int e_j);

__global__ void update_boundary_y(float* u1, float* v1, float* w1, float* s_u1, float* s_v1, float* s_w1, int rank, int flag);

__global__ void fvelxyz(float* u1, float* v1, float* w1, float* lam_mu, int xls, int NX, int rankx);

__global__ void dstrqc(float* xx, float* yy,    float* zz,    float* xy,    float* xz,     float* yz,
                       float* r1, float* r2,    float* r3,    float* r4,    float* r5,     float* r6,
                       float* u1, float* v1,    float* w1,    float* lam,   float* mu,     float* qp,
                       float* qs, float* dcrjx, float* dcrjy, float* dcrjz, float* lam_mu, int NX,
                       int rankx, int ranky,    int s_i,      int e_i,      int s_j);

__global__ void ComputeSGT_cu(float* xx,   float* yy,    float* zz,    float* xy,    float* xz,     float* yz,
                              float* sg1,  float* sg2,   float* mu,    int sgt_numsta, int* sgt_sta,float* sgtBuf,
                              int SGT_BLOCK_SIZE,        int SGT_NUMBLOCKS);

__global__ void addsrc_cu(int i,      int READ_STEP, int dim,    int npsrc, int* psrc,  int igreen, int nzt,
                          float* d1,  float* u1,     float* v1,  float* w1,
                          float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
                          float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz);
#endif
