#include <stdio.h>
#include "kernel.h"
#include "pmcl3d_cons.h"
#include "cuda_to_hip.h"

__constant__ float d_c1;
__constant__ float d_c2;
__constant__ float d_dth;
__constant__ float d_dt1;
__constant__ float d_dh1;
__constant__ float d_DT;
__constant__ float d_DH;
__constant__ int   d_nxt;
__constant__ int   d_nyt;
__constant__ int   d_nzt;
__constant__ int   d_slice_1;
__constant__ int   d_slice_2;
__constant__ int   d_yline_1;
__constant__ int   d_yline_2;

texture<float, 1, cudaReadModeElementType> p_vx1;
texture<float, 1, cudaReadModeElementType> p_vx2;

extern "C"
void SetDeviceConstValue(float DH, float DT, int nxt, int nyt, int nzt)
{
    float h_c1, h_c2, h_dth, h_dt1, h_dh1;
    int   slice_1,  slice_2,  yline_1,  yline_2;
    h_c1  = 9.0/8.0;
    h_c2  = -1.0/24.0;
    h_dth = DT/DH;
    h_dt1 = 1.0/DT;
    h_dh1 = 1.0/DH;
    slice_1  = (nyt+4+8*loop)*(nzt+2*align);
    slice_2  = (nyt+4+8*loop)*(nzt+2*align)*2;
    yline_1  = nzt+2*align;
    yline_2  = (nzt+2*align)*2;
  
    cudaMemcpyToSymbol(d_c1,      &h_c1,    sizeof(float));
    cudaMemcpyToSymbol(d_c2,      &h_c2,    sizeof(float));
    cudaMemcpyToSymbol(d_dth,     &h_dth,   sizeof(float));
    cudaMemcpyToSymbol(d_dt1,     &h_dt1,   sizeof(float));
    cudaMemcpyToSymbol(d_dh1,     &h_dh1,   sizeof(float));
    cudaMemcpyToSymbol(d_DT,      &DT,      sizeof(float));
    cudaMemcpyToSymbol(d_DH,      &DH,      sizeof(float));
    cudaMemcpyToSymbol(d_nxt,     &nxt,     sizeof(int));
    cudaMemcpyToSymbol(d_nyt,     &nyt,     sizeof(int));
    cudaMemcpyToSymbol(d_nzt,     &nzt,     sizeof(int));
    cudaMemcpyToSymbol(d_slice_1, &slice_1, sizeof(int));
    cudaMemcpyToSymbol(d_slice_2, &slice_2, sizeof(int));
    cudaMemcpyToSymbol(d_yline_1, &yline_1, sizeof(int));
    cudaMemcpyToSymbol(d_yline_2, &yline_2, sizeof(int));
    return;
}

extern "C"
void BindArrayToTexture(float* vx1, float* vx2, int memsize)
{
   cudaBindTexture(0, p_vx1,  vx1,  memsize);
   cudaBindTexture(0, p_vx2,  vx2,  memsize);
   cudaThreadSynchronize ();
   return;
}

extern "C"
void UnBindArrayFromTexture()
{
   cudaUnbindTexture(p_vx1);
   cudaUnbindTexture(p_vx2);
   return;
}

extern "C"
void dvelcx_H(float* u1,    float* v1,    float* w1,    float* xx,  float* yy, float* zz, float* xy,      float* xz, float* yz,
             float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, int nyt,   int nzt,   cudaStream_t St, int s_i,   int e_i, int rank)
{
    dim3 block (BLOCK_SIZE_Z, BLOCK_SIZE_Y, 1);
    dim3 grid ((nzt+BLOCK_SIZE_Z-1)/BLOCK_SIZE_Z, (nyt+BLOCK_SIZE_Y-1)/BLOCK_SIZE_Y,1);
    cudaFuncSetCacheConfig(dvelcx, cudaFuncCachePreferL1);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dvelcx, before kernel: %s\n",rank,cudaGetErrorString(err));
    dvelcx<<<grid, block, 0, St>>>(u1, v1, w1, xx, yy, zz, xy, xz, yz, dcrjx, dcrjy, dcrjz, d_1, s_i, e_i);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dvelcx: %s\n",rank,cudaGetErrorString(err));
    return;
}

extern "C"
void dvelcy_H(float* u1,       float* v1,    float* w1,    float* xx,  float* yy, float* zz, float* xy,   float* xz,   float* yz,
              float* dcrjx,    float* dcrjy, float* dcrjz, float* d_1, int nxt,   int nzt,   float* s_u1, float* s_v1, float* s_w1,  
              cudaStream_t St, int s_j,      int e_j,      int rank,   int rank_me)
{
    if(rank==-1) return;
    dim3 block (BLOCK_SIZE_Z, BLOCK_SIZE_Y, 1);
    dim3 grid ((nzt+BLOCK_SIZE_Z-1)/BLOCK_SIZE_Z, (nxt+BLOCK_SIZE_Y-1)/BLOCK_SIZE_Y,1);
    cudaFuncSetCacheConfig(dvelcy, cudaFuncCachePreferL1);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dvelcy, before kernel: %s\n",rank_me,cudaGetErrorString(err));
    dvelcy<<<grid, block, 0, St>>>(u1, v1, w1, xx, yy, zz, xy, xz, yz, dcrjx, dcrjy, dcrjz, d_1, s_u1, s_v1, s_w1, s_j, e_j);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dvelcy: %s\n",rank_me,cudaGetErrorString(err));
    return;
}

extern "C"
void update_bound_y_H(float* u1,   float* v1, float* w1, float* f_u1,      float* f_v1,      float* f_w1,  float* b_u1, float* b_v1, 
                      float* b_w1, int nxt,   int nzt,   cudaStream_t St1, cudaStream_t St2, int rank_f,  int rank_b, int rank)
{
     if(rank_f==-1 && rank_b==-1) return;
     dim3 block (BLOCK_SIZE_Z, BLOCK_SIZE_Y, 1);
     dim3 grid ((nzt+BLOCK_SIZE_Z-1)/BLOCK_SIZE_Z, (nxt+BLOCK_SIZE_Y-1)/BLOCK_SIZE_Y,1);
     cudaFuncSetCacheConfig(update_boundary_y, cudaFuncCachePreferL1);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d update_boundary, before kernel: %s\n",rank,cudaGetErrorString(err));
     update_boundary_y<<<grid, block, 0, St1>>>(u1, v1, w1, f_u1, f_v1, f_w1, rank_f, Front);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d update_boundary, intermediate: %s\n",rank,cudaGetErrorString(err));
     update_boundary_y<<<grid, block, 0, St2>>>(u1, v1, w1, b_u1, b_v1, b_w1, rank_b, Back);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! update_boundary: %s\n",cudaGetErrorString(err));
     return;
}

extern "C"
void dstrqc_H(float* xx,       float* yy,     float* zz,    float* xy,    float* xz, float* yz,
              float* r1,       float* r2,     float* r3,    float* r4,    float* r5, float* r6,
              float* u1,       float* v1,     float* w1,    float* lam,   float* mu, float* qp,
              float* qs,       float* dcrjx,  float* dcrjy, float* dcrjz, int nyt,   int nzt, 
              cudaStream_t St, float* lam_mu, int NX,       int rankx,    int ranky, int  s_i,  
              int e_i,         int s_j,       int e_j,      int rank)
{
    dim3 block (BLOCK_SIZE_Z, BLOCK_SIZE_Y, 1);
    dim3 grid ((nzt+BLOCK_SIZE_Z-1)/BLOCK_SIZE_Z, (e_j-s_j+1+BLOCK_SIZE_Y-1)/BLOCK_SIZE_Y,1);
    cudaFuncSetCacheConfig(dstrqc, cudaFuncCachePreferL1);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dstrqc, before kernel: %s\n",rank,cudaGetErrorString(err));
    dstrqc<<<grid, block, 0, St>>>(xx,    yy,    zz,  xy,  xz, yz, r1, r2,    r3,    r4,    r5,     r6, 
                                   u1,    v1,    w1,  lam, mu, qp, qs, dcrjx, dcrjy, dcrjz, lam_mu, NX, 
                                   rankx, ranky, s_i, e_i, s_j);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! rank=%d dstrqc: %s\n",rank,cudaGetErrorString(err));
    return;
}

extern "C"
void addsrc_H(int i,      int READ_STEP, int dim,    int* psrc,  int npsrc,  cudaStream_t St,
              int igreen, int nzt, 
              float* d1,     float* u1,  float* v1,  float* w1,
              float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
              float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz)
{
    dim3 grid, block;
    if(npsrc < 256)
    {
       block.x = npsrc;
       grid.x = 1;
    }
    else
    {
       block.x = 256;
       grid.x  = int((npsrc+255)/256);
    }
/*
    int nyt=420, nxt=310;
    int i_=2+4*loop,j=2+4*loop,k=align;
    int pos = i_*(nyt+4+8*loop)*(nzt+2*align)+j*(nzt+2*align)+k;
printf("xx,yy,zz,xy,xz,yz=%f,%f,%f,%f,%f,%f\n",
      xx[pos],yy[pos],zz[pos],xy[pos],xz[pos],yz[pos]);
*/
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! addsrc_H, before addsrc_cu: %s\n",cudaGetErrorString(err));
    addsrc_cu<<<grid, block, 0, St>>>(i,   READ_STEP, dim, npsrc, psrc, igreen, nzt, 
                                      d1,  u1,        v1,  w1,
                                      axx, ayy,       azz, axz,  ayz,    axy,
                                      xx,  yy,        zz,  xy,   yz,     xz);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! addsrc_H: %s\n",cudaGetErrorString(err));
    return;
}

extern "C"
void addkinsrc_H(int i,   int dim,    int* psrc,  int npsrc,  cudaStream_t St, float* mu,
              float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
              float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz,
              float* mom)
{
    dim3 grid, block;
    if(npsrc < 256)
    {
       block.x = npsrc;
       grid.x = 1;
    }
    else
    {
       block.x = 256;
       grid.x  = int((npsrc+255)/256);
    }
    cudaError_t cerr;
    cerr=cudaGetLastError();
    if(cerr!=cudaSuccess) printf("CUDA ERROR: addkinsrc before kernel: %s\n",cudaGetErrorString(cerr));
    //cudaPrintfInit();
	printf("grid.x=%d, block.x=%d\n", grid.x, block.x);
    addkinsrc_cu<<<grid, block, 0, St>>>(i, dim, psrc, npsrc, mu, axx, ayy, azz, axz, ayz, axy,
                                      xx, yy, zz,  xy,   yz,  xz, mom);
    cerr=cudaGetLastError();
    /*cudaPrintfDisplay(stdout, 1);
    cudaPrintfEnd();*/
    if(cerr!=cudaSuccess) printf("CUDA ERROR: addkinsrc after kernel: %s\n",cudaGetErrorString(cerr));
    return;
}

__device__ void compmt(float str, float dip, float rake,
        float *xx, float *yy, float *zz, float *xz, float *yz, float *xy ){

      //angles must be provided in rads

      *yy= -(sinf(dip)*cosf(rake)*sinf(2.*str)+
           sinf(2.*dip)*sinf(rake)*sinf(str)*sinf(str));

      *xy= sinf(dip)*cosf(rake)*cosf(2.*str)+
           0.5*(sinf(2.*dip)*sinf(rake)*sinf(2.*str));

      *yz= (cosf(dip)*cosf(rake)*cosf(str)+
           cosf(2.*dip)*sinf(rake)*sinf(str));

      *xx= sinf(dip)*cosf(rake)*sinf(2.*str)-
           sinf(2.*dip)*sinf(rake)*cosf(str)*cosf(str);

      *xz= (cosf(dip)*cosf(rake)*sinf(str)-
           cosf(2.*dip)*sinf(rake)*cosf(str));

      *zz= sinf(2.*dip)*sinf(rake);
}

__device__ float brune(float freq, float time){
   register float stf, omega;
   omega=freq * M_PI * 2.;
   if (time > 0.)
      stf = powf(omega, 2.) * time * expf(-omega*time);
   else
      stf = 0.;
   return(stf);
}

/* Liu et al. (2006) source time function.  tau = risetime */
__device__ float liu(float tau, float time){
   register float tau1, tau2, CN, stf;

   tau1 = 0.13 * tau;
   tau2 = tau-tau1;

   CN=M_PI / (1.4*M_PI*tau1 + 1.2*tau1 + 0.3*M_PI*tau2);

   if (time < tau1)
      stf = CN*(0.7 - 0.7*cosf(M_PI*time/tau1) + 0.6*sinf(0.5*M_PI*time/tau1));
   else if (time < 2*tau1)
      stf = CN*(1.0 - 0.7*cosf(M_PI*time/tau1) + 0.3*cosf(M_PI*(time-tau1)/tau2));
   else if (time < tau)
      stf = CN*(0.3 + 0.3*cosf(M_PI*(time-tau1) / tau2));
   else
      stf = 0.;

   return(stf);
}

/* GP 2010, revised based on Liu et al. (2006) source time function.  tau = risetime */
__device__ float gp10(float tau, float time){
   register float tau1, tau2, CN, stf;

   tau1 = 0.13 * tau;
   tau2 = tau-tau1;

   CN=M_PI / (1.5 * M_PI*tau1 + 1.2*tau1 + 0.2 * M_PI * tau2);
   if (time < 0.)
      stf = 0.;
   else if (time < tau1)
      stf = CN*(0.7 - 0.7*cosf(M_PI*time/tau1) + 0.6*sinf(0.5*M_PI*time/tau1));
   else if (time < 2*tau1)
      stf = CN*(1.0 - 0.8*cosf(M_PI*time/tau1) + 0.2*cosf(M_PI*(time-tau1)/tau2));
   else if (time < tau)
      stf = CN*(0.2 + 0.2*cosf(M_PI*(time-tau1) / tau2));
   else
      stf = 0.;

   return(stf);
}



__global__ void addkinsrc_cu(int i, int dim,    int* psrc,  int npsrc, float* mu,
                          float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
                          float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz,
                          float* mom)
{

        register float vtst;
        register int idx, idy, idz, j, pos;

        register float atime, freq, stf;
        register float axxt, ayyt, azzt, axzt, ayzt, axyt;
        register int stf_type;
        register float slip, ruptime, risetime, strike, dip, rake, area;

        register int READ_STEP = 2;
        register double *stff[MAXFILT];
        int n;

        j = blockIdx.x*blockDim.x+threadIdx.x;
        if(j >= npsrc) return;

        // For a kinematic source, the moment-rate is computed at run-time from given subfault parameters,
        // which are stored inside the arrays axx...axz 
        stf_type = (int) axx[j*READ_STEP]; // type of source time function.  1=Brune
        slip = ayy[j*READ_STEP];     // total slip
        ruptime = azz[j*READ_STEP];  // rupture time
        risetime = axz[j*READ_STEP]; // rise time

        atime = i*d_DT;

        if (atime > ruptime) {
       area = ayz[j*READ_STEP];   // subfault area

       strike = axx[j*READ_STEP+1] / 180. * M_PI;   // strike angle (given in degrees)
       dip = ayy[j*READ_STEP+1] / 180. * M_PI;      // dip angle
       rake = azz[j*READ_STEP+1] / 180. * M_PI;     // rake

           compmt(strike, dip, rake, &axxt, &ayyt, &azzt, &axzt, &ayzt, &axyt);

       if (stf_type == 1.0f){
          freq = 1./risetime;
          stf = brune(freq, atime - ruptime);
       }
       else if (stf_type == 2.0f)
          stf = liu(risetime, atime - ruptime);
       else if (stf_type == 3.0f)
        stf = gp10(risetime, atime - ruptime);
     else
          stf = 0.;

       vtst = (float)d_DT/(d_DH*d_DH*d_DH);


        idx = psrc[j*dim]   + 1 + 4*loop;
        idy = psrc[j*dim+1] + 1 + 4*loop;
        idz = psrc[j*dim+2] + align - 1;
        pos = idx*d_slice_1 + idy*d_yline_1 + idz;

       /*idx = psrc[j*dim]   + 1 + ngsl;
       idy = psrc[j*dim+1] + 1 + ngsl;
       idz = psrc[j*dim+2] + align - 1;
       pos = idx*d_slice_1 + idy*d_yline_1 + idz;*/

           stf *= slip*area/mu[pos];
           mom[j] += stf * d_DT;

           //cuPrintf("stf: %d %e %e %e %e\n", j, atime, stf, slip, area);
           //cuPrintf("mom: %d %e\n", j, mom[j]);

           stf *= vtst;

           /*if (j == 0)
          cuPrintf("addkinsrc_cu: (%d,%d,%d) (%e, %e,%e,%e,%e,%e,%e)\n", idx, idy, idz, 
             stf, axxt, ayyt, azzt, axzt, ayzt, axyt);
          cuPrintf("addkinsrc_cu: (%d,%d,%d) (%e, %e, %e m^2, %f m)\n", idx, idy, idz, 
             stf, 1./mu[pos], area, slip);*/

       xx[pos] = xx[pos] - stf*axxt;
       yy[pos] = yy[pos] - stf*ayyt;
       zz[pos] = zz[pos] - stf*azzt;
       xz[pos] = xz[pos] - stf*axzt;
       yz[pos] = yz[pos] - stf*ayzt;
       xy[pos] = xy[pos] - stf*axyt;

        }

        return;
}


extern "C"
void ComputeSGT(float* xx,   float* yy,        float* zz,    float* xy,    float* xz,     float* yz,
//                float* sgt1, float* sgt2,      float* sgt3,  float* sgt4,  float* sgt5,   float* sgt6,
                float* sg1,  float* sg2,       float* mu,    int sgt_numsta, int* sgt_sta,float* sgtBuf,     
                cudaStream_t St,  int nzt,     int SGT_BLOCK_SIZE, int SGT_NUMBLOCKS, float* qs, float* d1)

{
    //dim3 block (BLOCK_SIZE_Z, 1, 1);
    //dim3 grid ((sgt_numsta+BLOCK_SIZE_Z-1)/BLOCK_SIZE_Z, (e_j-s_j+1+BLOCK_SIZE_Y-1)/BLOCK_SIZE_Y,1);
    dim3 block (SGT_BLOCK_SIZE, 1, 1);
    dim3 grid (SGT_NUMBLOCKS, 1, 1);
/*
    int nyt=420, nxt=310;
    int i=2+4*loop,j=2+4*loop,k=align;
    int pos = i*(nyt+4+8*loop)*(nzt+2*align)+j*(nzt+2*align)+k;
printf("pos=%d\n",pos);
printf("xx,yy,zz,xy,xz,yz=%f,%f,%f,%f,%f,%f\n",
      xx[pos],yy[pos],zz[pos],xy[pos],xz[pos],yz[pos]);
printf("sg1,sg2,mu=%f,%f,%f\n",
      sg1[pos],sg2[pos],mu[pos]);
printf("sgt_n,sgt_st,sgtBuf=%d,(%d,%d,%d),%f\n",
      sgt_numsta, sgt_sta[0],sgt_sta[1],sgt_sta[2], sgtBuf[0]);
printf("BS,NB=%d,%d\n",
      SGT_BLOCK_SIZE, SGT_NUMBLOCKS);
*/
    size_t cmemfree, cmemtotal;
    cudaFuncSetCacheConfig(ComputeSGT_cu, cudaFuncCachePreferL1);
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! ComputeSGT before kernel: %s\n",cudaGetErrorString(err));
    //        cudaMemGetInfo(&cmemfree, &cmemtotal);
    //        printf("CUDA memory, inside ComputeSGT numsta=%d, free=%ld, total=%ld\n",sgt_numsta,cmemfree,cmemtotal);
    ComputeSGT_cu<<<grid, block, 0, St>>>(xx,    yy,    zz,   xy,   xz,   yz, 
//                                          sgt1,  sgt2,  sgt3, sgt4, sgt5, sgt6,
                                          sg1,   sg2,   mu,   sgt_numsta, sgt_sta, sgtBuf,
                                          SGT_BLOCK_SIZE, SGT_NUMBLOCKS, qs, d1);
    //        cudaMemGetInfo(&cmemfree, &cmemtotal);
    //        printf("CUDA memory, inside ComputeSGT, numsta=%d, free=%ld, total=%ld\n",sgt_numsta,cmemfree,cmemtotal);
    err = cudaGetLastError();
    if(err != cudaSuccess) printf("CUDA KERNEL ERROR! ComputeSGT: %s\n",cudaGetErrorString(err));
    return;
}


__global__ void dvelcx(float* u1,    float* v1,    float* w1,    float* xx, float* yy, float* zz, float* xy, float* xz, float* yz, 
                      float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, int s_i,   int e_i)
{
    register int   i, j, k, pos,     pos_im1, pos_im2;
    register int   pos_km2, pos_km1, pos_kp1, pos_kp2;
    register int   pos_jm2, pos_jm1, pos_jp1, pos_jp2;
    register int   pos_ip1, pos_jk1, pos_ik1, pos_ijk;
    register float f_xx,    xx_im1,  xx_ip1,  xx_im2;
    register float f_xy,    xy_ip1,  xy_ip2,  xy_im1;
    register float f_xz,    xz_ip1,  xz_ip2,  xz_im1;
    register float f_d1,    f_d2,    f_d3,    f_dcrj, f_dcrjy, f_dcrjz, f_yz;

    k    = blockIdx.x*BLOCK_SIZE_Z+threadIdx.x+align;
    j    = blockIdx.y*BLOCK_SIZE_Y+threadIdx.y+2+4*loop;
    i    = e_i;
    pos  = i*d_slice_1+j*d_yline_1+k;

    f_xx    = xx[pos+d_slice_1];
    xx_im1  = xx[pos];
    xx_im2  = xx[pos-d_slice_1]; 
    xy_ip1  = xy[pos+d_slice_2];
    f_xy    = xy[pos+d_slice_1];
    xy_im1  = xy[pos];
    xz_ip1  = xz[pos+d_slice_2];
    f_xz    = xz[pos+d_slice_1];
    xz_im1  = xz[pos];
    f_dcrjz = dcrjz[k];
    f_dcrjy = dcrjy[j]; 
    for(i=e_i;i>=s_i;i--)   
    {
        pos_km2  = pos-2;
        pos_km1  = pos-1;
        pos_kp1  = pos+1;
        pos_kp2  = pos+2;
        pos_jm2  = pos-d_yline_2;
        pos_jm1  = pos-d_yline_1;
        pos_jp1  = pos+d_yline_1;
        pos_jp2  = pos+d_yline_2;
        pos_im1  = pos-d_slice_1;
        pos_im2  = pos-d_slice_2;
        pos_ip1  = pos+d_slice_1;
        pos_jk1  = pos-d_yline_1-1;
        pos_ik1  = pos+d_slice_1-1;
        pos_ijk  = pos+d_slice_1-d_yline_1;

        xx_ip1   = f_xx;
        f_xx     = xx_im1;
        xx_im1   = xx_im2;
        xx_im2   = xx[pos_im2];
        xy_ip2   = xy_ip1;
        xy_ip1   = f_xy;
        f_xy     = xy_im1;
        xy_im1   = xy[pos_im1];
        xz_ip2   = xz_ip1;
        xz_ip1   = f_xz;
        f_xz     = xz_im1;
        xz_im1   = xz[pos_im1];
        f_yz     = yz[pos];

        f_dcrj   = dcrjx[i]*f_dcrjy*f_dcrjz;
        f_d1     = 0.25*(d_1[pos] + d_1[pos_jm1] + d_1[pos_km1] + d_1[pos_jk1]);
        f_d2     = 0.25*(d_1[pos] + d_1[pos_ip1] + d_1[pos_km1] + d_1[pos_ik1]);
        f_d3     = 0.25*(d_1[pos] + d_1[pos_ip1] + d_1[pos_jm1] + d_1[pos_ijk]);

        f_d1     = d_dth/f_d1;
        f_d2     = d_dth/f_d2;
	f_d3     = d_dth/f_d3;

    	u1[pos]  = (u1[pos] + f_d1*( d_c1*(f_xx        - xx_im1)      + d_c2*(xx_ip1      - xx_im2) 
                                   + d_c1*(f_xy        - xy[pos_jm1]) + d_c2*(xy[pos_jp1] - xy[pos_jm2])
                                   + d_c1*(f_xz        - xz[pos_km1]) + d_c2*(xz[pos_kp1] - xz[pos_km2]) ))*f_dcrj; 
        v1[pos]  = (v1[pos] + f_d2*( d_c1*(xy_ip1      - f_xy)        + d_c2*(xy_ip2      - xy_im1)
                                   + d_c1*(yy[pos_jp1] - yy[pos])     + d_c2*(yy[pos_jp2] - yy[pos_jm1])
                                   + d_c1*(f_yz        - yz[pos_km1]) + d_c2*(yz[pos_kp1] - yz[pos_km2]) ))*f_dcrj;

        w1[pos]  = (w1[pos] + f_d3*( d_c1*(xz_ip1      - f_xz)        + d_c2*(xz_ip2      - xz_im1)
                                   + d_c1*(f_yz        - yz[pos_jm1]) + d_c2*(yz[pos_jp1] - yz[pos_jm2])
                                   + d_c1*(zz[pos_kp1] - zz[pos])     + d_c2*(zz[pos_kp2] - zz[pos_km1]) ))*f_dcrj;
        pos      = pos_im1;
    }

    return;
}


__global__ void dvelcy(float* u1,    float* v1,    float* w1,    float* xx,  float* yy,   float* zz,   float* xy, float* xz, float* yz,
                       float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, float* s_u1, float* s_v1, float* s_w1, int s_j,   int e_j)
{
    register int   i, j, k, pos,     j2,      pos2, pos_jm1, pos_jm2;
    register int   pos_km2, pos_km1, pos_kp1, pos_kp2;
    register int   pos_im2, pos_im1, pos_ip1, pos_ip2;
    register int   pos_jk1, pos_ik1, pos_ijk;
    register float f_xy,    xy_jp1,  xy_jm1,  xy_jm2;
    register float f_yy,    yy_jp2,  yy_jp1,  yy_jm1;
    register float f_yz,    yz_jp1,  yz_jm1,  yz_jm2;
    register float f_d1,    f_d2,    f_d3,    f_dcrj, f_dcrjx, f_dcrjz, f_xz;

    k     = blockIdx.x*BLOCK_SIZE_Z+threadIdx.x+align;
    i     = blockIdx.y*BLOCK_SIZE_Y+threadIdx.y+2+4*loop;
    j     = e_j;
    j2    = 4*loop-1;
    pos   = i*d_slice_1+j*d_yline_1+k;
    pos2  = i*4*loop*d_yline_1+j2*d_yline_1+k; 

    f_xy    = xy[pos+d_yline_1];
    xy_jm1  = xy[pos];
    xy_jm2  = xy[pos-d_yline_1];
    yy_jp1  = yy[pos+d_yline_2];
    f_yy    = yy[pos+d_yline_1];
    yy_jm1  = yy[pos];
    f_yz    = yz[pos+d_yline_1];
    yz_jm1  = yz[pos];
    yz_jm2  = yz[pos-d_yline_1];
    f_dcrjz = dcrjz[k];
    f_dcrjx = dcrjx[i];
    for(j=e_j; j>=s_j; j--)
    {
        pos_km2  = pos-2;
        pos_km1  = pos-1;
        pos_kp1  = pos+1;
        pos_kp2  = pos+2;
        pos_jm2  = pos-d_yline_2;
        pos_jm1  = pos-d_yline_1;
        pos_im1  = pos-d_slice_1;
        pos_im2  = pos-d_slice_2;
        pos_ip1  = pos+d_slice_1;
        pos_ip2  = pos+d_slice_2;
        pos_jk1  = pos-d_yline_1-1;
        pos_ik1  = pos+d_slice_1-1;
        pos_ijk  = pos+d_slice_1-d_yline_1;

        xy_jp1   = f_xy;
        f_xy     = xy_jm1;
        xy_jm1   = xy_jm2;
        xy_jm2   = xy[pos_jm2];
        yy_jp2   = yy_jp1;
        yy_jp1   = f_yy;
        f_yy     = yy_jm1;
        yy_jm1   = yy[pos_jm1];
        yz_jp1   = f_yz;
        f_yz     = yz_jm1;
        yz_jm1   = yz_jm2;
        yz_jm2   = yz[pos_jm2];
        f_xz     = xz[pos];

        f_dcrj   = f_dcrjx*dcrjy[j]*f_dcrjz;
        f_d1     = 0.25*(d_1[pos] + d_1[pos_jm1] + d_1[pos_km1] + d_1[pos_jk1]);
        f_d2     = 0.25*(d_1[pos] + d_1[pos_ip1] + d_1[pos_km1] + d_1[pos_ik1]);
        f_d3     = 0.25*(d_1[pos] + d_1[pos_ip1] + d_1[pos_jm1] + d_1[pos_ijk]);

        f_d1     = d_dth/f_d1;
        f_d2     = d_dth/f_d2;
        f_d3     = d_dth/f_d3;

        s_u1[pos2] = (u1[pos] + f_d1*( d_c1*(xx[pos]     - xx[pos_im1]) + d_c2*(xx[pos_ip1] - xx[pos_im2])
                                     + d_c1*(f_xy        - xy_jm1)      + d_c2*(xy_jp1      - xy_jm2)
                                     + d_c1*(f_xz        - xz[pos_km1]) + d_c2*(xz[pos_kp1] - xz[pos_km2]) ))*f_dcrj;
        s_v1[pos2] = (v1[pos] + f_d2*( d_c1*(xy[pos_ip1] - f_xy)        + d_c2*(xy[pos_ip2] - xy[pos_im1])
                                     + d_c1*(yy_jp1      - f_yy)        + d_c2*(yy_jp2      - yy_jm1)
                                     + d_c1*(f_yz        - yz[pos_km1]) + d_c2*(yz[pos_kp1] - yz[pos_km2]) ))*f_dcrj;
        s_w1[pos2] = (w1[pos] + f_d3*( d_c1*(xz[pos_ip1] - f_xz)        + d_c2*(xz[pos_ip2] - xz[pos_im1])
                                     + d_c1*(f_yz        - yz_jm1)      + d_c2*(yz_jp1      - yz_jm2)
                                     + d_c1*(zz[pos_kp1] - zz[pos])     + d_c2*(zz[pos_kp2] - zz[pos_km1]) ))*f_dcrj;

        pos        = pos_jm1;
        pos2       = pos2 - d_yline_1;
    }
    return;
}

__global__ void update_boundary_y(float* u1, float* v1, float* w1, float* s_u1, float* s_v1, float* s_w1, int rank, int flag)
{
    register int i, j, k, pos, posj;
    k     = blockIdx.x*BLOCK_SIZE_Z+threadIdx.x+align;
    i     = blockIdx.y*BLOCK_SIZE_Y+threadIdx.y+2+4*loop;

    if(flag==Front && rank!=-1){
	j     = 2;
    	pos   = i*d_slice_1+j*d_yline_1+k;
        posj  = i*4*loop*d_yline_1+k;
	for(j=2;j<2+4*loop;j++){
		u1[pos] = s_u1[posj];
		v1[pos] = s_v1[posj];
		w1[pos] = s_w1[posj];
		pos	= pos  + d_yline_1;
  		posj	= posj + d_yline_1;	
	}
    }

    if(flag==Back && rank!=-1){
    	j     = d_nyt+4*loop+2;
    	pos   = i*d_slice_1+j*d_yline_1+k;
        posj  = i*4*loop*d_yline_1+k;
	for(j=d_nyt+4*loop+2;j<d_nyt+8*loop+2;j++){
	        u1[pos] = s_u1[posj];
                v1[pos] = s_v1[posj];
                w1[pos] = s_w1[posj];
                pos     = pos  + d_yline_1;
                posj    = posj + d_yline_1;
	}
    }
    return;
}

__global__ void ComputeSGT_cu(float* xx,   float* yy,    float* zz,    float* xy,    float* xz,     float* yz,
//                              float* sgt1, float* sgt2,  float* sgt3,  float* sgt4,  float* sgt5,   float* sgt6,
                              float* sg1,  float* sg2,   float* mu,    int sgt_numsta, int* sgt_sta,float* sgtBuf,
                              int SGT_BLOCK_SIZE,        int SGT_NUMBLOCKS, float* qs, float* d1)
{
    register int ind, i, j, k, pos;

    //if(SGT_NUMBLOCKS==1)
      //ind = threadIdx.x;
    //else
      ind = blockIdx.x*SGT_BLOCK_SIZE + threadIdx.x;
    if(sgt_numsta <= ind) return;
    i = sgt_sta[ind*3]   + 1 + 4*loop;
    j = sgt_sta[ind*3+1] + 1 + 4*loop;
    k = sgt_sta[ind*3+2] + align - 1;
    pos  = i*d_slice_1+j*d_yline_1+k;
    ind  = ind*6;
//    for(i=e_i;i>=s_i;i--)
//    {
/*
        (*sgtBuf)[ind] = sg1[pos]*xx[pos] + sg2[pos]*yy[pos] + sg2[pos]*zz[pos];
        ind++;
        (*sgtBuf)[ind] = sg2[pos]*xx[pos] + sg1[pos]*yy[pos] + sg2[pos]*zz[pos];
        ind++;
        (*sgtBuf)[ind] = sg2[pos]*xx[pos] + sg2[pos]*yy[pos] + sg1[pos]*zz[pos];
        ind++;
        (*sgtBuf)[ind] = mu[pos]*xy[pos];
        ind++;
        (*sgtBuf)[ind] = mu[pos]*xz[pos];
        ind++;
        (*sgtBuf)[ind] = mu[pos]*yz[pos];
*/
		//Recompute sg1, sg2, mu to use uncorrected values
		float FP = 1.0;
		float FL = 0.01;
		float FH = 25.0;
		float pi = 3.14159;
		float w0 = 2*pi*FP;
		float ww1 = 2*pi*FL;
		float w2 = 2*pi*FH;
		float taumax = 1.0/ww1;
		float taumin = 1.0/w2;
		float tmp1 = 2.0/pi*log(taumax/taumin);
		float tmp2 = 2.0/pi*log(w0*taumin);
		//Assumes Qp = 2Qs
		float orig_qs = tmp1/qs[pos]+tmp2;
		float orig_qp = 2*orig_qs;
		float vs_scalefac = 1 + log(w2/w0)/(pi*orig_qs);
		float vp_scalefac = 1 + log(w2/w0)/(pi*orig_qp);
		float orig_vs = sqrt(1.0/(mu[pos]*d1[pos]))/vs_scalefac;
		float lam = 0.5*(((-1.0*mu[pos])/(2.0*sg2[pos])-3.0)*mu[pos]);
		float orig_vp = sqrt(1.0/(d1[pos]*lam)+2*vs_scalefac*orig_vs*vs_scalefac*orig_vs)/vp_scalefac;
		float local_mu = mu[pos]*vs_scalefac*vs_scalefac;
		float local_lam = 1.0/(d1[pos]*(orig_vp*orig_vp-2.0*orig_vs*orig_vs));
		float local_sg1 = (1.0/local_lam+1.0/local_mu)/(1.0/local_mu*(3.0/local_lam+2.0/local_mu));
		float local_sg2 = -1.0/local_lam/(2.0/local_mu*(3.0/local_lam+2.0/local_mu));


        /*sgtBuf[ind] = sg1[pos]*xx[pos] + sg2[pos]*yy[pos] + sg2[pos]*zz[pos];
        ind++;
        sgtBuf[ind] = sg2[pos]*xx[pos] + sg1[pos]*yy[pos] + sg2[pos]*zz[pos];
        ind++;
        sgtBuf[ind] = sg2[pos]*xx[pos] + sg2[pos]*yy[pos] + sg1[pos]*zz[pos];
        ind++;
        sgtBuf[ind] = mu[pos]*xy[pos];
        ind++;
        sgtBuf[ind] = mu[pos]*xz[pos];
        ind++;
        sgtBuf[ind] = mu[pos]*yz[pos];*/

		sgtBuf[ind] = local_sg1*xx[pos] + local_sg2*yy[pos] + local_sg2*zz[pos];
        ind++;
        sgtBuf[ind] = local_sg2*xx[pos] + local_sg1*yy[pos] + local_sg2*zz[pos];
        ind++;
        sgtBuf[ind] = local_sg2*xx[pos] + local_sg2*yy[pos] + local_sg1*zz[pos];
        ind++;
        sgtBuf[ind] = local_mu*xy[pos];
        ind++;
        sgtBuf[ind] = local_mu*xz[pos];
        ind++;
        sgtBuf[ind] = local_mu*yz[pos];
		
		//To output strains
        /*sgtBuf[ind] = xx[pos];
        ind++;
        sgtBuf[ind] = yy[pos];
        ind++;
        sgtBuf[ind] = zz[pos];
        ind++;
        sgtBuf[ind] = xy[pos];
        ind++;
        sgtBuf[ind] = xz[pos];
        ind++;
        sgtBuf[ind] = yz[pos];*/
//        pos       = pos - d_slice_1;
//    }
    return;
}


__global__ void dstrqc(float* xx, float* yy,    float* zz,    float* xy,    float* xz,     float* yz,
                       float* r1, float* r2,    float* r3,    float* r4,    float* r5,     float* r6,
                       float* u1, float* v1,    float* w1,    float* lam,   float* mu,     float* qp,
                       float* qs, float* dcrjx, float* dcrjy, float* dcrjz, float* lam_mu, int NX,    
                       int rankx, int ranky,    int s_i,      int e_i,      int s_j)
{
    register int   i,  j,  k,  g_i;
    register int   pos,     pos_ip1, pos_im2, pos_im1;
    register int   pos_km2, pos_km1, pos_kp1, pos_kp2;
    register int   pos_jm2, pos_jm1, pos_jp1, pos_jp2;
    register int   pos_ik1, pos_jk1, pos_ijk, pos_ijk1;
    register float vs1, vs2, vs3, a1, tmp, vx1;
    register float xl,  xm,  xmu1, xmu2, xmu3;
    register float qpa, h,   h1,   h2,   h3;
    register float f_vx1, f_vx2,  f_dcrj, f_r,  f_dcrjy, f_dcrjz;
    register float f_u1, u1_ip1, u1_ip2, u1_im1;
    register float f_v1, v1_im1, v1_ip1, v1_im2;
    register float f_w1, w1_im1, w1_im2, w1_ip1;
    
    k    = blockIdx.x*BLOCK_SIZE_Z+threadIdx.x+align;
    j    = blockIdx.y*BLOCK_SIZE_Y+threadIdx.y+s_j;
    i    = e_i;
    pos  = i*d_slice_1+j*d_yline_1+k;

    u1_ip1 = u1[pos+d_slice_2];
    f_u1   = u1[pos+d_slice_1];
    u1_im1 = u1[pos];    
    f_v1   = v1[pos+d_slice_1];
    v1_im1 = v1[pos];
    v1_im2 = v1[pos-d_slice_1];
    f_w1   = w1[pos+d_slice_1];
    w1_im1 = w1[pos];
    w1_im2 = w1[pos-d_slice_1];
    f_dcrjz = dcrjz[k];
    f_dcrjy = dcrjy[j];
    for(i=e_i;i>=s_i;i--)
    {
        f_vx1    = tex1Dfetch(p_vx1, pos);
        f_vx2    = tex1Dfetch(p_vx2, pos);
        f_dcrj   = dcrjx[i]*f_dcrjy*f_dcrjz;

        pos_km2  = pos-2;
        pos_km1  = pos-1;
        pos_kp1  = pos+1;
        pos_kp2  = pos+2;
        pos_jm2  = pos-d_yline_2;
        pos_jm1  = pos-d_yline_1;
        pos_jp1  = pos+d_yline_1;
        pos_jp2  = pos+d_yline_2;
        pos_im2  = pos-d_slice_2;
        pos_im1  = pos-d_slice_1;
        pos_ip1  = pos+d_slice_1;
        pos_jk1  = pos-d_yline_1-1;
        pos_ik1  = pos+d_slice_1-1;
        pos_ijk  = pos+d_slice_1-d_yline_1;
        pos_ijk1 = pos+d_slice_1-d_yline_1-1;

        xl       = 8.0/(  lam[pos]      + lam[pos_ip1] + lam[pos_jm1] + lam[pos_ijk]
                        + lam[pos_km1]  + lam[pos_ik1] + lam[pos_jk1] + lam[pos_ijk1] );
        xm       = 16.0/( mu[pos]       + mu[pos_ip1]  + mu[pos_jm1]  + mu[pos_ijk]
                        + mu[pos_km1]   + mu[pos_ik1]  + mu[pos_jk1]  + mu[pos_ijk1] );
        xmu1     = 2.0/(  mu[pos]       + mu[pos_km1] );
        xmu2     = 2.0/(  mu[pos]       + mu[pos_jm1] );
        xmu3     = 2.0/(  mu[pos]       + mu[pos_ip1] );
        xl       = xl  +  xm;

        qpa      = 0.0625*( qp[pos]     + qp[pos_ip1] + qp[pos_jm1] + qp[pos_ijk]
                          + qp[pos_km1] + qp[pos_ik1] + qp[pos_jk1] + qp[pos_ijk1] );
		//qpa		= 0.0625*(qp[pos] + qp[pos_ip1] + qp[pos_jm1] + qp[pos_ijk]
		//			+ qp[pos_kp1] + qp[pos_ip1+1] + qp[pos_jm1+1] + qp[pos_ijk+1]);
        h        = 0.0625*( qs[pos]     + qs[pos_ip1] + qs[pos_jm1] + qs[pos_ijk]
                          + qs[pos_km1] + qs[pos_ik1] + qs[pos_jk1] + qs[pos_ijk1] );
		//h		= 0.0625*(qs[pos] + qs[pos_ip1] + qs[pos_jm1] + qs[pos_ijk]
		//				+ qs[pos_kp1] + qs[pos_ip1+1] + qs[pos_jm1+1] + qs[pos_ijk+1]);
        h1       = 0.250*(  qs[pos]     + qs[pos_km1] );
		//h1		= 0.250*(qs[pos] + qs[pos_kp1]);
        h2       = 0.250*(  qs[pos]     + qs[pos_jm1] );
        h3       = 0.250*(  qs[pos]     + qs[pos_ip1] );

        h        = -xm*h*d_dh1;
        h1       = -xmu1*h1*d_dh1;
        h2       = -xmu2*h2*d_dh1;
        h3       = -xmu3*h3*d_dh1;
        qpa      = -qpa*xl*d_dh1;
        xm       = xm*d_dth;
        xmu1     = xmu1*d_dth;
        xmu2     = xmu2*d_dth;
        xmu3     = xmu3*d_dth;
        xl       = xl*d_dth;
        f_vx2    = f_vx2*f_vx1;
        h        = h*f_vx1;
        h1       = h1*f_vx1;
        h2       = h2*f_vx1;
        h3       = h3*f_vx1;
        qpa      = qpa*f_vx1;

        xm       = xm+d_DT*h;
        xmu1     = xmu1+d_DT*h1;
        xmu2     = xmu2+d_DT*h2;
        xmu3     = xmu3+d_DT*h3;
        vx1      = d_DT*(1+f_vx2);
        
        u1_ip2   = u1_ip1;
        u1_ip1   = f_u1;
        f_u1     = u1_im1;
        u1_im1   = u1[pos_im1];
        v1_ip1   = f_v1;
        f_v1     = v1_im1;
        v1_im1   = v1_im2;
        v1_im2   = v1[pos_im2];
        w1_ip1   = f_w1;
        f_w1     = w1_im1;
        w1_im1   = w1_im2;
        w1_im2   = w1[pos_im2];

        if(k == d_nzt+align-1)
        {
		u1[pos_kp1] = f_u1 - (f_w1        - w1_im1);
    		v1[pos_kp1] = f_v1 - (w1[pos_jp1] - f_w1);

                g_i  = d_nxt*rankx + i - 4*loop - 1;

                //Interpretation: do not consider gradients over w1 going over domain boundary
                //for u1 and v1, set value that would be outside of boundary to zero
 
                /*original implementation*/
    		/*if(g_i<NX)
        		vs1	= u1_ip1 - (w1_ip1    - f_w1);
    		else
        		vs1	= 0.0;

                g_i  = d_nyt*ranky + j - 4*loop - 1;
    		if(g_i>1)
        		vs2	= v1[pos_jm1] - (f_w1 - w1[pos_jm1]);
    		else
        		vs2	= 0.0;

    		w1[pos_kp1]	= w1[pos_km1] - lam_mu[i*(d_nyt+4+8*loop) + j]*((vs1         - u1[pos_kp1]) + (u1_ip1 - f_u1)
                                      +     			                (v1[pos_kp1] - vs2)         + (f_v1   - v1[pos_jm1]) );*/

//              w1[i,j,k+1] =
//                w1[i,j,k-1] - lambda/(lambda + 2 mu) *
//                u1[i+1,j,k] - u1[i,j,k+1] + v1[i,j,k+1] - u1[i,j,k] - w1[i+1,j,k] + w1[i,j,k] + u1[i+1,j,k]
//		  + w1[i,j,k] - w1[i,j-1,k] + v1(i,j,k)

                /* new implementation */

    		/*if(g_i<NX)
        		vs1	= u1_ip1 - (w1_ip1    - f_w1);
    		else
        		vs1	= 0.0;

                g_i  = d_nyt*ranky + j - 4*loop - 1;
    		if(g_i>1)
        		vs2	= v1[pos_jm1] - (f_w1 - w1[pos_jm1]);
    		else
        		vs2	= 0.0;*/

                float al, amu, a, b, a5, b5;
                register int pos_jm1_kp1, pos_ip1_kp1;
                pos_jm1_kp1 = pos_jm1 + 1;
                pos_ip1_kp1 = pos_ip1 + 1;

                al=1./(0.5*(lam[pos] + lam[pos_km1]));
                amu=1./(0.5*(mu[pos] + mu[pos_km1]));

                a=al;
                b=al+2*amu;

                al=1./(0.5*(lam[pos] + lam[pos_km1]));
                amu=1./(0.5*(mu[pos] + mu[pos_km1]));

                a5=al;
                b5=al+2*amu;

                //may still need workaround for i,j close to boundary
    		w1[pos_kp1] =  w1[pos] - a/b *(u1[pos_ip1_kp1] - u1[pos_kp1] + v1[pos_kp1] - v1[pos_jm1_kp1]) 
                              -(w1[pos] - w1[pos_km1])
                              -a5/b5 * (u1[pos_ip1] - u1[pos] + v1[pos] -v1[pos_jm1]);

        }
	else if(k == d_nzt+align-2)
	{
                u1[pos_kp2] = u1[pos_kp1] - (w1[pos_kp1]   - w1[pos_im1+1]);
                v1[pos_kp2] = v1[pos_kp1] - (w1[pos_jp1+1] - w1[pos_kp1]);
	}
 
    	vs1      = d_c1*(u1_ip1 - f_u1)        + d_c2*(u1_ip2      - u1_im1);
        vs2      = d_c1*(f_v1   - v1[pos_jm1]) + d_c2*(v1[pos_jp1] - v1[pos_jm2]);
        vs3      = d_c1*(f_w1   - w1[pos_km1]) + d_c2*(w1[pos_kp1] - w1[pos_km2]);
 
        tmp      = xl*(vs1+vs2+vs3);
        a1       = qpa*(vs1+vs2+vs3);
        tmp      = tmp+d_DT*a1;

        f_r      = r1[pos];
        xx[pos]  = (xx[pos]  + tmp - xm*(vs2+vs3) + vx1*f_r)*f_dcrj;
        r1[pos]  = f_vx2*f_r - h*(vs2+vs3)        + a1;
        f_r      = r2[pos];
        yy[pos]  = (yy[pos]  + tmp - xm*(vs1+vs3) + vx1*f_r)*f_dcrj;
        r2[pos]  = f_vx2*f_r - h*(vs1+vs3)        + a1;
        f_r      = r3[pos];
        zz[pos]  = (zz[pos]  + tmp - xm*(vs1+vs2) + vx1*f_r)*f_dcrj;
        r3[pos]  = f_vx2*f_r - h*(vs1+vs2)        + a1;

        vs1      = d_c1*(u1[pos_jp1] - f_u1)   + d_c2*(u1[pos_jp2] - u1[pos_jm1]);
        vs2      = d_c1*(f_v1        - v1_im1) + d_c2*(v1_ip1      - v1_im2);
        f_r      = r4[pos];
        xy[pos]  = (xy[pos]  + xmu1*(vs1+vs2) + vx1*f_r)*f_dcrj;
        r4[pos]  = f_vx2*f_r + h1*(vs1+vs2);
  
        if(k == d_nzt+align-1)
        {
                zz[pos+1] = -zz[pos];
        	xz[pos]   = 0.0;
                yz[pos]   = 0.0;
        }
        else
        {
        	vs1     = d_c1*(u1[pos_kp1] - f_u1)   + d_c2*(u1[pos_kp2] - u1[pos_km1]);
        	vs2     = d_c1*(f_w1        - w1_im1) + d_c2*(w1_ip1      - w1_im2);
        	f_r     = r5[pos];
        	xz[pos] = (xz[pos]  + xmu2*(vs1+vs2) + vx1*f_r)*f_dcrj;
        	r5[pos] = f_vx2*f_r + h2*(vs1+vs2);
	 

        	vs1     = d_c1*(v1[pos_kp1] - f_v1) + d_c2*(v1[pos_kp2] - v1[pos_km1]);
        	vs2     = d_c1*(w1[pos_jp1] - f_w1) + d_c2*(w1[pos_jp2] - w1[pos_jm1]);
        	f_r     = r6[pos];
        	yz[pos] = (yz[pos]  + xmu3*(vs1+vs2) + vx1*f_r)*f_dcrj;
        	r6[pos] = f_vx2*f_r + h3*(vs1+vs2);

                if(k == d_nzt+align-2)
                {
                    zz[pos+3] = -zz[pos];
                    xz[pos+2] = -xz[pos];
                    yz[pos+2] = -yz[pos];                                               
		}
		else if(k == d_nzt+align-3)
		{
                    xz[pos+4] = -xz[pos];
                    yz[pos+4] = -yz[pos];
		}
 	}
        pos     = pos_im1;
    }
    return;
}


__global__ void addsrc_cu(int i,      int READ_STEP, int dim,    int npsrc, int* psrc,  int igreen, int nzt,
                          float* d1,  float* u1,     float* v1,  float* w1,
                          float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
                          float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz)
{
        register float vtst, vtst1, tmp;
        register int idx, idy, idz, j, pos, pos_nzt;
        j = blockIdx.x*blockDim.x+threadIdx.x;
        if(j >= npsrc) return;
        vtst = d_DT/(d_DH*d_DH*d_DH);
        vtst1 = 1.0/(d_DH*d_DH);

        i   = i - 1;
        idx = psrc[j*dim]   + 1 + 4*loop;
        idy = psrc[j*dim+1] + 1 + 4*loop;
        idz = psrc[j*dim+2] + align - 1;
        pos = idx*d_slice_1 + idy*d_yline_1 + idz;
        pos_nzt = pos - idz + nzt+align-1;

        if(igreen == -1)
        {
          xx[pos] = xx[pos] - vtst*axx[j*READ_STEP+i];
          yy[pos] = yy[pos] - vtst*ayy[j*READ_STEP+i];
          zz[pos] = zz[pos] - vtst*azz[j*READ_STEP+i];
          xz[pos] = xz[pos] - vtst*axz[j*READ_STEP+i];
          yz[pos] = yz[pos] - vtst*ayz[j*READ_STEP+i];
          xy[pos] = xy[pos] - vtst*axy[j*READ_STEP+i];
        }
        else if(igreen == 1)
        {
          u1[pos] = u1[pos] + vtst*axx[j*READ_STEP+i]/d1[pos];
        }
        else if(igreen == 2)
        {
          v1[pos] = v1[pos] + vtst*ayy[j*READ_STEP+i]/d1[pos];
        }
        else if(igreen == 3)
        {
          w1[pos] = w1[pos] + vtst*azz[j*READ_STEP+i]/d1[pos];
        }
        else if(igreen == -2)
        {
          u1[pos] = u1[pos] + vtst*axx[j*READ_STEP+i]/d1[pos];
          v1[pos] = v1[pos] + vtst*ayy[j*READ_STEP+i]/d1[pos];
          w1[pos] = w1[pos] + vtst*azz[j*READ_STEP+i]/d1[pos];
        }
        else if(igreen == 4){
          tmp = vtst1*axx[j*READ_STEP+i];
          xz[pos_nzt] = tmp;
          tmp = tmp*2;
          xz[pos_nzt+1] = tmp - xz[pos_nzt-1];
          xz[pos_nzt+2] = tmp - xz[pos_nzt-2];
        }
        else if(igreen == 5){
          tmp = vtst1*ayy[j*READ_STEP+i];
          yz[pos_nzt] = tmp;
          tmp = tmp*2;
          yz[pos_nzt+1] = tmp - yz[pos_nzt-1];
          yz[pos_nzt+2] = tmp - yz[pos_nzt-2];
        }
        else if(igreen == 6){
          tmp = 2.0*vtst1*azz[j*READ_STEP+i];
          zz[pos_nzt+1] = tmp - zz[pos_nzt];
          zz[pos_nzt+2] = tmp - zz[pos_nzt-1];
        }


        return;
}
