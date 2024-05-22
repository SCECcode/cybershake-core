/*  
********************************************************************************
* pmcl3d.c                                                                     *
* programming in C&CUDA language                                                    *
* Author: Jun Zhou                                                             * 
* First Version: Cerjan Mode and Homogenous                                    *
********************************************************************************
*/  
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pmcl3d.h"
#include "cuda_to_hip.h"
#include <hip/hip_runtime.h>
const double   micro = 1.0e-6;

void SetDeviceConstValue(float DH, float DT, int nxt, int nyt, int nzt);
void BindArrayToTexture(float* vx1, float* vx2, int memsize);
void UnBindArrayFromTexture();
void dvelcx_H(float* u1,    float* v1,    float* w1,    float* xx,  float* yy, float* zz, float* xy,       float* xz, float* yz,
              float* dcrjx, float* dcrjy, float* dcrjz, float* d_1, int nyt,   int nzt,   cudaStream_t St, int s_i,   int e_i, int rank);
void dvelcy_H(float* u1,       float* v1,    float* w1,    float* xx,  float* yy, float* zz, float* xy,   float* xz,   float* yz,
              float* dcrjx,    float* dcrjy, float* dcrjz, float* d_1, int nxt,   int nzt,   float* s_u1, float* s_v1, float* s_w1,
              cudaStream_t St, int s_j,      int e_j,      int rank,   int rank_me);
void dstrqc_H(float* xx,       float* yy,     float* zz,    float* xy,    float* xz, float* yz,
              float* r1,       float* r2,     float* r3,    float* r4,    float* r5, float* r6,
              float* u1,       float* v1,     float* w1,    float* lam,   float* mu, float* qp,
              float* qs,       float* dcrjx,  float* dcrjy, float* dcrjz, int nyt,   int nzt, 
              cudaStream_t St, float* lam_mu, int NX,       int rankx,    int ranky, int s_i,  
              int e_i,         int s_j,       int e_j,      int rank, float* d_vx1, float* d_vx2);
void addsrc_H(int i,      int READ_STEP, int dim,    int* psrc,  int npsrc,  cudaStream_t St,
              int igreen, int nzt, 
              float* d1,     float* u1,  float* v1,  float* w1,
              float* axx, float* ayy,    float* azz, float* axz, float* ayz, float* axy,
              float* xx,  float* yy,     float* zz,  float* xy,  float* yz,  float* xz);
void ComputeSGT(float* xx,   float* yy,        float* zz,    float* xy,    float* xz,     float* yz,
//                float* sgt1, float* sgt2,      float* sgt3,  float* sgt4,  float* sgt5,   float* sgt6,
                float* sg1,  float* sg2,       float* mu,    int sgt_numsta, int* sgt_sta,float* sgtBuf,     
                cudaStream_t St,  int nzt,     int SGT_BLOCK_SIZE, int SGT_NUMBLOCKS, float* qs, float* d1);

void calcRecordingPoints(int *rec_nbgx, int *rec_nedx, 
  int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz, 
  int *rec_nxt, int *rec_nyt, int *rec_nzt, MPI_Offset *displacement,
  long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
  int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY, 
  int NBGZ, int NEDZ, int NSKPZ, int *coord);

double gethrtime()
{
    struct timeval TV;
    int RC = gettimeofday(&TV,NULL);

    if (RC == -1){
       printf("Bad call to gettimeofday\n");
       return(-1);
    }

    return ( ((double)TV.tv_sec ) + micro * ((double)  TV.tv_usec));
}

int main(int argc,char **argv)
{
//  variable definition begins
    float TMAX, DH, DT, ARBC, PHT;
    int   NPC, ND, NSRC, NST;
    int   NVE, NVAR, MEDIASTART, IFAULT, READ_STEP, READ_STEP_GPU;
    int   NX, NY, NZ, PX, PY, IDYNA, SoCalQ;
    int   NBGX, NEDX, NSKPX, NBGY, NEDY, NSKPY, NBGZ, NEDZ, NSKPZ;
    int   IGREEN;
    int   nxt, nyt, nzt;
    MPI_Offset displacement, *displacement_sgt;
    float FL, FH, FP;
    char  INSRC[50], INVEL[50], OUT[50], CHKFILE[50];
    double GFLOPS = 1.0;
    double GFLOPS_SUM = 0.0;
    Grid3D u1=NULL, v1=NULL, w1=NULL;
    Grid3D d1=NULL, mu=NULL, lam=NULL;
    Grid3D xx=NULL, yy=NULL, zz=NULL, xy=NULL, yz=NULL, xz=NULL;
    Grid3D r1=NULL, r2=NULL, r3=NULL, r4=NULL, r5=NULL, r6=NULL;
    Grid3D qp=NULL, qs=NULL;
    PosInf tpsrc=NULL;
    Grid1D taxx=NULL, tayy=NULL, tazz=NULL, taxz=NULL, tayz=NULL, taxy=NULL; 
    Grid1D Bufx=NULL;
    Grid1D Bufy=NULL, Bufz=NULL;
    Grid3D vx1=NULL,   vx2=NULL,   lam_mu=NULL;
    Grid1D dcrjx=NULL, dcrjy=NULL, dcrjz=NULL;
    float vse[2], vpe[2], dde[2];
    FILE *fchk;
//  SGT CPU Define
    int SGT_BLOCK_SIZE, SGT_NUMBLOCKS;
    Grid3D sg1 =NULL, sg2 =NULL; 
    float *sgtBuf;
    char INSGT[50];
    int sgt_numsta_all;
    int sgt_numsta;
    int *sgt_sta;
    int *sgt_sta_indices;
    int sgtproc;
    int sgtmaster;  // master rank among SGT writers
    float *tmp_sgtBuf;
    int sgt_numstaPerWrite;
    int sgt_numstaLastWrite;
    int sgt_max_numsta;
    int sgt_nfiletypes;
    MPI_Datatype *sgt_filetypes;
    int sgt_max_nftypes;
//  GPU variables
    long int num_bytes;
    float* d_d1;
    float* d_u1;
    float* d_v1;
    float* d_w1;
    float* d_f_u1;
    float* d_f_v1;
    float* d_f_w1;
    float* d_b_u1;
    float* d_b_v1;
    float* d_b_w1;
    float* d_dcrjx;
    float* d_dcrjy;
    float* d_dcrjz;
    float* d_lam;
    float* d_mu;
    float* d_qp;
    float* d_qs;
    float* d_vx1;
    float* d_vx2;
    float* d_xx;
    float* d_yy;
    float* d_zz;
    float* d_xy;
    float* d_xz;
    float* d_yz;
    float* d_r1;
    float* d_r2;
    float* d_r3;
    float* d_r4;
    float* d_r5;
    float* d_r6;
    float* d_lam_mu;
    int*   d_tpsrc;
    float* d_taxx;
    float* d_tayy;
    float* d_tazz;
    float* d_taxz;
    float* d_tayz;
    float* d_taxy;
//  SGT GPU Define
    float* d_sg1;
    float* d_sg2;
    float* d_sgtBuf;
    int*   d_sgt_sta;
//  end of GPU variables
    int i,j,k,idx,idy,idz,ind;
    long int idtmp;
    long int tmpInd;
    const int maxdim = 3;
    float taumax, taumin, tauu;
    Grid3D tau=NULL, tau1=NULL, tau2=NULL;
    int npsrc;
    long int nt, cur_step, source_step;
    double time_un = 0.0;
    // time_src and time_mesh measures the time spent 
    // in source and mesh reading
    double time_src = 0.0, time_mesh = 0.0;
    double time_inisgt = 0.0;
    // time_fileio and time_gpuio measures the time spent
    // in file system IO and gpu memory copying for IO
    double time_fileio = 0.0, time_gpuio = 0.0;
    double time_gpuio_tmp = 0.0;
//  MPI+CUDA variables
    cudaError_t cerr;
    cudaStream_t stream_1, stream_2, stream_i;
    int   rank, size, err, srcproc, rank_gpu;
    int   dim[2], period[2], coord[2], reorder;
    //int   fmtype[3], fptype[3], foffset[3];
    int   x_rank_L  = -1,  x_rank_R  = -1,  y_rank_F = -1,  y_rank_B = -1;
    MPI_Comm MCW, MC1, MC_SGT;
    MPI_Request  request_x[4], request_y[4];
    MPI_Status   status_x[4],  status_y[4], filestatus;
    MPI_Datatype filetype;
    MPI_File fh;
    int   msg_v_size_x, msg_v_size_y, count_x = 0, count_y = 0;
    int   xls, xre, xvs, xve, xss1, xse1, xss2, xse2, xss3, xse3;
    int   yfs, yfe, ybs, ybe, yls,  yre;
    float* SL_vel;     // Velocity to be sent to   Left  in x direction (u1,v1,w1)
    float* SR_vel;     // Velocity to be Sent to   Right in x direction (u1,v1,w1)
    float* RL_vel;     // Velocity to be Recv from Left  in x direction (u1,v1,w1)
    float* RR_vel;     // Velocity to be Recv from Right in x direction (u1,v1,w1)
    float* SF_vel;     // Velocity to be sent to   Front in y direction (u1,v1,w1)
    float* SB_vel;     // Velocity to be Sent to   Back  in y direction (u1,v1,w1)
    float* RF_vel;     // Velocity to be Recv from Front in y direction (u1,v1,w1)
    float* RB_vel;     // Velocity to be Recv from Back  in y direction (u1,v1,w1)
//  variable definition ends    

    int tmpSize;
    int WRITE_STEP;
    int NTISKP;
    int NTISKP_SGT;
    int rec_NX;
    int rec_NY;
    int rec_NZ;
    int rec_nxt;
    int rec_nyt;
    int rec_nzt;
    int rec_nbgx;   // 0-based indexing, however NBG* is 1-based
    int rec_nedx;   // 0-based indexing, however NED* is 1-based
    int rec_nbgy;   // 0-based indexing
    int rec_nedy;   // 0-based indexing
    int rec_nbgz;   // 0-based indexing
    int rec_nedz;   // 0-based indexing
    char filename[50];
    char filenamebasex[50];
    char filenamebasey[50];
    char filenamebasez[50];
    char filenamebase_sgt[50];

    /*computation of moment and magnitude for kinematic source */
    float *mom, *d_mom, tmom=0.0f, gmom, mag;
    int n;

//  variable initialization begins 
    command(argc,argv,
        &TMAX,&DH,&DT,&ARBC,&PHT,&NPC,&ND,&NSRC,&NST,
        &NVAR,&NVE,&MEDIASTART,&IFAULT,&READ_STEP,&READ_STEP_GPU,
        &NTISKP,&WRITE_STEP,&NX,&NY,&NZ,&PX,&PY,
        &NBGX,&NEDX,&NSKPX,&NBGY,&NEDY,&NSKPY,&NBGZ,&NEDZ,&NSKPZ,
        &FL,&FH,&FP,&IDYNA,&SoCalQ,
        INSRC,INVEL,OUT,CHKFILE,&IGREEN,&NTISKP_SGT,INSGT);
    sprintf(filenamebasex,"%s/SX",OUT);
    sprintf(filenamebasey,"%s/SY",OUT);
    sprintf(filenamebasez,"%s/SZ",OUT);
    sprintf(filenamebase_sgt,"%s/SGT",OUT);

    //printf("After command.\n");
    // Below 12 lines are NOT for HPGPU4 machine!
/*    int local_rank;
    char* str;
    if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK")) != NULL) {
      local_rank = atoi(str);
      rank_gpu = local_rank%3;
    }
    else{
      printf("CANNOT READ LOCAL RANK!\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
*/
    //printf("%d) After rank_gpu calc.\n",local_rank);
//    rank_gpu = 0;
//    cudaSetDevice(rank_gpu);
    //if(local_rank==0) printf("after cudaSetDevice\n");
    hipInit(0); 
    MPI_Init(&argc,&argv);
     MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_dup(MPI_COMM_WORLD, &MCW );
    MPI_Barrier(MCW);
    nxt       = NX/PX;
    nyt       = NY/PY;
    nzt       = NZ;
    nt        = (int)(TMAX/DT) + 1;
    dim[0]    = PX;
    dim[1]    = PY;
    period[0] = 0;
    period[1] = 0;
    reorder   = 0;
    err       = MPI_Cart_create(MCW, 2, dim, period, reorder, &MC1);
    err       = MPI_Cart_shift(MC1, 0,  1,  &x_rank_L, &x_rank_R );
    err       = MPI_Cart_shift(MC1, 1,  1,  &y_rank_F, &y_rank_B ); 
    err       = MPI_Cart_coords(MC1, rank, 2, coord);
    err       = MPI_Barrier(MCW);
    
    // If any neighboring rank is out of bounds, then MPI_Cart_shift sets the
    // destination argument to a negative number. We use the convention that -1 
    // denotes ranks out of bounds.
    if (x_rank_L < 0) {
            x_rank_L = -1;
    }    
    if (x_rank_R < 0 ) {
            x_rank_R = -1;
    }    
    if (y_rank_F < 0) {
            y_rank_F = -1;
    }    
    if (y_rank_B < 0) {
            y_rank_B = -1;
    }    

    // Below 2 lines are only for HPGPU4 machine!
    rank_gpu = rank%8;
//    rank_gpu = 0;
    cudaSetDevice(rank_gpu);

printf("\n\nrank=%d) RS=%d, RSG=%d, NST=%d, IF=%d\n\n\n", 
rank, READ_STEP, READ_STEP_GPU, NST, IFAULT);

    // same for each processor:
    if(NEDX==-1) NEDX = NX;
    if(NEDY==-1) NEDY = NY;
    if(NEDZ==-1) NEDZ = NZ;
    // make NED's a record point
    // for instance if NBGX:NSKPX:NEDX = 1:3:9
    // then we have 1,4,7 but NEDX=7 is better
    NEDX = NEDX-(NEDX-NBGX)%NSKPX;
    NEDY = NEDY-(NEDY-NBGY)%NSKPY;
    NEDZ = NEDZ-(NEDZ-NBGZ)%NSKPZ;
    // number of recording points in total
    rec_NX = (NEDX-NBGX)/NSKPX+1;
    rec_NY = (NEDY-NBGY)/NSKPY+1;
    rec_NZ = (NEDZ-NBGZ)/NSKPZ+1;

    // specific to each processor:
    calcRecordingPoints(&rec_nbgx, &rec_nedx, &rec_nbgy, &rec_nedy, 
      &rec_nbgz, &rec_nedz, &rec_nxt, &rec_nyt, &rec_nzt, &displacement,
      (long int)nxt,(long int)nyt,(long int)nzt, rec_NX, rec_NY, rec_NZ, 
      NBGX,NEDX,NSKPX, NBGY,NEDY,NSKPY, NBGZ,NEDZ,NSKPZ, coord);
    printf("%d = (%d,%d)) NX,NY,NZ=%d,%d,%d\nnxt,nyt,nzt=%d,%d,%d\nrec_N=(%d,%d,%d)\nrec_nxt,=%d,%d,%d\nNBGX,SKP,END=(%d:%d:%d),(%d:%d:%d),(%d:%d:%d)\nrec_nbg,ed=(%d,%d),(%d,%d),(%d,%d)\ndisp=%ld\n",
        rank,coord[0],coord[1],NX,NY,NZ,nxt,nyt,nzt,
        rec_NX, rec_NY, rec_NZ, rec_nxt, rec_nyt, rec_nzt,
        NBGX,NSKPX,NEDX,NBGY,NSKPY,NEDY,NBGZ,NSKPZ,NEDZ,
        rec_nbgx,rec_nedx,rec_nbgy,rec_nedy,rec_nbgz,rec_nedz,(long int)displacement);

    //Check NZ against BLOCK_SIZE_Z
    if (NZ % BLOCK_SIZE_Z !=0) {
	fprintf(stderr, "Error in configuration: BLOCK_SIZE_Z (%d) must be a factor of NZ (%d), aborting.\n", BLOCK_SIZE_Z, NZ);
	exit(1);
    }
    
    int maxNX_NY_NZ_WS = (rec_NX>rec_NY?rec_NX:rec_NY);
    maxNX_NY_NZ_WS = (maxNX_NY_NZ_WS>rec_NZ?maxNX_NY_NZ_WS:rec_NZ);
    maxNX_NY_NZ_WS = (maxNX_NY_NZ_WS>WRITE_STEP?maxNX_NY_NZ_WS:WRITE_STEP);
    int ones[maxNX_NY_NZ_WS];
    MPI_Aint dispArray[maxNX_NY_NZ_WS];
    for(i=0;i<maxNX_NY_NZ_WS;++i){
      ones[i] = 1;
    }

    err = MPI_Type_contiguous(rec_nxt, MPI_FLOAT, &filetype);
    err = MPI_Type_commit(&filetype);
    for(i=0;i<rec_nyt;i++){
      dispArray[i] = sizeof(float);
      dispArray[i] = dispArray[i]*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(rec_nyt, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    for(i=0;i<rec_nzt;i++){
      dispArray[i] = sizeof(float);
      dispArray[i] = dispArray[i]*rec_NY*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(rec_nzt, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    for(i=0;i<WRITE_STEP;i++){
      dispArray[i] = sizeof(float);
      dispArray[i] = dispArray[i]*rec_NZ*rec_NY*rec_NX*i;
    }
    err = MPI_Type_create_hindexed(WRITE_STEP, ones, dispArray, filetype, &filetype);
    err = MPI_Type_commit(&filetype);
    MPI_Type_size(filetype, &tmpSize);
    if(rank==0) printf("filetype size (supposedly=rec_nxt*nyt*nzt*WS*4=%d) =%d\n", rec_nxt*rec_nyt*rec_nzt*WRITE_STEP*4,tmpSize);

/*
    fmtype[0]  = WRITE_STEP;
    //fmtype[1]  = NZ;
    fmtype[1]  = NY;
    fmtype[2]  = NX;
    fptype[0]  = WRITE_STEP;
    //fptype[1]  = nzt;
    fptype[1]  = nyt;
    fptype[2]  = nxt;
    foffset[0] = 0;
    //foffset[1] = 0;
    foffset[1] = nyt*coord[1];
    foffset[2] = nxt*coord[0];
    err = MPI_Type_create_subarray(3, fmtype, fptype, foffset, MPI_ORDER_C, MPI_FLOAT, &filetype);
    err = MPI_Type_commit(&filetype);
*/     
    if(rank==0) printf("rank=%d, x_rank_L=%d, x_rank_R=%d, y_rank_F=%d, y_rank_B=%d\n", rank, x_rank_L, x_rank_R, y_rank_F, y_rank_B);

    if(x_rank_L<0)
       xls = 2+4*loop;
    else
       xls = 4*loop;

    if(x_rank_R<0)
       xre = nxt+4*loop+1;
    else
       xre = nxt+4*loop+3;

    xvs   = 2+4*loop;
    xve   = nxt+4*loop+1;

    xss1  = xls;
    xse1  = 4*loop+3;
    xss2  = 4*loop+4;
    xse2  = nxt+4*loop-1;
    xss3  = nxt+4*loop;
    xse3  = xre;

    if(y_rank_F<0)
       yls = 2+4*loop;
    else
       yls = 4*loop;

    if(y_rank_B<0)
       yre = nyt+4*loop+1;
    else
       yre = nyt+4*loop+3;

    yfs  = 2+4*loop;
    yfe  = 2+8*loop-1;   
    ybs  = nyt+2;
    ybe  = nyt+4*loop+1;   

    if(IGREEN > -1){
      time_inisgt -= gethrtime();
      if(rank==0) printf("Before inisgt\n");
      err = inisgt(rank, INSGT, &sgt_numsta_all, &sgt_numsta, &sgt_max_numsta,
          &sgt_sta, &sgt_sta_indices, &sgtproc, MCW, NZ, nxt, nyt, nzt, coord);
      if(err){
        printf("ERROR: sgt initialization failed!\n"); return -1;
      }
      err = MPI_Comm_split(MCW, sgt_numsta>0, rank, &MC_SGT);
      if(rank==sgtproc){
        MPI_Allreduce(&rank, &sgtmaster, 1, MPI_INT, MPI_MIN, MC_SGT);
        int nIntBits = sizeof(int)*8-1;
        long int maxInt = pow(2,nIntBits);
        if((long long int)sizeof(float)*sgt_numsta*6*WRITE_STEP < maxInt){
          sgt_nfiletypes = 1;
          sgt_numstaPerWrite = sgt_numsta;
          sgt_numstaLastWrite = sgt_numsta;
        }
        else{
          sgt_nfiletypes = (int)((double)sgt_numsta*(double)sizeof(float)*(double)6*(double)WRITE_STEP/maxInt)+1;
          sgt_numstaPerWrite = (int)(maxInt/(double)sizeof(float)/(double)6/(double)WRITE_STEP);
          sgt_numstaLastWrite = (int)(sgt_numsta - (long long int)sgt_numstaPerWrite*(sgt_nfiletypes-1));
        }
        MPI_Allreduce(&sgt_nfiletypes, &sgt_max_nftypes, 1, MPI_INT, MPI_MAX, MC_SGT);
        sgt_filetypes = (MPI_Datatype*)malloc(sizeof(MPI_Datatype)*sgt_nfiletypes);
        num_bytes = sizeof(int)*3*sgt_numsta;
        cudaMalloc((void**)&d_sgt_sta, num_bytes);
        cudaMemset(d_sgt_sta, 0, num_bytes);
        cudaMemcpy(d_sgt_sta, sgt_sta, num_bytes, cudaMemcpyHostToDevice);
        num_bytes = sizeof(float)*6*sgt_numsta;
        cudaMalloc((void**)&d_sgtBuf,  num_bytes);
        cudaMemset(d_sgtBuf, 0, num_bytes);
        if(sgt_numsta < BLOCK_SIZE_Z){
          SGT_NUMBLOCKS = 1;
          SGT_BLOCK_SIZE = sgt_numsta;
        }
        else{
          SGT_BLOCK_SIZE = BLOCK_SIZE_Z;
          SGT_NUMBLOCKS = (sgt_numsta-1)/SGT_BLOCK_SIZE+1;
        }
        printf("sgtrank=%d NUMBLOCKS=%d, BLOCKSIZE=%d, n-filetypes=%d, numstaPerWrite=%d, numstaLastWrite=%d\n",
            rank,SGT_NUMBLOCKS,SGT_BLOCK_SIZE,
            sgt_nfiletypes, sgt_numstaPerWrite, sgt_numstaLastWrite);
        displacement_sgt = (MPI_Offset*)malloc(sizeof(MPI_Offset)*sgt_nfiletypes);
        // FAST-T
        for(i=0;i<sgt_nfiletypes;i++){
          err = MPI_Type_contiguous(WRITE_STEP*6, MPI_FLOAT, &sgt_filetypes[i]);
          MPI_Type_commit(&sgt_filetypes[i]);
          displacement_sgt[i] = 0;
        }
        int *ones_sgt = (int*)malloc(sizeof(int)*sgt_numstaPerWrite);
        MPI_Aint *dispArray_sgt = (MPI_Aint*)malloc(sizeof(MPI_Aint)*sgt_numsta);
        for(i=0;i<sgt_numsta;i++){
          ind = i/sgt_numstaPerWrite;
          // FAST-T
          dispArray_sgt[i] = (long long int)sizeof(float)*6*WRITE_STEP*sgt_sta_indices[i] - displacement_sgt[ind];
          if(i%sgt_numstaPerWrite==0){
            displacement_sgt[ind] = dispArray_sgt[i];
            dispArray_sgt[i] = 0;
          }
          if(i<sgt_numstaPerWrite)
            ones_sgt[i] = 1;
        }
        /*printf("sgt: rank=%d displacements are set: index:%d,%d,%d, dispArray:%d,%d,%d, displacement:%ld\n",
            rank,sgt_sta_indices[0],sgt_sta_indices[1],sgt_sta_indices[2],
            dispArray_sgt[0],dispArray_sgt[1],dispArray_sgt[2],displacement_sgt);*/
        tmpInd = 0;
        for(i=0;i<sgt_nfiletypes;i++){
          if(i==sgt_nfiletypes-1)
            err = MPI_Type_create_hindexed(sgt_numstaLastWrite, ones_sgt, &dispArray_sgt[tmpInd], sgt_filetypes[i], &sgt_filetypes[i]);
          else
            err = MPI_Type_create_hindexed(sgt_numstaPerWrite, ones_sgt, &dispArray_sgt[tmpInd], sgt_filetypes[i], &sgt_filetypes[i]);
          MPI_Type_commit(&sgt_filetypes[i]);
          tmpInd += sgt_numstaPerWrite;
        }
        MPI_Type_size(sgt_filetypes[0], &tmpSize);
        printf("SGT: rank=%d disp=%lld filetype_sgt size (supposedly=sgt_numstaPerWrite*6*WRITE_STEP*4=%d)=%d\n", 
            rank,displacement_sgt[0],sgt_numstaPerWrite*6*WRITE_STEP*4,tmpSize);
        if(sgt_nfiletypes>1){
          MPI_Type_size(sgt_filetypes[sgt_nfiletypes-1], &tmpSize);
          printf("SGT: rank=%d disp=%lld filetype_sgt_last size (supposedly=sgt_numstaLastWrite*6*WRITE_STEP*4=%d)=%d\n", 
            rank,displacement_sgt[0],sgt_numstaLastWrite*6*WRITE_STEP*4,tmpSize);
        }
        free(ones_sgt);
        free(dispArray_sgt);
      }
      time_inisgt += gethrtime();
      if(rank == sgtmaster)
        printf("SGT is initialized. Time elapsed (sec): %lf\n", time_inisgt);
    }

    if (IFAULT == 5){
       if (rank==0) fprintf(stdout, "Using IFAULT=5: kinematic source.\n");
       if ((NST != 2) || (READ_STEP != 2)) {
             if (rank==0) fprintf(stderr, "IFAULT=5 requires NST = READ_STEP =2.\nQuitting.");
          MPI_Finalize();
          return(-1);
          }
		mom = (float*)calloc(npsrc, sizeof(float));
		CUCHK(cudaMalloc((void**) &d_mom, num_bytes));
		CUCHK(cudaMemcpy(d_mom, mom, num_bytes, cudaMemcpyHostToDevice));
    }

    time_src -= gethrtime();
    if(rank==0) printf("Before inisource\n");
    err = inisource(rank,   IFAULT, NSRC,  READ_STEP, NST,   &srcproc, NZ, MCW, nxt, nyt, nzt, coord, maxdim, &npsrc,
                    &tpsrc, &taxx,  &tayy, &tazz,     &taxz, &tayz,    &taxy, INSRC);
    if(err)
    {
       printf("source initialization failed\n");
       return -1;
    }
    time_src += gethrtime();
    if(rank==0) printf("After inisource. Time elapsed (seconds): %lf\n", time_src);

    if(rank==srcproc)
    {
       printf("rank=%d, source rank, npsrc=%d\n", rank, npsrc);
       num_bytes = sizeof(float)*npsrc*READ_STEP_GPU;
       cudaMalloc((void**)&d_taxx, num_bytes);
       cudaMemset(d_taxx, 0, num_bytes);
       cudaMalloc((void**)&d_tayy, num_bytes);
       cudaMemset(d_tayy, 0, num_bytes);
       cudaMalloc((void**)&d_tazz, num_bytes);
       cudaMemset(d_tazz, 0, num_bytes);
       cudaMalloc((void**)&d_taxz, num_bytes);
       cudaMemset(d_taxz, 0, num_bytes);
       cudaMalloc((void**)&d_tayz, num_bytes);
       cudaMemset(d_tayz, 0, num_bytes);
       cudaMalloc((void**)&d_taxy, num_bytes);
       cudaMemset(d_taxy, 0, num_bytes);
       cudaMemcpy(d_taxx,taxx,num_bytes,cudaMemcpyHostToDevice);
       cudaMemcpy(d_tayy,tayy,num_bytes,cudaMemcpyHostToDevice);
       cudaMemcpy(d_tazz,tazz,num_bytes,cudaMemcpyHostToDevice);
       cudaMemcpy(d_taxz,taxz,num_bytes,cudaMemcpyHostToDevice);
       cudaMemcpy(d_tayz,tayz,num_bytes,cudaMemcpyHostToDevice);
       cudaMemcpy(d_taxy,taxy,num_bytes,cudaMemcpyHostToDevice);
       num_bytes = sizeof(int)*npsrc*maxdim;
       cudaMalloc((void**)&d_tpsrc, num_bytes);
       cudaMemset(d_tpsrc, 0, num_bytes);
       cudaMemcpy(d_tpsrc,tpsrc,num_bytes,cudaMemcpyHostToDevice);
    }

    d1     = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align); 
    mu     = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    lam    = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    lam_mu = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, 1); 

    if(NVE==1)
    { 
       qp   = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
       qs   = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    }

    time_mesh -= gethrtime();
    if(rank==0) printf("Before inimesh\n");
    inimesh(MEDIASTART, d1, mu, lam, qp, qs, &taumax, &taumin, NVAR, FP, FL, FH, 
            nxt, nyt, nzt, PX, PY, NX, NY, NZ, coord, MCW, IDYNA, NVE, SoCalQ, INVEL,
            vse, vpe, dde);
    time_mesh += gethrtime();
    if(rank==0) printf("After inimesh. Time elapsed (seconds): %lf\n", time_mesh);
    if(rank==0)
      writeCHK(CHKFILE, NTISKP, DT, DH, nxt, nyt, nzt,
        nt, ARBC, NPC, NVE, FL, FH, FP, vse, vpe, dde);

    mediaswap(d1, mu, lam, qp, qs, rank, x_rank_L, x_rank_R, y_rank_F, y_rank_B, nxt, nyt, nzt, MCW);

    for(i=xls;i<xre+1;i++)
      for(j=yls;j<yre+1;j++)
      {
         float t_xl, t_xl2m;
         t_xl             = 1.0/lam[i][j][nzt+awp_align-1];
         t_xl2m           = 2.0/mu[i][j][nzt+awp_align-1] + t_xl; 
         lam_mu[i][j][0]  = t_xl/t_xl2m;
      }

    num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop);
    cudaMalloc((void**)&d_lam_mu, num_bytes);
    cudaMemset(d_lam_mu, 0, num_bytes);
    cudaMemcpy(d_lam_mu,&lam_mu[0][0][0],num_bytes,cudaMemcpyHostToDevice);

    vx1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    vx2  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    if(NPC==0)
    {
	dcrjx = Alloc1D(nxt+4+8*loop);
        dcrjy = Alloc1D(nyt+4+8*loop);
        dcrjz = Alloc1D(nzt+2*awp_align);

        for(i=0;i<nxt+4+8*loop;i++)
	   dcrjx[i]  = 1.0;
        for(j=0;j<nyt+4+8*loop;j++)
           dcrjy[j]  = 1.0;
        for(k=0;k<nzt+2*awp_align;k++)
           dcrjz[k]  = 1.0;

        inicrj(ARBC, coord, nxt, nyt, nzt, NX, NY, ND, dcrjx, dcrjy, dcrjz);
    }

    if(NVE==1)
    {
        tau  = Alloc3D(2, 2, 2);
        tau1 = Alloc3D(2, 2, 2);
        tau2 = Alloc3D(2, 2, 2);
        tausub(tau, taumin, taumax);           
        float dt1 = 1.0/DT;
        for(i=0;i<2;i++)
          for(j=0;j<2;j++)
            for(k=0;k<2;k++)
            {
               tauu          = tau[i][j][k];
               tau1[i][j][k] = 1.0/((tauu*dt1)+(1.0/2.0));
               tau2[i][j][k] = (tauu*dt1)-(1.0/2.0);
            }

    	init_texture(nxt, nyt, nzt, tau1, tau2, vx1, vx2, xls, xre, yls, yre);

        Delloc3D(tau);
        Delloc3D(tau1);
        Delloc3D(tau2);
    }

    if(rank==0) printf("Allocate device media pointers and copy.\n");
    num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
    cudaMalloc((void**)&d_d1, num_bytes);
    cudaMemset(d_d1, 0, num_bytes);
    cudaMemcpy(d_d1,&d1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_lam, num_bytes);
    cudaMemset(d_lam, 0, num_bytes);
    cudaMemcpy(d_lam,&lam[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_mu, num_bytes);
    cudaMemset(d_mu, 0, num_bytes);
    cudaMemcpy(d_mu,&mu[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_qp, num_bytes);
    cudaMemset(d_qp, 0, num_bytes);
    cudaMemcpy(d_qp,&qp[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_qs, num_bytes);
    cudaMemset(d_qs, 0, num_bytes);
    cudaMemcpy(d_qs,&qs[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_vx1, num_bytes);
    cudaMemset(d_vx1, 0, num_bytes);
    cudaMemcpy(d_vx1,&vx1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_vx2, num_bytes);
    cudaMemset(d_vx2, 0, num_bytes);
    cudaMemcpy(d_vx2,&vx2[0][0][0],num_bytes,cudaMemcpyHostToDevice);

    /*BindArrayToTexture(d_vx1, d_vx2, num_bytes);*/

    if(NPC==0)
    {
    	num_bytes = sizeof(float)*(nxt+4+8*loop);
    	cudaMalloc((void**)&d_dcrjx, num_bytes);
    	cudaMemset(d_dcrjx, 0, num_bytes);
    	cudaMemcpy(d_dcrjx,dcrjx,num_bytes,cudaMemcpyHostToDevice);
        num_bytes = sizeof(float)*(nyt+4+8*loop);
        cudaMalloc((void**)&d_dcrjy, num_bytes);
        cudaMemset(d_dcrjy, 0, num_bytes);
        cudaMemcpy(d_dcrjy,dcrjy,num_bytes,cudaMemcpyHostToDevice);
        num_bytes = sizeof(float)*(nzt+2*awp_align);
        cudaMalloc((void**)&d_dcrjz, num_bytes);
        cudaMemset(d_dcrjz, 0, num_bytes);
        cudaMemcpy(d_dcrjz,dcrjz,num_bytes,cudaMemcpyHostToDevice);
    }

    if(rank==0) printf("Allocate host velocity and stress pointers.\n");
    u1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    v1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    w1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    xx  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    yy  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    zz  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    xy  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    yz  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    xz  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    if(NVE==1)
    {
        r1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        r2  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        r3  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        r4  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        r5  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        r6  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
    }

    if(IGREEN != -1)
    {
        sg1  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        sg2  = Alloc3D(nxt+4+8*loop, nyt+4+8*loop, nzt+2*awp_align);
        for(i=2+4*loop;i<nxt+2+4*loop;i++)
          for(j=2+4*loop;j<nyt+2+4*loop;j++)
            for(k=awp_align;k<nzt+awp_align;k++)
            {
              sg1[i][j][k]=(1./lam[i][j][k]+1./mu[i][j][k])/(1./mu[i][j][k]*(3.0/lam[i][j][k]+2.0/mu[i][j][k]));   
              sg2[i][j][k]=-1./lam[i][j][k]/(2.0/mu[i][j][k]*(3.0/lam[i][j][k]+2.0/mu[i][j][k]));
            }

        num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
        cudaMalloc((void**)&d_sg1, num_bytes);
        cudaMemset(d_sg1, 0, num_bytes);
        cudaMemcpy(d_sg1,&sg1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
        cudaMalloc((void**)&d_sg2, num_bytes);
        cudaMemset(d_sg2, 0, num_bytes);
        cudaMemcpy(d_sg2,&sg2[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    }

    source_step = 1;
    if(rank==srcproc)
    {
       printf("%d) add initial src\n", rank);
       addsrc(source_step, DH, DT, NST, npsrc, READ_STEP, maxdim, 
       IGREEN, nzt,  
       tpsrc, d1, u1, v1, w1, 
       taxx, tayy, tazz, taxz, tayz, taxy, xx, yy, zz, xy, yz, xz);
    }

    if(rank==0) printf("Allocate device velocity and stress pointers and copy.\n");
    num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
    cudaMalloc((void**)&d_u1, num_bytes);
    cudaMemset(d_u1, 0, num_bytes);
    cudaMemcpy(d_u1,&u1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_v1, num_bytes);
    cudaMemset(d_v1, 0, num_bytes);
    cudaMemcpy(d_v1,&v1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_w1, num_bytes);
    cudaMemset(d_w1, 0, num_bytes);
    cudaMemcpy(d_w1,&w1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_xx, num_bytes);
    cudaMemset(d_xx, 0, num_bytes);
    cudaMemcpy(d_xx,&xx[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_yy, num_bytes);
    cudaMemset(d_yy, 0, num_bytes);
    cudaMemcpy(d_yy,&yy[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_zz, num_bytes);
    cudaMemset(d_zz, 0, num_bytes);
    cudaMemcpy(d_zz,&zz[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_xy, num_bytes);
    cudaMemset(d_xy, 0, num_bytes);
    cudaMemcpy(d_xy,&xy[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_xz, num_bytes);
    cudaMemset(d_xz, 0, num_bytes);
    cudaMemcpy(d_xz,&xz[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    cudaMalloc((void**)&d_yz, num_bytes);
    cudaMemset(d_yz, 0, num_bytes);
    cudaMemcpy(d_yz,&yz[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    if(NVE==1)
    {
      if(rank==0) printf("Allocate additional device pointers (r) and copy.\n");
    	cudaMalloc((void**)&d_r1, num_bytes);
    	cudaMemset(d_r1, 0, num_bytes);
    	cudaMemcpy(d_r1,&r1[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&d_r2, num_bytes);
    	cudaMemset(d_r2, 0, num_bytes);
    	cudaMemcpy(d_r2,&r2[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&d_r3, num_bytes);
    	cudaMemset(d_r3, 0, num_bytes);
    	cudaMemcpy(d_r3,&r3[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&d_r4, num_bytes);
    	cudaMemset(d_r4, 0, num_bytes);
    	cudaMemcpy(d_r4,&r4[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&d_r5, num_bytes);
    	cudaMemset(d_r5, 0, num_bytes);
    	cudaMemcpy(d_r5,&r5[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    	cudaMalloc((void**)&d_r6, num_bytes);
    	cudaMemset(d_r6, 0, num_bytes);
    	cudaMemcpy(d_r6,&r6[0][0][0],num_bytes,cudaMemcpyHostToDevice);
    }
//  variable initialization ends
    if(rank==0) printf("Allocate buffers of #elements: %d\n",rec_nxt*rec_nyt*rec_nzt*WRITE_STEP);
    Bufx  = Alloc1D(rec_nxt*rec_nyt*rec_nzt*WRITE_STEP);
    Bufy  = Alloc1D(rec_nxt*rec_nyt*rec_nzt*WRITE_STEP);
    Bufz  = Alloc1D(rec_nxt*rec_nyt*rec_nzt*WRITE_STEP);
    if(IGREEN != -1){
      cudaMallocHost((void**)&tmp_sgtBuf,sgt_numsta*6*sizeof(float));
      cudaMallocHost((void**)&sgtBuf,sgt_numsta*6*sizeof(float)*WRITE_STEP);
    }
    num_bytes = sizeof(float)*3*(4*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
    cudaMallocHost((void**)&SL_vel, num_bytes);
    cudaMemset(SL_vel, 0, num_bytes);
    cudaMallocHost((void**)&SR_vel, num_bytes);
    cudaMemset(SR_vel, 0, num_bytes);
    cudaMallocHost((void**)&RL_vel, num_bytes);
    cudaMemset(RL_vel, 0, num_bytes);
    cudaMallocHost((void**)&RR_vel, num_bytes);
    cudaMemset(RR_vel, 0, num_bytes);
    num_bytes = sizeof(float)*3*(4*loop)*(nxt+4+8*loop)*(nzt+2*awp_align);
    cudaMallocHost((void**)&SF_vel, num_bytes);
    cudaMemset(SF_vel, 0, num_bytes);
    cudaMallocHost((void**)&SB_vel, num_bytes);
    cudaMemset(SB_vel, 0, num_bytes);
    cudaMallocHost((void**)&RF_vel, num_bytes);
    cudaMemset(RF_vel, 0, num_bytes);
    cudaMallocHost((void**)&RB_vel, num_bytes);
    cudaMemset(RB_vel, 0, num_bytes);
    num_bytes = sizeof(float)*(4*loop)*(nxt+4+8*loop)*(nzt+2*awp_align);
    cudaMalloc((void**)&d_f_u1, num_bytes);
    cudaMemset(d_f_u1, 0, num_bytes);
    cudaMalloc((void**)&d_f_v1, num_bytes);
    cudaMemset(d_f_v1, 0, num_bytes);
    cudaMalloc((void**)&d_f_w1, num_bytes);
    cudaMemset(d_f_w1, 0, num_bytes);
    cudaMalloc((void**)&d_b_u1, num_bytes);
    cudaMemset(d_b_u1, 0, num_bytes);
    cudaMalloc((void**)&d_b_v1, num_bytes);
    cudaMemset(d_b_v1, 0, num_bytes);
    cudaMalloc((void**)&d_b_w1, num_bytes);
    cudaMemset(d_b_w1, 0, num_bytes);
    msg_v_size_x = 3*(4*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
    msg_v_size_y = 3*(4*loop)*(nxt+4+8*loop)*(nzt+2*awp_align);
    SetDeviceConstValue(DH, DT, nxt, nyt, nzt);
    cudaStreamCreate(&stream_1);
    cudaStreamCreate(&stream_2);
    cudaStreamCreate(&stream_i);
    if(rank==0)
      fchk = fopen(CHKFILE,"a+");
    size_t cmemfree, cmemtotal;
    cudaMemGetInfo(&cmemfree, &cmemtotal);
    if(rank==0) printf("CUDA memory: master free=%ld\ttotal=%ld\n",cmemfree,cmemtotal);
    //if(rank==32) printf("CUDA memory: rank=32-largest sgt free=%ld\ttotal=%ld\n",cmemfree,cmemtotal);
    //if(rank==26) printf("CUDA memory: rank=26-largest source free=%ld\ttotal=%ld\n",cmemfree,cmemtotal);
//  Main Loop Starts
    if(NPC==0 && NVE==1)
    {
       time_un  -= gethrtime();
       //if(READ_STEP > nt) READ_STEP = nt;	
       //This loop has no loverlapping because there is source input
       for(cur_step=1;cur_step<=nt;cur_step++)
       {
         if(rank==0) printf("Time Step =                   %ld    OF  Total Timesteps = %ld\n", cur_step, nt); 
         if(cur_step==100 && rank==0) printf("Time per timestep:\t%f seconds\n",(gethrtime()+time_un)/100);
         cerr = cudaGetLastError();
         if(cerr!=cudaSuccess) printf("CUDA ERROR! rank=%d before timestep: %s\n",rank,cudaGetErrorString(cerr));
	 //pre-post MPI Message
         PostRecvMsg_Y(RF_vel, RB_vel, MCW, request_y, &count_y, msg_v_size_y, y_rank_F, y_rank_B);
 	 PostRecvMsg_X(RL_vel, RR_vel, MCW, request_x, &count_x, msg_v_size_x, x_rank_L, x_rank_R);
         //velocity computation in y boundary, two ghost cell regions
         dvelcy_H(d_u1, d_v1, d_w1, d_xx,   d_yy,   d_zz,   d_xy,       d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz,
                  d_d1, nxt,  nzt,  d_f_u1, d_f_v1, d_f_w1, stream_i,   yfs,  yfe, y_rank_F, rank);
         dvelcy_H(d_u1, d_v1, d_w1, d_xx,   d_yy,   d_zz,   d_xy,       d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz,
                  d_d1, nxt,  nzt,  d_b_u1, d_b_v1, d_b_w1, stream_i,   ybs,  ybe, y_rank_B, rank);
         Cpy2Host_VY(d_f_u1, d_f_v1, d_f_w1,  SF_vel, nxt, nzt, stream_i, y_rank_F);
         Cpy2Host_VY(d_b_u1, d_b_v1, d_b_w1,  SB_vel, nxt, nzt, stream_i, y_rank_B);
         cudaThreadSynchronize();
         //velocity communication in y direction
         PostSendMsg_Y(SF_vel, SB_vel, MCW, request_y, &count_y, msg_v_size_y, y_rank_F, y_rank_B, rank, Both);
         MPI_Waitall(count_y, request_y, status_y);
         Cpy2Device_VY(d_u1,     d_v1,     d_w1,     d_f_u1, d_f_v1, d_f_w1, d_b_u1, d_b_v1, d_b_w1, RF_vel, RB_vel, nxt, nyt, nzt, 
                       stream_i, stream_i, y_rank_F, y_rank_B, rank);
         //velocity computation whole 3D Grid (nxt, nyt, nzt)
         dvelcx_H(d_u1, d_v1, d_w1, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz, 
                  d_d1, nyt,  nzt,  stream_i,   xvs,  xve, rank);
         Cpy2Host_VX(d_u1, d_v1, d_w1, SL_vel, nxt, nyt, nzt, stream_i, x_rank_L, Left);
         Cpy2Host_VX(d_u1, d_v1, d_w1, SR_vel, nxt, nyt, nzt, stream_i, x_rank_R, Right);
	 cudaThreadSynchronize();
         //velocity communication in x direction
	 PostSendMsg_X(SL_vel, SR_vel, MCW, request_x, &count_x, msg_v_size_x, x_rank_L, x_rank_R, rank, Both);
	 MPI_Waitall(count_x, request_x, status_x);
         Cpy2Device_VX(d_u1, d_v1, d_w1, RL_vel, RR_vel, nxt, nyt, nzt, stream_i, stream_i, x_rank_L, x_rank_R);
	 //stress computation whole 3D Grid (nxt+4, nyt+4, nzt)
         dstrqc_H(d_xx, d_yy, d_zz, d_xy,    d_xz,    d_yz,    d_r1, d_r2, d_r3,     d_r4,     d_r5, d_r6,     d_u1, d_v1, d_w1, d_lam,
                  d_mu, d_qp, d_qs, d_dcrjx, d_dcrjy, d_dcrjz, nyt,  nzt,  stream_i, d_lam_mu, NX,   coord[0], coord[1],   xls,  xre,
                  yls,  yre, rank, d_vx1, d_vx2);
         //update source input
        //**Below was added in 2021 as part of BBP verification work**
		if (IFAULT==5) {
                CUCHK(cudaDeviceSynchronize());
                addkinsrc_H(cur_step, maxdim, d_tpsrc, npsrc, stream_i, d_mu,
            d_taxx, d_tayy, d_tazz, d_taxz, d_tayz, d_taxy,
            d_xx, d_yy, d_zz, d_xy, d_yz, d_xz, d_mom);
		}

         if(rank==srcproc && cur_step<NST && IFAULT<4)
         {
            ++source_step;
/*            //printf("call addsrc_H rank=%d cur_step=%ld\n",rank,cur_step);
int i_,j_,k_,pos;
    i_=2+4*loop;
    j_=2+4*loop;
    k_=awp_align;
    pos = i_*(nyt+4+8*loop)*(nzt+2*awp_align)+j_*(nzt+2*awp_align)+k_;
printf("xx,yy,zz,xy,xz,yz=%f,%f,%f,%f,%f,%f\naxx,yy,zz,xy,xz,yz=%f,%f,%f,%f,%f,%f\n",
      xx[i_][j_][k_],yy[i_][j_][k_],zz[i_][j_][k_],xy[i_][j_][k_],xz[i_][j_][k_],yz[i_][j_][k_],
      taxx[cur_step],tayy[cur_step],tazz[cur_step],taxy[cur_step],taxz[cur_step],tayz[cur_step]);
*/
	            addsrc_H(source_step, READ_STEP_GPU, maxdim, d_tpsrc, npsrc, stream_i, 
                    IGREEN, nzt, d_d1, d_u1, d_v1, d_w1, 
                    d_taxx, d_tayy, d_tazz, d_taxz, d_tayz, d_taxy,
                    d_xx,       d_yy,      d_zz,   d_xy,    d_yz,  d_xz);
         }
         if((IGREEN != -1) && (cur_step%NTISKP_SGT==0) && (sgtproc == rank))
         { 
            ComputeSGT(d_xx,          d_yy,      d_zz,   d_xy,      d_xz,         d_yz,
                       d_sg1,         d_sg2,     d_mu,   sgt_numsta, d_sgt_sta, d_sgtBuf,
                       stream_i,      nzt,       SGT_BLOCK_SIZE,    SGT_NUMBLOCKS, d_qs, d_d1);
            cudaThreadSynchronize();
            cerr = cudaGetLastError();
            if(cerr!=cudaSuccess) printf("CUDA ERROR! rank=%d after threadSync: %s\n",rank,cudaGetErrorString(cerr));
            num_bytes = sizeof(float)*sgt_numsta*6;
            idtmp = (cur_step/NTISKP_SGT-1)%WRITE_STEP;
            // FAST-T - 2 lines
            tmpInd = idtmp;
            if(rank == sgtmaster) time_gpuio_tmp = -gethrtime();
            cerr = cudaMemcpy(tmp_sgtBuf,d_sgtBuf,num_bytes,cudaMemcpyDeviceToHost);
            if(cerr!=cudaSuccess) printf("CUDA ERROR! rank=%d copying d_sgtBuf to sgtBuf: %s\n",rank,cudaGetErrorString(cerr));
            // FAST-T
            for(i=0;i<sgt_numsta;i++){
              for(j=0;j<6;j++){
                sgtBuf[tmpInd] = tmp_sgtBuf[i*6+j];
                tmpInd += WRITE_STEP;
              }
            }
            if(rank == sgtmaster){
              time_gpuio_tmp += gethrtime();
              time_gpuio += time_gpuio_tmp;
              printf("Output data buffered in (sec): %lf\n",time_gpuio_tmp);
            }
            // FAST-T end
            if((cur_step/NTISKP_SGT)%WRITE_STEP == 0){
              sprintf(filename, "%s%07ld", filenamebase_sgt, cur_step);
              time_fileio = -gethrtime();
              tmpInd = 0;
              err = MPI_File_open(MC_SGT,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
              for(i=0;i<sgt_max_nftypes;i++){
                printf("SGT rank=%d (%d/%d) tmpInd=%ld\n",rank,i,sgt_nfiletypes,tmpInd);
                if(i<sgt_nfiletypes){
                  err = MPI_File_set_view(fh, displacement_sgt[i], MPI_FLOAT, sgt_filetypes[i], "native", MPI_INFO_NULL);
                  printf("SGT rank=%d (%d/%d) File_set_view: disp=%lld\n",rank,i,sgt_nfiletypes,displacement_sgt[i]);
                  if(i == sgt_nfiletypes-1)
                    err = MPI_File_write(fh, &sgtBuf[tmpInd], sgt_numstaLastWrite*6*WRITE_STEP, MPI_FLOAT, &filestatus);
                  else
                    err = MPI_File_write(fh, &sgtBuf[tmpInd], sgt_numstaPerWrite*6*WRITE_STEP, MPI_FLOAT, &filestatus);
                  printf("SGT rank=%d (%d/%d) File_write\n",rank,i,sgt_nfiletypes);
                  tmpInd += sgt_numstaPerWrite*6*WRITE_STEP;
                }
                else{
                  err = MPI_File_set_view(fh, rank*sizeof(float), MPI_FLOAT, MPI_FLOAT, "native", MPI_INFO_NULL);
                  //printf("SGT rank=%d (%d/%d) File_set_view: disp=%lld\n",rank,i,sgt_nfiletypes,displacement_sgt[i]);
                  //err = MPI_File_write(fh, &sgtBuf[0], 0, MPI_FLOAT, &filestatus);
                  printf("SGT rank=%d (%d/%d) No file_write\n",rank,i,sgt_nfiletypes);
                }
              }
              err = MPI_File_close(&fh);
              time_fileio += gethrtime();
              printf("SGT rank=%d Output data written in (sec): %lf\n",rank,time_fileio);
            }
         }
         cudaThreadSynchronize();
 
         if(cur_step%NTISKP == 0){
          num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
          cudaMemcpy(&u1[0][0][0],d_u1,num_bytes,cudaMemcpyDeviceToHost);
          cudaMemcpy(&v1[0][0][0],d_v1,num_bytes,cudaMemcpyDeviceToHost);
          cudaMemcpy(&w1[0][0][0],d_w1,num_bytes,cudaMemcpyDeviceToHost);
          idtmp = ((cur_step/NTISKP+WRITE_STEP-1)%WRITE_STEP);
          idtmp = idtmp*rec_nxt*rec_nyt*rec_nzt;
          tmpInd = idtmp;
          //if(rank==0) printf("idtmp=%ld\n", idtmp);
          // surface: k=nzt+awp_align-1;
          for(k=nzt+awp_align-1 - rec_nbgz; k>=nzt+awp_align-1 - rec_nedz; k=k-NSKPZ)
            for(j=2+4*loop + rec_nbgy; j<=2+4*loop + rec_nedy; j=j+NSKPY)
              for(i=2+4*loop + rec_nbgx; i<=2+4*loop + rec_nedx; i=i+NSKPX)
              {
                //idx = (i-2-4*loop)/NSKPX;
                //idy = (j-2-4*loop)/NSKPY;
                //idz = ((nzt+awp_align-1) - k)/NSKPZ;
                //tmpInd = idtmp + idz*rec_nxt*rec_nyt + idy*rec_nxt + idx;
                //if(rank==0) printf("%ld:%d,%d,%d\t",tmpInd,i,j,k);
                printf("%d) Writing value for local point (%d, %d, %d), which corresponds to array index (%d, %d, %d)\n", rank, i-(2+4*loop), j-(2+4*loop), (-1*k)+(nzt+awp_align-1), i, j, k);
                Bufx[tmpInd] = u1[i][j][k];
                Bufy[tmpInd] = v1[i][j][k];
                Bufz[tmpInd] = w1[i][j][k];
                tmpInd++;
              }
          if((cur_step/NTISKP)%WRITE_STEP == 0){
            cudaThreadSynchronize();
            sprintf(filename, "%s%07ld", filenamebasex, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
	    if (err!=MPI_SUCCESS) {
			printf("Error on open: %d\n", err);
			exit(1);
		}
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            if (err!=MPI_SUCCESS) {
                        printf("Error on set view: %d\n", err);
                        exit(1);
                }
            err = MPI_File_write_all(fh, Bufx, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
	    if (err!=MPI_SUCCESS) {
                        printf("Error on write_all: %d\n", err);
			printf("Buffer size: %ld\n", rec_nxt*rec_nyt*rec_nzt*WRITE_STEP*sizeof(float));
                        exit(1);
                }
            err = MPI_File_close(&fh);
            if (err!=MPI_SUCCESS) {
                        printf("Error on close: %d\n", err);
                        exit(1);
                }
            sprintf(filename, "%s%07ld", filenamebasey, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(fh, Bufy, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
            err = MPI_File_close(&fh);
            sprintf(filename, "%s%07ld", filenamebasez, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(fh, Bufz, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
            err = MPI_File_close(&fh);
          }
          //else 
            //cudaThreadSynchronize();
          // write-statistics to chk file:
          if(rank==0){
            i = ND+2+4*loop;
            j = i;
            k = nzt+awp_align-1-ND;
            fprintf(fchk,"%ld :\t%e\t%e\t%e\n",cur_step,u1[i][j][k],v1[i][j][k],w1[i][j][k]);
            fflush(fchk);
          }
         }
         //else
          //cudaThreadSynchronize();

          /*
          if((cur_step<NST) && (cur_step%25==0) && (rank==srcproc)){
            printf("%d) SOURCE: taxx,xy,xz:%e,%e,%e\n",rank,
                taxx[cur_step],taxy[cur_step],taxz[cur_step]);
          }*/
       }
/*
       //This loop is purly stencil computation, overlapping are utilized here
       //time_un  -= gethrtime();
       for(cur_step=cur_step;cur_step<=nt;cur_step++)
       {
         if(rank==0) printf("Time Step =                   %ld    OF  Total Timesteps = %ld\n", cur_step, nt);
         //pre-post MPI Message
         PostRecvMsg_Y(RF_vel, RB_vel, MCW, request_y, &count_y, msg_v_size_y, y_rank_F, y_rank_B);
         PostRecvMsg_X(RL_vel, RR_vel, MCW, request_x, &count_x, msg_v_size_x, x_rank_L, x_rank_R);
         //velocity computation in y boundary, two ghost cell region, two different streams to control
         dvelcy_H(d_u1, d_v1, d_w1, d_xx,   d_yy,   d_zz,   d_xy,       d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz,
                  d_d1, nxt,  nzt,  d_f_u1, d_f_v1, d_f_w1, stream_1,   yfs,  yfe,  y_rank_F, rank);
         Cpy2Host_VY(d_f_u1, d_f_v1, d_f_w1,  SF_vel, nxt, nzt, stream_1, y_rank_F);
         dvelcy_H(d_u1, d_v1, d_w1, d_xx,   d_yy,   d_zz,   d_xy,       d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz,
                  d_d1, nxt,  nzt,  d_b_u1, d_b_v1, d_b_w1, stream_2,   ybs,  ybe,  y_rank_B, rank);
         Cpy2Host_VY(d_b_u1, d_b_v1, d_b_w1,  SB_vel, nxt, nzt, stream_2, y_rank_B);
         //Memory copy from GPU to CPU, and velocity computation at the same time
         dvelcx_H(d_u1, d_v1, d_w1, d_xx, d_yy, d_zz, d_xy, d_xz, d_yz, d_dcrjx, d_dcrjy, d_dcrjz,
                  d_d1, nyt,  nzt,  stream_i,   xvs,  xve, rank);
         //MPI overlapping velocity computation
         cudaStreamSynchronize(stream_1);
         PostSendMsg_Y(SF_vel, SB_vel, MCW, request_y, &count_y, msg_v_size_y, y_rank_F, y_rank_B, rank, Front);
         cudaStreamSynchronize(stream_2);
         PostSendMsg_Y(SF_vel, SB_vel, MCW, request_y, &count_y, msg_v_size_y, y_rank_F, y_rank_B, rank, Back);
         MPI_Waitall(count_y, request_y, status_y);
         Cpy2Device_VY(d_u1,     d_v1,     d_w1,     d_f_u1, d_f_v1, d_f_w1, d_b_u1, d_b_v1, d_b_w1, RF_vel, RB_vel, nxt, nyt, nzt,
                       stream_1, stream_2, y_rank_F, y_rank_B, rank);
         cudaThreadSynchronize();
         //start stress computation in insider part
         dstrqc_H(d_xx, d_yy, d_zz, d_xy,    d_xz,    d_yz,    d_r1, d_r2, d_r3,     d_r4,     d_r5, d_r6,     d_u1, d_v1, d_w1, d_lam,
                  d_mu, d_qp, d_qs, d_dcrjx, d_dcrjy, d_dcrjz, nyt,  nzt,  stream_i, d_lam_mu, NX,   coord[0], coord[1],   xss2, xse2,
                  yls,  yre, rank);
         Cpy2Host_VX(d_u1, d_v1, d_w1, SL_vel, nxt, nyt, nzt, stream_1, x_rank_L, Left);
         Cpy2Host_VX(d_u1, d_v1, d_w1, SR_vel, nxt, nyt, nzt, stream_2, x_rank_R, Right);
         //MPI overlapping stress computation
         cudaStreamSynchronize(stream_1);
         PostSendMsg_X(SL_vel, SR_vel, MCW, request_x, &count_x, msg_v_size_x, x_rank_L, x_rank_R, rank, Left);
         cudaStreamSynchronize(stream_2);
         PostSendMsg_X(SL_vel, SR_vel, MCW, request_x, &count_x, msg_v_size_x, x_rank_L, x_rank_R, rank, Right);
         MPI_Waitall(count_x, request_x, status_x);
         Cpy2Device_VX(d_u1, d_v1, d_w1, RL_vel, RR_vel, nxt, nyt, nzt, stream_1, stream_2, x_rank_L, x_rank_R);
	 //stress computation in ghost cells
         dstrqc_H(d_xx, d_yy, d_zz, d_xy,    d_xz,    d_yz,    d_r1, d_r2, d_r3,     d_r4,     d_r5, d_r6,     d_u1, d_v1, d_w1, d_lam,
                  d_mu, d_qp, d_qs, d_dcrjx, d_dcrjy, d_dcrjz, nyt,  nzt,  stream_1, d_lam_mu, NX,   coord[0], coord[1],   xss1, xse1,
		  yls,  yre, rank);
         dstrqc_H(d_xx, d_yy, d_zz, d_xy,    d_xz,    d_yz,    d_r1, d_r2, d_r3,     d_r4,     d_r5, d_r6,     d_u1, d_v1, d_w1, d_lam,
                  d_mu, d_qp, d_qs, d_dcrjx, d_dcrjy, d_dcrjz, nyt,  nzt,  stream_2, d_lam_mu, NX,   coord[0], coord[1],   xss3, xse3,
                  yls,  yre, rank);
         if((IGREEN != -1) && (cur_step%NTISKP_SGT==0) && (sgtproc == rank))
         { 
            ComputeSGT(d_xx,          d_yy,      d_zz,   d_xy,      d_xz,         d_yz,
//                       d_sgt1,        d_sgt2,    d_sgt3, d_sgt4,    d_sgt5,       d_sgt6,
                       d_sg1,         d_sg2,     d_mu,   sgt_numsta, d_sgt_sta, d_sgtBuf,
                       stream_i,      nzt,       SGT_BLOCK_SIZE,    SGT_NUMBLOCKS);
            cudaThreadSynchronize();
            //num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
            num_bytes = sizeof(float)*sgt_numsta*6;
            idtmp = (cur_step/NTISKP_SGT-1)%WRITE_STEP;
            idtmp *= sgt_numsta*6;
            cudaMemcpy(sgtBuf+idtmp,d_sgtBuf,num_bytes,cudaMemcpyDeviceToHost);
            //memcpy(sgtBuf[idtmp],(void*)tmp_sgtBuf,num_bytes);

            if((cur_step/NTISKP_SGT)%WRITE_STEP == 0){
              sprintf(filename, "%s%07ld", filenamebase_sgt, cur_step);
              err = MPI_File_open(MC_SGT,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
              err = MPI_File_set_view(fh, displacement_sgt, MPI_FLOAT, filetype_sgt, "native", MPI_INFO_NULL);
              err = MPI_File_write_all(fh, sgtBuf, sgt_numsta*6*WRITE_STEP, MPI_FLOAT, &filestatus);
              err = MPI_File_close(&fh);
            }
         }
         else cudaThreadSynchronize();
         
         if(cur_step%NTISKP == 0){
          num_bytes = sizeof(float)*(nxt+4+8*loop)*(nyt+4+8*loop)*(nzt+2*awp_align);
          cudaMemcpy(&u1[0][0][0],d_u1,num_bytes,cudaMemcpyDeviceToHost);
          cudaMemcpy(&v1[0][0][0],d_v1,num_bytes,cudaMemcpyDeviceToHost);
          cudaMemcpy(&w1[0][0][0],d_w1,num_bytes,cudaMemcpyDeviceToHost);
          idtmp = ((cur_step/NTISKP+WRITE_STEP-1)%WRITE_STEP);
          idtmp = idtmp*rec_nxt*rec_nyt*rec_nzt;
          tmpInd = idtmp;
          // surface: k=nzt+awp_align-1;
          for(k=nzt+awp_align-1 - rec_nbgz; k>=nzt+awp_align-1 - rec_nedz; k=k-NSKPZ)
            for(j=2+4*loop + rec_nbgy; j<=2+4*loop + rec_nedy; j=j+NSKPY)
              for(i=2+4*loop + rec_nbgx; i<=2+4*loop + rec_nedx; i=i+NSKPX)
              {
                //idx = (i-2-4*loop)/NSKPX;
                //idy = (j-2-4*loop)/NSKPY;
                //idz = ((nzt+awp_align-1) - k)/NSKPZ;
                //tmpInd = idtmp + idz*rec_nxt*rec_nyt + idy*rec_nxt + idx;
                Bufx[tmpInd] = u1[i][j][k];
                Bufy[tmpInd] = v1[i][j][k];
                Bufz[tmpInd] = w1[i][j][k];
                tmpInd++;
              }
          if((cur_step/NTISKP)%WRITE_STEP == 0){
            cudaThreadSynchronize();
            //printf("I'm %d, my disp=%ld\n", rank, displacement);
            sprintf(filename, "%s%07ld", filenamebasex, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(fh, Bufx, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
            err = MPI_File_close(&fh);
            sprintf(filename, "%s%07ld", filenamebasey, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(fh, Bufy, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
            err = MPI_File_close(&fh);
            sprintf(filename, "%s%07ld", filenamebasez, cur_step);
            err = MPI_File_open(MCW,filename,MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL,&fh);
            err = MPI_File_set_view(fh, displacement, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
            err = MPI_File_write_all(fh, Bufz, rec_nxt*rec_nyt*rec_nzt*WRITE_STEP, MPI_FLOAT, &filestatus);
            err = MPI_File_close(&fh);
          }
          else
            cudaThreadSynchronize();
          if(rank==0){
            i = ND+2+4*loop;
            j = i;
            k = nzt+awp_align-1-ND;
            fprintf(fchk,"%ld :\t%e\t%e\t%e\n",cur_step,u1[i][j][k],v1[i][j][k],w1[i][j][k]);
            fflush(fchk);
          }
         }
         else
           cudaThreadSynchronize();

       }*/
       time_un += gethrtime();
    } 
    if(rank==0){
      fprintf(fchk,"END\n");
      fclose(fchk);
    }

    cudaStreamDestroy(stream_1);
    cudaStreamDestroy(stream_2);
    cudaStreamDestroy(stream_i);
    cudaFreeHost(SL_vel);
    cudaFreeHost(SR_vel);
    cudaFreeHost(RL_vel);
    cudaFreeHost(RR_vel);
    cudaFreeHost(SF_vel);
    cudaFreeHost(SB_vel);
    cudaFreeHost(RF_vel);
    cudaFreeHost(RB_vel);
    GFLOPS  = 1.0;
    GFLOPS  = GFLOPS*307.0*(xre - xls)*(yre-yls)*nzt;
    GFLOPS  = GFLOPS/(1000*1000*1000);
    //time_un = time_un/(cur_step-READ_STEP);
    time_un = time_un/cur_step;
    GFLOPS  = GFLOPS/time_un;
    MPI_Allreduce( &GFLOPS, &GFLOPS_SUM, 1, MPI_DOUBLE, MPI_SUM, MCW );
    double time_src_max;
    MPI_Allreduce( &time_src, &time_src_max, 1, MPI_DOUBLE, MPI_MAX, MCW);
    double time_fileio_max;
    MPI_Allreduce( &time_fileio, &time_fileio_max, 1, MPI_DOUBLE, MPI_MAX, MCW);
    if(rank==0)
    {
        printf("GPU benchmark size NX=%d, NY=%d, NZ=%d, ReadStep=%d\n", NX, NY, NZ, READ_STEP);
    	  printf("GPU computing flops=%1.18f GFLOPS, time = %1.18f secs per timestep\n", GFLOPS_SUM, time_un);
        printf("GPU total I/O buffering time (memcpy+buffer)=%lf secs\n", time_gpuio);
        printf("GPU maximum I/O time (file system)=%lf secs\n", time_fileio_max);
        printf("GPU source reading time=%lf secs\nGPU mesh reading time=%lf secs\n", time_src_max, time_mesh);
    }	
//  Main Loop Ends
 
//  program ends, free all memories
    /*UnBindArrayFromTexture();*/
    Delloc3D(u1);
    Delloc3D(v1);
    Delloc3D(w1); 
    Delloc3D(xx);
    Delloc3D(yy);
    Delloc3D(zz);
    Delloc3D(xy);
    Delloc3D(yz);
    Delloc3D(xz);
    Delloc3D(vx1);
    Delloc3D(vx2);

    cudaFree(d_u1);
    cudaFree(d_v1);
    cudaFree(d_w1);
    cudaFree(d_f_u1);
    cudaFree(d_f_v1);
    cudaFree(d_f_w1);
    cudaFree(d_b_u1);
    cudaFree(d_b_v1);
    cudaFree(d_b_w1);
    cudaFree(d_xx);
    cudaFree(d_yy);
    cudaFree(d_zz);
    cudaFree(d_xy);
    cudaFree(d_yz);
    cudaFree(d_xz);
    cudaFree(d_vx1);
    cudaFree(d_vx2);

    if(NVE==1)
    {
       Delloc3D(r1);
       Delloc3D(r2);
       Delloc3D(r3);
       Delloc3D(r4);
       Delloc3D(r5);
       Delloc3D(r6);
       cudaFree(d_r1);
       cudaFree(d_r2);
       cudaFree(d_r3);
       cudaFree(d_r4);
       cudaFree(d_r5);
       cudaFree(d_r6);

       Delloc3D(qp);
       Delloc3D(qs);
       cudaFree(d_qp);
       cudaFree(d_qs);
    }

    if(IGREEN != -1)
    {
       Delloc3D(sg1);
       Delloc3D(sg2);
/*
       Delloc3D(sgt1);
       Delloc3D(sgt2);
       Delloc3D(sgt3);
       Delloc3D(sgt4);
       Delloc3D(sgt5);
       Delloc3D(sgt6);
*/
       cudaFree(d_sg1);
       cudaFree(d_sg2);
/*
       cudaFree(d_sgt1);
       cudaFree(d_sgt2);
       cudaFree(d_sgt3);
       cudaFree(d_sgt4);
       cudaFree(d_sgt5);
       cudaFree(d_sgt6);
*/
       cudaFreeHost(tmp_sgtBuf);
       cudaFreeHost(sgtBuf);
    }

    if(NPC==0)
    {
        Delloc1D(dcrjx);
        Delloc1D(dcrjy);
        Delloc1D(dcrjz);
        cudaFree(d_dcrjx);
        cudaFree(d_dcrjy);
        cudaFree(d_dcrjz);
    }

    Delloc3D(d1);
    Delloc3D(mu);
    Delloc3D(lam);
    Delloc3D(lam_mu);
    cudaFree(d_d1);
    cudaFree(d_mu);
    cudaFree(d_lam);
    cudaFree(d_lam_mu);

    if(rank==srcproc)
    {
       Delloc1D(taxx);
       Delloc1D(tayy);
       Delloc1D(tazz);
       Delloc1D(taxz);
       Delloc1D(tayz);
       Delloc1D(taxy); 
       cudaFree(d_taxx);
       cudaFree(d_tayy);
       cudaFree(d_tazz);
       cudaFree(d_taxz);
       cudaFree(d_tayz);
       cudaFree(d_taxy);
       Delloc1P(tpsrc);
       cudaFree(d_tpsrc);
    }

    MPI_Comm_free( &MC1 );
    if (IGREEN!=-1) {
	    MPI_Comm_free( &MC_SGT );
    }
    MPI_Finalize();
    return (0);
}

// Calculates recording points for each core
// rec_nbgxyz rec_nedxyz...
// WARNING: Assumes NPZ = 1! Only surface outputs are needed!
void calcRecordingPoints(int *rec_nbgx, int *rec_nedx, 
  int *rec_nbgy, int *rec_nedy, int *rec_nbgz, int *rec_nedz, 
  int *rec_nxt, int *rec_nyt, int *rec_nzt, MPI_Offset *displacement,
  long int nxt, long int nyt, long int nzt, int rec_NX, int rec_NY, int rec_NZ,
  int NBGX, int NEDX, int NSKPX, int NBGY, int NEDY, int NSKPY, 
  int NBGZ, int NEDZ, int NSKPZ, int *coord){

  *displacement = 0;

  if(NBGX > nxt*(coord[0]+1))     *rec_nxt = 0;
  else if(NEDX < nxt*coord[0]+1)  *rec_nxt = 0;
  else{
    if(nxt*coord[0] >= NBGX){
      *rec_nbgx = (nxt*coord[0]+NBGX-1)%NSKPX;
      //*rec_nbgx = NSKPX-(nxt*coord[0]-NBGX+1)%NSKPX;
      *displacement += (nxt*coord[0]-NBGX)/NSKPX+1;
    }
    else
      *rec_nbgx = NBGX-nxt*coord[0]-1;  // since rec_nbgx is 0-based
    if(nxt*(coord[0]+1) <= NEDX)
      *rec_nedx = (nxt*(coord[0]+1)+NBGX-1)%NSKPX-NSKPX+nxt;
      //*rec_nedx = nxt-(nxt*(coord[0]+1)-1-NBGX+1)%NSKPX;
    else
      *rec_nedx = NEDX-nxt*coord[0]-1;
    *rec_nxt = (*rec_nedx-*rec_nbgx)/NSKPX+1;
  }

  if(NBGY > nyt*(coord[1]+1))     *rec_nyt = 0;
  else if(NEDY < nyt*coord[1]+1)  *rec_nyt = 0;
  else{
    if(nyt*coord[1] >= NBGY){
      *rec_nbgy = (nyt*coord[1]+NBGY-1)%NSKPY;
      //*rec_nbgy = NSKPY-(nyt*coord[1]+NBGY-1)%NSKPY;
      *displacement += ((nyt*coord[1]-NBGY)/NSKPY+1)*rec_NX;
    }
    else
      *rec_nbgy = NBGY-nyt*coord[1]-1;  // since rec_nbgy is 0-based
    if(nyt*(coord[1]+1) <= NEDY)
      *rec_nedy = (nyt*(coord[1]+1)+NBGY-1)%NSKPY-NSKPY+nyt;
      //*rec_nedy = nyt-(nxt*(coord[1]+1)-1-NBGY+1)%NSKPY;
    else
      *rec_nedy = NEDY-nyt*coord[1]-1;
    *rec_nyt = (*rec_nedy-*rec_nbgy)/NSKPY+1;
  }

  if(NBGZ > nzt) *rec_nzt = 0;
  else{
    *rec_nbgz = NBGZ-1;  // since rec_nbgz is 0-based
    *rec_nedz = NEDZ-1;
    *rec_nzt = (*rec_nedz-*rec_nbgz)/NSKPZ+1;
  }

  if(*rec_nxt == 0 || *rec_nyt == 0 || *rec_nzt == 0){
    *rec_nxt = 0;
    *rec_nyt = 0;
    *rec_nzt = 0;
  }

  // displacement assumes NPZ=1!
  *displacement *= sizeof(float);

  return;
}
