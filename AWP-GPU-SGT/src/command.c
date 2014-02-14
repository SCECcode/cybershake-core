/*
**************************************************************************************************************** 
*  command.c						                                                       *
*  Process Command Line	                                                                                       *
*                                                                                                              *
*  Name         Type        Command             Description 	                                                 *
*  TMAX         <FLOAT>       -T              propagation time	                	                             *
*  DH           <FLOAT>       -H              spatial step for x, y, z (meters)                                *
*  DT           <FLOAT>       -t              time step (seconds)                                              *
*  ARBC         <FLOAT>       -A              oefficient for PML (3-4), or Cerjan (0.90-0.96)                  *
*  PHT          <FLOAT>       -P                                                                               *
*  NPC          <INTEGER>     -M              PML or Cerjan ABC (1=PML, 0=Cerjan)                              *
*  ND           <INTEGER>     -D              ABC thickness (grid-points) PML <= 20, Cerjan >= 20              *
*  NSRC         <INTEGER>     -S              number of source nodes on fault                                  *
*  NST          <INTEGER>     -N              number of time steps in rupture functions                        *
*  NVAR         <INTEGER>     -n              number of variables in a grid point                              *
*  NVE          <INTEGER>     -V              visco or elastic scheme (1=visco, 0=elastic)                     *
*  MEDIASTART   <INTEGER>     -B              initial media restart option(0=homogenous)                       *
*  IFAULT       <INTEGER>     -I              mode selection and fault or initial stress setting (1 or 2)      *
*  READ_STEP    <INTEGER>     -R                                                                               *
*  READ_STEP_GPU<INTEGER>     -Q              CPU reads larger chunks and sends to GPU at every READ_STEP_GPU  *
*                                               (IFAULT=2) READ_STEP must be divisible by READ_STEP_GPU        *
*  NX           <INTEGER>     -X              x model dimension in nodes                                       *      
*  NY           <INTEGER>     -Y              y model dimension in nodes                                       *
*  NZ           <INTEGER>     -Z              z model dimension in nodes                                       *
*  PX           <INTEGER>     -x              number of procs in the x direction                               *
*  PY           <INTEGER>     -y              number of procs in the y direction                               * 
*  NBGX         <INTEGER>                     index (starts with 1) to start recording points in X             *
*  NEDX         <INTEGER>                     index to end recording points in X (-1 for all)                  *
*  NSKPX        <INTEGER>                     #points to skip in recording points in X                         *
*  NBGY         <INTEGER>                     index to start recording points in Y                             *
*  NEDY         <INTEGER>                     index to end recording points in Y (-1 for all)                  *
*  NSKPY        <INTEGER>                     #points to skip in recording points in Y                         *
*  NBGZ         <INTEGER>                     index to start recording points in Z                             *
*  NEDZ         <INTEGER>                     index to end recording points in Z (-1 for all)                  *
*  NSKPZ        <INTEGER>                     #points to skip in recording points in Z                         *
*  IDYNA        <INTEGER>     -i              mode selection of dynamic rupture model                          *
*  SoCalQ       <INTEGER>     -s              Southern California Vp-Vs Q relationship enabling flag           *
*  FL           <FLOAT>       -l              Q bandwidth low frequency                                        *
*  FH           <FLOAT>       -h              Q bandwidth high frequency                                       *
*  FP           <FLOAT>       -p              Q bandwidth central frequency                                    *
*  NTISKP       <INTEGER>     -r              # timesteps to skip to copy velocities from GPU to CPU           *
*  WRITE_STEP   <INTEGER>     -W              # timesteps to write the buffer to the files                     *
*                                               (written timesteps are n*NTISKP*WRITE_STEP for n=1,2,...)      *
*  INSRC        <STRING>                      source input file (if IFAULT=2, then this is prefix of tpsrc)    *
*  INVEL        <STRING>                      mesh input file                                                  *
*  CHKFILE      <STRING>      -c              Checkpoint statistics file to write to                           *
*  IGREEN       <INTEGER>     -G              IGREEN mode for SGT. IGREEN=-1 is AWPODC                         *
*  NTISKP_SGT   <INTEGER>                     NTISKP for SGT recording and outputs                             *
*  INSGT        <STRING>                      SGT input file to read receivers                                 *
****************************************************************************************************************
*/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <string.h>

// Default IN3D Values
const float def_TMAX       = 20.00;
const float def_DH         = 200.0;
const float def_DT         = 0.01; 
const float def_ARBC       = 0.92;
const float def_PHT        = 0.1;

const int   def_NPC        = 0;
const int   def_ND         = 20;
const int   def_NSRC       = 1;
const int   def_NST        = 91;
const int   def_NVAR       = 3;

const int   def_NVE        = 1;
const int   def_MEDIASTART = 0;
const int   def_IFAULT     = 1; 
const int   def_READ_STEP  = 91;
const int   def_READ_STEP_GPU = 91;

const int   def_NTISKP     = 10;
const int   def_WRITE_STEP = 10;

const int   def_NX         = 224;
const int   def_NY         = 224; 
const int   def_NZ         = 1024;

const int   def_PX         = 1;
const int   def_PY         = 1;

const int   def_NBGX       = 1;
const int   def_NEDX       = -1;   // use -1 for all
const int   def_NSKPX      = 1;
const int   def_NBGY       = 1;
const int   def_NEDY       = -1;   // use -1 for all
const int   def_NSKPY      = 1;
const int   def_NBGZ       = 1;
const int   def_NEDZ       = 1;    // only surface
const int   def_NSKPZ      = 1;

const int   def_IDYNA      = 0;
const int   def_SoCalQ     = 1;

const float def_FL         = 0.005;
const float def_FH         = 5.0;
const float def_FP         = 2.5;

const char  def_INSRC[50]  = "input/FAULTPOW";
const char  def_INVEL[50]  = "input/media";

const char  def_OUT[50] = "output_sfc";

const char  def_CHKFILE[50]   = "output_sfc/CHKP";

const int   def_IGREEN     = -1;
const int   def_NTISKP_SGT = 1;
const char  def_INSGT[50]  = "input/SGT";

void command(int argc,    char **argv,
	     float *TMAX, float *DH,       float *DT,   float *ARBC,    float *PHT,
             int *NPC,    int *ND,         int *NSRC,   int *NST,       int *NVAR,
             int *NVE,    int *MEDIASTART, int *IFAULT, int *READ_STEP, int *READ_STEP_GPU,
             int *NTISKP, int *WRITE_STEP,
    	       int *NX,     int *NY,         int *NZ,     int *PX,        int *PY,
             int *NBGX,   int *NEDX,       int *NSKPX, 
             int *NBGY,   int *NEDY,       int *NSKPY, 
             int *NBGZ,   int *NEDZ,       int *NSKPZ, 
             float *FL,   float *FH,       float *FP,   int *IDYNA,     int *SoCalQ,
             char *INSRC, char *INVEL,     char *OUT,   char *CHKFILE,
             int *IGREEN, int *NTISKP_SGT, char *INSGT)
{

   // Fill in default values
   *TMAX       = def_TMAX;
   *DH         = def_DH;
   *DT         = def_DT;
   *ARBC       = def_ARBC;
   *PHT        = def_PHT;

   *NPC        = def_NPC;
   *ND         = def_ND;
   *NSRC       = def_NSRC;
   *NST        = def_NST;
   
   *NVE        = def_NVE;
   *MEDIASTART = def_MEDIASTART;
   *NVAR       = def_NVAR;
   *IFAULT     = def_IFAULT; 
   *READ_STEP  = def_READ_STEP;
   *READ_STEP_GPU = def_READ_STEP_GPU;

   *NTISKP     = def_NTISKP;
   *WRITE_STEP = def_WRITE_STEP;

   *NX         = def_NX;
   *NY         = def_NY;
   *NZ         = def_NZ;
   *PX         = def_PX;
   *PY         = def_PY;

   *NBGX       = def_NBGX;
   *NEDX       = def_NEDX;
   *NSKPX      = def_NSKPX;
   *NBGY       = def_NBGY;
   *NEDY       = def_NEDY;
   *NSKPY      = def_NSKPY;
   *NBGZ       = def_NBGZ;
   *NEDZ       = def_NEDZ;
   *NSKPZ      = def_NSKPZ;

   *IDYNA      = def_IDYNA;
   *SoCalQ     = def_SoCalQ;
   *FL         = def_FL;
   *FH         = def_FH;
   *FP         = def_FP;

    strcpy(INSRC, def_INSRC);
    strcpy(INVEL, def_INVEL);
    strcpy(OUT, def_OUT);
    strcpy(CHKFILE, def_CHKFILE);

   *IGREEN     = def_IGREEN;
   *NTISKP_SGT = def_NTISKP_SGT;
    strcpy(INSGT, def_INSGT);

    extern char *optarg;
    static const char *optstring = "-T:H:t:A:P:M:D:S:N:V:B:n:I:R:Q:X:Y:Z:x:y:z:i:l:h:p:s:r:W:1:2:3:11:12:13:21:22:23:100:101:o:c:G:200:201:";
    static struct option long_options[] = {
        {"TMAX", required_argument, NULL, 'T'},
        {"DH", required_argument, NULL, 'H'},
        {"DT", required_argument, NULL, 't'},
        {"ARBC", required_argument, NULL, 'A'},
        {"PHT", required_argument, NULL, 'P'},
        {"NPC", required_argument, NULL, 'M'},
        {"ND", required_argument, NULL, 'D'},
        {"NSRC", required_argument, NULL, 'S'},
        {"NST", required_argument, NULL, 'N'},
        {"NVE", required_argument, NULL, 'V'},
        {"MEDIASTART", required_argument, NULL, 'B'},
        {"NVAR", required_argument, NULL, 'n'},
        {"IFAULT", required_argument, NULL, 'I'},
        {"READ_STEP", required_argument, NULL, 'R'},
        {"READ_STEP_GPU", required_argument, NULL, 'Q'},
        {"NX", required_argument, NULL, 'X'},
        {"NY", required_argument, NULL, 'Y'},
        {"NZ", required_argument, NULL, 'Z'},
        {"PX", required_argument, NULL, 'x'},
        {"PY", required_argument, NULL, 'y'},
        {"NBGX", required_argument, NULL, 1},
        {"NEDX", required_argument, NULL, 2},
        {"NSKPX", required_argument, NULL, 3},
        {"NBGY", required_argument, NULL, 11},
        {"NEDY", required_argument, NULL, 12},
        {"NSKPY", required_argument, NULL, 13},
        {"NBGZ", required_argument, NULL, 21},
        {"NEDZ", required_argument, NULL, 22},
        {"NSKPZ", required_argument, NULL, 23},
        {"IDYNA", required_argument, NULL, 'i'},
        {"SoCalQ", required_argument, NULL, 's'},
        {"FL", required_argument, NULL, 'l'},
        {"FH", required_argument, NULL, 'h'},
        {"FP", required_argument, NULL, 'p'},
        {"NTISKP", required_argument, NULL, 'r'},
        {"WRITE_STEP", required_argument, NULL, 'W'},
        {"INSRC", required_argument, NULL, 100},
        {"INVEL", required_argument, NULL, 101},
        {"OUT", required_argument, NULL, 'o'},
        {"CHKFILE", required_argument, NULL, 'c'},
        {"IGREEN", required_argument, NULL, 'G'},
        {"NTISKP_SGT", required_argument, NULL, 200},
        {"INSGT", required_argument, NULL, 201},
    };

    // If IFAULT=1 and READ_STEP_GPU is not set, it should be = READ_STEP
    int readstepGpuIsSet = 0;
    int c;
    while ((c=getopt_long(argc, argv, optstring, long_options, NULL)) != -1)
    {
        switch (c) {
            case 'T':
                *TMAX       = atof(optarg); break;
            case 'H':
                *DH         = atof(optarg); break;
            case 't':
                *DT         = atof(optarg); break;
            case 'A':
                *ARBC       = atof(optarg); break;
            case 'P':
                *PHT        = atof(optarg); break;
            case 'M':
                *NPC        = atoi(optarg); break;
            case 'D':
                *ND         = atoi(optarg); break;
            case 'S':
                *NSRC       = atoi(optarg); break;
            case 'N':
                *NST        = atoi(optarg); break;
            case 'V':
                *NVE        = atoi(optarg); break;
            case 'B':
                *MEDIASTART = atoi(optarg); break;
            case 'n':
                *NVAR       = atoi(optarg); break;
            case 'I':
	        *IFAULT     = atoi(optarg); break;
            case 'R':
                *READ_STEP  = atoi(optarg); break;		
            case 'Q':
                readstepGpuIsSet = 1;
                *READ_STEP_GPU  = atoi(optarg); break;		
            case 'X':
                *NX         = atoi(optarg); break;
            case 'Y':
                *NY         = atoi(optarg); break;
            case 'Z':
                *NZ         = atoi(optarg); break;
            case 'x':
                *PX         = atoi(optarg); break;
            case 'y':
                *PY         = atoi(optarg); break;
            case 1:
                *NBGX       = atoi(optarg); break;
            case 2:
                *NEDX       = atoi(optarg); break;
            case 3:
                *NSKPX      = atoi(optarg); break;
            case 11:
                *NBGY       = atoi(optarg); break;
            case 12:
                *NEDY       = atoi(optarg); break;
            case 13:
                *NSKPY      = atoi(optarg); break;
            case 21:
                *NBGZ       = atoi(optarg); break;
            case 22:
                *NEDZ       = atoi(optarg); break;
            case 23:
                *NSKPZ      = atoi(optarg); break;
            case 'i':
                *IDYNA      = atoi(optarg); break;
            case 's':
                *SoCalQ     = atoi(optarg); break;
            case 'l':
                *FL         = atof(optarg); break;
            case 'h':
                *FH         = atof(optarg); break;
            case 'p':
                *FP         = atof(optarg); break;
            case 'r':
                *NTISKP     = atoi(optarg); break;
            case 'W':
                *WRITE_STEP = atoi(optarg); break;
            case 100:
                strcpy(INSRC, optarg); break;
            case 101:
                strcpy(INVEL, optarg); break;
            case 'o':
                strcpy(OUT, optarg); break;
            case 'c':
                strcpy(CHKFILE, optarg); break;
            case 'G':
                *IGREEN     = atoi(optarg); break;
            case 200:
                *NTISKP_SGT = atoi(optarg); break;
            case 201:
                strcpy(INSGT, optarg); break;
            default:
                printf("Usage: %s \nOptions:\n\t[(-T | --TMAX) <TMAX>]\n\t[(-H | --DH) <DH>]\n\t[(-t | --DT) <DT>]\n\t[(-A | --ARBC) <ARBC>]\n\t[(-P | --PHT) <PHT>]\n\t[(-M | --NPC) <NPC>]\n\t[(-D | --ND) <ND>]\n\t[(-S | --NSRC) <NSRC>]\n\t[(-N | --NST) <NST>]\n",argv[0]);
                printf("\n\t[(-V | --NVE) <NVE>]\n\t[(-B | --MEDIASTART) <MEDIASTART>]\n\t[(-n | --NVAR) <NVAR>]\n\t[(-I | --IFAULT) <IFAULT>]\n\t[(-R | --READ_STEP) <x READ_STEP]\n");
                printf("\n\t[(-X | --NX) <x length]\n\t[(-Y | --NY) <y length>]\n\t[(-Z | --NZ) <z length]\n\t[(-x | --NPX) <x processors]\n\t[(-y | --NPY) <y processors>]\n\t[(-z | --NPZ) <z processors>]\n");
                printf("\n\t[(-1 | --NBGX) <starting point to record in X>]\n\t[(-2 | --NEDX) <ending point to record in X>]\n\t[(-3 | --NSKPX) <skipping points to record in X>]\n\t[(-11 | --NBGY) <starting point to record in Y>]\n\t[(-12 | --NEDY) <ending point to record in Y>]\n\t[(-13 | --NSKPY) <skipping points to record in Y>]\n\t[(-21 | --NBGZ) <starting point to record in Z>]\n\t[(-22 | --NEDZ) <ending point to record in Z>]\n\t[(-23 | --NSKPZ) <skipping points to record in Z>]\n");
                printf("\n\t[(-i | --IDYNA) <i IDYNA>]\n\t[(-s | --SoCalQ) <s SoCalQ>]\n\t[(-l | --FL) <l FL>]\n\t[(-h | --FH) <i FH>]\n\t[(-p | --FP) <p FP>]\n\t[(-r | --NTISKP) <time skipping in writing>]\n\t[(-W | --WRITE_STEP) <time aggregation in writing>]\n");
                printf("\n\t[(-100 | --INSRC) <source file>]\n\t[(-101 | --INVEL) <mesh file>]\n\t[(-o | --OUT) <output file>]\n\t[(-c | --CHKFILE) <checkpoint file to write statistics>]\n\n");
                printf("\n\t[(-G | --IGREEN) <IGREEN for SGT>]\n\t[(-200 | --NTISKP_SGT) <NTISKP for SGT>]\n\t[(-201 | --INSGT) <SGT input file>]\n\n");
                exit(-1);
        }
    }
    if(!readstepGpuIsSet){
      *READ_STEP_GPU = *READ_STEP;
    }
    return;
}
