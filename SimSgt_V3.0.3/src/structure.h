#include "defs.h"
 
struct nodeinfo
   {
   int nproc;
   int nproc_x;
   int nproc_y;
   int nproc_z;
   int min_nproc;
   int segmentId;
   int minusId_x;
   int plusId_x;
   int minusId_y;
   int plusId_y;
   int minusId_z;
   int plusId_z;
   int procId_x;
   int procId_y;
   int procId_z;
   int globnx;
   int nx1;
   int nx2;
   int ixminus;
   int ixplus;
   int loc_nx;
   int globny;
   int ny1;
   int ny2;
   int iyminus;
   int iyplus;
   int loc_ny;
   int globnz;
   int nz1;
   int nz2;
   int izminus;
   int izplus;
   int loc_nz;
   char procname[512];
   };
 
struct runparamsP3
   {
   struct nodeinfo ni;
   int nx;
   int ny;
   int nz;
   int nt;
   int freesurf;
   float h;
   float dt;
   int geoproj;
   float modelrot;
   float modellon;
   float modellat;
   int center_origin;
   float xshift;
   float yshift;
   float kmlon;
   float kmlat;
   float cosR;
   float sinR;
   float erad;
   float xmom;
   float ymom;
   float zmom;
   float tdelay;
   double amat[9];
   double ainv[9];
   double g0;
   double b0;
   };
 
struct runparams
   {
   int segmentId;
   int nodeType;
   int nx;
   int globny;
   int nz;
   int nt;
   int ny1;
   int ny2;
   int span;
   int freesurf;
   float h;
   float dt;
   int geoproj;
   float modelrot;
   float modellon;
   float modellat;
   int center_origin;
   float xshift;
   float yshift;
   float kmlon;
   float kmlat;
   float cosR;
   float sinR;
   float erad;
   float xmom;
   float ymom;
   float zmom;
   float tdelay;
   double amat[9];
   double ainv[9];
   double g0;
   double b0;
   };
 
struct restartinfo
   {
   int enable_flag;
   int read_flag;
   int itinc;
   int it;
   int fdr;
   int fdw;
   struct runparams rpars;
   char dir[DIR_STR_LEN];
   char name[FILE_STR_LEN];
   char readfile[FILE_STR_LEN];
   char dumpfile[FILE_STR_LEN];
   char bkupfile[FILE_STR_LEN];
   };
 
struct restartinfoP3
   {
   int enable_flag;
   int read_flag;
   int itinc;
   int it;
   int fdr;
   int fdw;
   struct runparamsP3 rpars_p3;
   char dir[DIR_STR_LEN];
   char name[FILE_STR_LEN];
   char readfile[FILE_STR_LEN];
   char dumpfile[FILE_STR_LEN];
   char bkupfile[FILE_STR_LEN];
   };
 
struct dump_output_info
   {
   int enable_flag;
   int itinc;
   char local_restartdir[DIR_STR_LEN];
   char local_outputdir[DIR_STR_LEN];
   char main_dir[DIR_STR_LEN];
   char command[FILE_STR_LEN];
   char name[FILE_STR_LEN];
   };

struct endian_type     /* global structure for byte-ordering information */
   {
   int num_nodes;
   int server_swap;
   int model_swap;
   int output_swap;
   int left_swap;
   int right_swap;
   char server_endian[ENDIAN_BUF_LEN];
   char my_endian[ENDIAN_BUF_LEN];
   char model_endian[ENDIAN_BUF_LEN];
   char output_endian[ENDIAN_BUF_LEN];
   char *cbuf;
   char **node_endian;
   };

struct seiscords          /* structure for output coordinates */
   {
   int x;
   int y;
   int z;
   };

struct tsheader    /* structure for time slice header information */
   {
   int ix0;          /* starting x grid location for output */
   int iy0;          /* starting y grid location for output */
   int iz0;          /* starting z grid location for output */
   int it0;          /* starting time step for output */
   int nx;          /* number of x points                                */
   int ny;          /* total number of y planes for entire model output */
   int nz;          /* number of z points                                */
   int nt;          /* number of time points                                */
   float dx;         /* grid spacing                                         */
   float dy;         /* grid spacing                                         */
   float dz;         /* grid spacing                                         */
   float dt;        /* time sampling                                        */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   };

struct tsheader_proc  /* structure for individual processor time slice header */
   {
   int ix0;          /* starting x grid location for output */
   int iy0;          /* starting y grid location for output */
   int iz0;          /* starting z grid location for output */
   int it0;          /* starting time step for output */
   int iyleft;          /* global iy of 1st plane in this file */
   int iyright;          /* global iy of last plane in this file */
   int localny;          /* local # of iy planes */
   int nx;          /* number of x points                                */
   int ny;          /* total number of y planes for entire model output */
   int nz;          /* number of z points                                */
   int nt;          /* number of time points                                */
   float dx;         /* grid spacing                                         */
   float dy;         /* grid spacing                                         */
   float dz;         /* grid spacing                                         */
   float dt;        /* time sampling                                        */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   };

struct tsheaderP3    /* structure for P3 time slice header information */
   {
   int loc_np;      /* total of points in local time slice, loc_np=loc_na*loc_nb */
   int loc_na;      /* local number of points in fast direction */
   int loc_nb;      /* local number of points in slow direction */
   int globnp;      /* global total of points for entire time slice, globnp=globna*globnb */
   int globna;      /* global number of points in fast direction */
   int globnb;      /* global number of points in slow direction */
   int na1;         /* global index of first point in fast direction */
   int nb1;         /* global index of first point in slow direction */
   int nt;          /* number of time points in time slice                  */
   float da;        /* grid spacing in fast direction                       */
   float db;        /* grid spacing in slow direction                       */
   float dt;        /* time sampling                                        */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   };
 
struct outputparams       /* structure for output information */
   {
   int iflag;
   struct seiscords pnt;  /* only used for individual points */
   float *v;
   int ix0;
   int iy0;
   int iz0;
   int it0;
   int idx;
   int idy;
   int idz;
   int idt;
   int nx;
   int ny;
   int nz;
   int fdw;
   int iyleft;     /* global iy of first 'responsible' plane in this node */
   int iyright;    /* global iy of last  'responsible' plane in this node */
   off_t head_off;  /* length in bytes of header */
   off_t cur_off;
   };
 
struct tsoutputparams       /* structure for time slice output with buffering */
   {
   int iflag;
   int buflen;
   int nbufrite;
   int flushbuf;
   float *vbuf;
   int ix0;
   int iy0;
   int iz0;
   int it0;
   int idx;
   int idy;
   int idz;
   int idt;
   int nx;
   int ny;
   int nz;
   int np;
   int ntout;
   int fdw;
   int iyleft;     /* global iy of first 'responsible' plane in this node */
   int iyright;    /* global iy of last  'responsible' plane in this node */
   off_t head_off;  /* length in bytes of header */
   off_t cur_off;
   char local_file[FILE_STR_LEN];
   char main_file[FILE_STR_LEN];
   };

struct seisheader
   {
   int indx;        /* numerical index of this location in statcords file */
   int ix;          /* x grid location for output */
   int iy;          /* y grid location for output */
   int iz;          /* z grid location for output */
   int nt;          /* number of time points                                */
   float dt;        /* time sampling                                        */
   float h;         /* grid spacing                                         */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float slat;      /* station latitude (as input from modelcords) */
   float slon;      /* station longitude (as input from modelcords) */
   char name[NBYTE_STATNAME]; /* station name */
   };

struct seisoutputparams  /* structure for time history output to single file */
   {
   int iflag;
   int buflen;
   int nbufrite;
   int flushbuf;
   struct seisheader *shead;
   float *s;
   int np;
   int ntout;
   int tinc;
   int fdw;
   off_t cur_off;
   char local_file[FILE_STR_LEN];
   char main_file[FILE_STR_LEN];
   };
 
struct sgtmaster
   {
   int geoproj;     /* =0: RWG local flat earth; =1: RWG great circle arcs; =2: UTM */
   float modellon;  /* longitude of geographic origin */
   float modellat;  /* latitude of geographic origin */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float xshift;    /* xshift of cartesian origin from geographic origin */
   float yshift;    /* yshift of cartesian origin from geographic origin */
   int globnp;      /* total number of SGT locations (entire model) */
   int localnp;     /* local number of SGT locations (this file only) */
   int nt;          /* number of time points                                */
   };
 
struct sgtindex   /* indices for all 'globnp' SGT locations */
   {
   long long indx; /* indx= xsgt*100000000 + ysgt*10000 + zsgt */
   int xsgt;     /* x grid location */
   int ysgt;     /* y grid location */
   int zsgt;     /* z grid location */
   float h;         /* grid spacing                                         */
   };
 
struct sgtheader
   {
   long long indx;  /* index of this SGT */
   int geoproj;     /* =0: RWG local flat earth; =1: RWG great circle arcs; =2: UTM */
   float modellon;  /* longitude of geographic origin */
   float modellat;  /* latitude of geographic origin */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float xshift;    /* xshift of cartesian origin from geographic origin */
   float yshift;    /* yshift of cartesian origin from geographic origin */
   int nt;          /* number of time points                                */
   float xazim;     /* azimuth of X-axis in FD model (clockwise from north) */
   float dt;        /* time sampling                                        */
   float tst;       /* start time of 1st point in GF                        */
   float h;         /* grid spacing                                         */
   float src_lat;   /* site latitude */
   float src_lon;   /* site longitude */
   float src_dep;   /* site depth */
   int xsrc;        /* x grid location for source (station in recip. exp.)  */
   int ysrc;        /* y grid location for source (station in recip. exp.)  */
   int zsrc;        /* z grid location for source (station in recip. exp.)  */
   float sgt_lat;   /* SGT location latitude */
   float sgt_lon;   /* SGT location longitude */
   float sgt_dep;   /* SGT location depth */
   int xsgt;        /* x grid location for output (source in recip. exp.)   */
   int ysgt;        /* y grid location for output (source in recip. exp.)   */
   int zsgt;        /* z grid location for output (source in recip. exp.)   */
   float cdist;     /* straight-line distance btw site and SGT location */
   float lam;       /* lambda [in dyne/(cm*cm)] at output point             */
   float mu;        /* rigidity [in dyne/(cm*cm)] at output point           */
   float rho;       /* density [in gm/(cm*cm*cm)] at output point           */
   float xmom;      /* moment strength of x-oriented force in this run      */
   float ymom;      /* moment strength of y-oriented force in this run      */
   float zmom;      /* moment strength of z-oriented force in this run      */
   };
 
struct sgtoutputparams    /* structure for strain green tensor output information */
   {
   int iflag;
   int buflen;
   int nbufrite;
   int flushbuf;
   struct sgtmaster sgtmast;
   struct sgtindex *sgtindx;
   struct sgtheader *sgthead;
   float *sgt;
   int tinc;
   int fdw;
   off_t cur_off;
   char local_file[FILE_STR_LEN];
   char main_file[FILE_STR_LEN];
   };
 
struct outputfields       /* structure for output fields */
   {
   char seisdir[DIR_STR_LEN];
   int nseis;
   char seiscords[FILE_STR_LEN];
   struct seisoutputparams spnts;  /* seismos at arbitrary points */
   int ts_xy;
   int ts_xz;
   int ts_yz;
   int ix_ts;
   int iy_ts;
   int iz_ts;
   int dxts;
   int dyts;
   int dzts;
   int dtts;
   struct outputparams xyts;    /* time slices on xy plane */
   struct outputparams xzts;    /* time slices on xz plane */
   struct outputparams yzts;    /* time slices on yz plane */
   int sgtout;
   char sgtcords[256];
   struct sgtoutputparams sgtpnts;  /* strain green tensor */
   struct tsoutputparams xyslice;  /* xy time slice */
   struct tsoutputparams xzslice;  /* xz time slice */
   struct tsoutputparams yzslice;  /* yz time slice */
   };
 
struct statseis
   {
   int fdw;
   int ixs;
   int iys;
   int izs;
   char sname[8];
   };
 
struct medstencil
   {
   char med_name[4];
   int idep0;
   int idep1;
   };
 
struct medprof
   {
   char prof_name[4];
   int ndep;
   int mxdp;
   float *medptr;
   float *dep;
   };
 
struct media_input
   {
   int model_style;
   struct medstencil *medsten;
   struct medprof *medprofs;
   int nprof;
   float xorg;
   float yorg;
   char modfile[128];
   char pmodfile[128];
   char smodfile[128];
   char dmodfile[128];
   char qmodfile[128];
   int media_boundary;
   int dampwidth;
   float qbndmax;
   float qbndmin;
   };
 
struct qvalues
   {
   int n;
   float qpfrac;
   float qsfrac;
   float qpqs_factor;
   float *vp;
   float *qp;
   float *vs;
   float *qs;
   };
 
struct fdcoefs
   {
   int order;
   int ordx;
   float vmin;
   float vmax;
   int izord2;
   float dtoh;
   float c0;
   float c1;
   };

struct doublecouple
   {
   float *strike;
   float *dip;
   float *rake;
   float moment;
   float *momwgts;
   int *itsrc;
   float h;
   };

struct momenttensorOLD
   {
   float *mxx;
   float *myy;
   float *mzz;
   float *mxy;
   float *mxz;
   float *myz;
   };

struct momenttensor
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float rake;
   float area;
   float tinit;
   float slip1;
   float slip2;
   float slip3;
   float dt_stf;
   int nt_stf;
   float *stf1;
   float *stf2;
   float *stf3;
   int nt;
   float *mxx;
   float *myy;
   float *mzz;
   float *mxy;
   float *mxz;
   float *myz;
   };

struct adjointsource
   {
   int ix;
   int iy;
   int iz;
   int nt;
   float *sx;
   float *sy;
   float *sz;
   };

struct pntsrcs
   {
   int psrc;
   int expl;
   int eqsrc;
   int bforce;
   int dblcpl;
   int pointmt;
   int ffault;
   int adjoint;
   float *fxsrc;
   float *fysrc;
   float *fzsrc;
   struct momenttensor *momten;
   struct adjointsource *adjsrc;
   struct doublecouple doublec;
   int nsource;
   int nstf;
   int *stfindx;
   float *rtime;
   float *xs;
   float *ys;
   float *zs;
   int *ix;
   int *iy;
   int *iz;
   int *it;
   int relative_slip;
   int absolute_slip;
   float area;
   float modelrot;
   char faultfile[512];
   char adjointfile[512];
   };
 
struct planesrc
   {
   float incidence;
   float azimuth;
   int intercept;
   int bnd_zero_pad;
   float srcl2m;
   float srclam;
   float srcmu;
   float svel;
   float pvel;
   float p;
   float sh;
   float sv;
   int ixpzero;
   int iypzero;
   int izpzero;
   };

struct complex
   {
   float re;
   float im;
   };

struct sourceparams
   {
   float ts;
   float vxamp;
   float vyamp;
   float vzamp;
   };

struct gfindexpar
   {
   unsigned short igf[8];
   float cfac[8];
   };

struct gfsourcepar
   {
   /*
   int isrc;
   int vxf;
   int vyf;
   int vzf;
   int itgfst;
   int npts;
   */
   unsigned short isrc;
   unsigned short vflag;
   unsigned short itgfst;
   unsigned short npts;
   };

struct gftable
   {
   int ix;
   int iy;
   int iz;
   float *ts;
   int *npts;
   float **gfv;
   };

struct pointsource
   {
   float *xs;
   float *ys;
   float *zs;
   float *rtime;
   float *vsrc;
   float *den;
   float *srcamp;
   float h;
   float dt;
   int expl;
   float strike;
   float dip;
   float rake;
   float cosL;
   float sinL;
   float cosD;
   float sinD;
   float cos2D;
   float sin2D;
   float cosA;
   float sinA;
   float *st;
   int itshft;
   float tstrt;
   float tend;
   int itlen;
   int nsrc;
   int isrc;
   };

struct boundflags
   {
   int x0bnd;
   int xnbnd;
   int y0bnd;
   int ynbnd;
   int z0bnd;
   int znbnd;
   };

struct finitefault
   {
   float xf;
   float yf;
   float zf;
   float moment;
   float strike;
   float dip;
   float length;
   float width;
   int nstrk;
   int ndown;
   int nsrc;
   int npntsrc;
   int npnt2;
   float shypo;
   float dhypo;
   float rvel;
   float tsmin;
   char slipmodel[256];
   float *slip;
   float *rake;
   int *goslip;
   char srcmodel[256];
   int nlay;
   float *vp;
   float *inva2;
   float *vs;
   float *invb2;
   float *vs2;
   float *den;
   float *th;
   struct pointsource *pntsrc;
   };

struct interface
   {
   char gfname[256];
   int go;
   int *goplane;
   int *newgo;
   int mode;
   int ny;
   int iy;
   int yshft;
   int it;
   int nt;
   int min_tspan;
   float *stf;
   float dt;
   struct boundflags *bflags;
   struct finitefault ffault;
   float **gfvbuf;
   int itend;
   char *buf;
   char *buf_int;
   char *buf_ext;
   int ngfsrc;
   int ngfmem;
   float dz;
   float dh;
   struct gftable *gftable;
   };

struct modelstorage
   {
   size_t maxmem;
   int intmem;
   size_t ln_pvslice;
   size_t sz_pvslice;
   int pv_fdr;
   int pv_fdw;
   float *pv_buf;
   float *pv_bufptr;
   size_t ln_medslice;
   size_t sz_medslice;
   int med_fdr;
   int med_fdw;
   float *med_buf;
   float *med_bufptr;
   };
