/*
 * structure.h
 *
 *  Created on: Dec 12, 2014
 *      Author: scott
 */

#ifndef STRUCTURE_H_
#define STRUCTURE_H_

struct statdata
   {
   char stat[STATCHAR];
   char comp[COMPCHAR];
   char stitle[TITLCHAR];
   int nt;
   float dt;
   int hr;
   int min;
   float sec;
   float edist;
   float az;
   float baz;
   };


struct geoprojection
   {
   int geoproj;
   float modellon;
   float modellat;
   float modelrot;
   float xshift;
   float yshift;
   int center_origin;
   double rperd;
   float erad;
   float fc;
   float g2;
   float radc;
   float cosR;
   float sinR;
   float kmlon;
   float kmlat;
   double g0;
   double b0;
   double amat[9];
   double ainv[9];
   };

struct mechparam
   {
   int nmech;
   int flag[3];
   float stk;
   float dip;
   float rak;
   };

struct sgtfileparams
   {
   FILE* xfp;
   FILE* yfp;
   FILE* zfp;
   char xfile[256];
   char yfile[256];
   char zfile[256];
   off_t head_off;
   off_t xcur_off;
   off_t ycur_off;
   off_t zcur_off;
   char xfile_header[256];
   char yfile_header[256];
   char zfile_header[256];
   FILE* x_head_fp;
   FILE* y_head_fp;
   FILE* z_head_fp;
   off_t x_head_off;
   off_t y_head_off;
   off_t z_head_off;
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

/*Scott modified on 6/7/16 to allot more digits per direction */
struct sgtindex   /* indices for all 'globnp' SGT locations */
   {
   long long indx; /* indx= xsgt*1,000,000,000 + ysgt*10,000 + zsgt */
   int xsgt;     /* x grid location */
   int ysgt;     /* y grid location */
   int zsgt;     /* z grid location */
   float h;         /* grid spacing                                         */
   };

struct sgtheader    /* sgt header for v1.14 */
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

struct sgtparams
   {
   int nsgt;
   long long indx[4];
   float wt[4];
   int master_ip[4];
   //Pointer to the sgtheaders for processing
   struct sgtheader* sgt_head_ptrs[4];
   //Pointer to the SGTs in sgtbuf
   float* sgtbuf_ptrs[4];
   };


//All the information needed to regenerate a seismogram
//For version 12.10, try avoiding CyberShake-specific data
struct seisheader {
  char version[8];
  char site_name[8];
  //in case we think of something later
  char padding[8];
  int source_id;
  int rupture_id;
  int rup_var_id;
  float dt;
  int nt;
  int comps;
  float det_max_freq;
  float stoch_max_freq;
};


struct rotD_record {
  float period;
  float rotD100;
  int rotD100_angle;
  float rotD50;
};

struct vertical_rsp_record {
  float period;
  float rsp;
};

struct rup_geom_point {
  float lon;
  float lat;
  float dep;
};

typedef struct master_message
	{
	int msg_src;
	int msg_type;
	int msg_data;
} master_msg;

typedef struct data_file_metadata {
	int src_id;
	int rup_id;
	int starting_rv;
	//exclusive
	int ending_rv;
	char filename[256];
} data_file_metadata;

typedef struct fp_cache_entry {
	char filename[256];
	FILE* fp;
} fp_cache_entry;

typedef struct handler_message {
	int msg_src;
	int msg_type;
	int num_sgts_requested;
} handler_msg;

typedef struct rupture_variation {
  int rup_var_id;
  int slip_id;
  int hypo_id;
} rup_var_struct;

typedef struct worker_task {
	char rupture_filename[256];
	int source_id;
	int rupture_id;
	int num_slips;
	int num_hypos;
	int starting_rv_id;
	int ending_rv_id;
	int rup_var_spacing;
} worker_task;

typedef struct manager_message {
	int msg_src;
	int msg_type;
} manager_msg;

typedef struct worker_message {
	int msg_type;
	worker_task task;
} worker_msg;

typedef struct task_info {
	//Assemble all the info we need to work with a task, so we can pass less around
	int argc;
	char** argv;
	struct sgtmaster* sgtmast;
	struct sgtindex* sgtindx;
	struct sgtfileparams* sgtfilepar;
	worker_task* task;
} task_info;

typedef struct rv_info {
	int source_id;
	int rupture_id;
	int rup_var_id;
	float rvfrac;
	int seed;
} rv_info;

#endif /* STRUCTURE_H_ */
