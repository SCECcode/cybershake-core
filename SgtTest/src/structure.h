#define STATCHAR 8
#define COMPCHAR 4
#define TITLCHAR 64

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

struct complex
   {
   float re;
   float im;
   };

struct gftable
   {
   char filename[128];
   float xp;
   float yp;
   float zp;
   float rng;
   float tst;
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
   long long indx; /* indx= xsgt*10000000 + ysgt*1000 + zsgt */
   int xsgt;     /* x grid location */
   int ysgt;     /* y grid location */
   int zsgt;     /* z grid location */
   float h;         /* grid spacing                                         */
   };

struct sgtheader114    /* sgt header for v1.14 */
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

struct sgtheaderout
   {
   int swap_flag;/* set to -12345 -> can check for needed byte-swapping  */
   float olon;   /* lon. of observation point (source in recip. exp.)    */
   float olat;   /* lat. of observation point (source in recip. exp.)    */
   float slon;   /* lon. of source point (observation in recip. exp.)    */
   float slat;   /* lat. of source point (observation in recip. exp.)    */
   float north;  /* northing (km) of observation pt. relative to source  */
   float east;   /* easting (km) of observation pt. relative to source   */
   float depth;  /* depth (km) of source point                           */
   float xazim;  /* azimuth of X-axis in FD model (clockwise from north) */
   float lam;    /* lambda [in dyne/(cm*cm)] at source point             */
   float mu;     /* rigidity [in dyne/(cm*cm)] at source point           */
   float rho;    /* density [in gm/(cm*cm*cm)] at source point           */
   float stime;  /* approximate S travel time source to observation      */
   float xmom;   /* moment strength of x-oriented force in this run      */
   float ymom;   /* moment strength of y-oriented force in this run      */
   float zmom;   /* moment strength of z-oriented force in this run      */
   float tst;    /* start time of 1st point in GF                        */
   int nt;       /* number of time points                                */
   float dt;     /* time sampling                                        */
   };

struct sgtheader
   {
   int nt;          /* number of time points                                */
   float dt;        /* time sampling                                        */
   float tst;       /* start time of 1st point in GF                        */
   float h;         /* grid spacing                                         */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   float lam;       /* lambda [in dyne/(cm*cm)] at output point             */
   float mu;        /* rigidity [in dyne/(cm*cm)] at output point           */
   float rho;       /* density [in gm/(cm*cm*cm)] at output point           */
   int xsrc;        /* x grid location for source (station in recip. exp.)  */
   int ysrc;        /* y grid location for source (station in recip. exp.)  */
   int zsrc;        /* z grid location for source (station in recip. exp.)  */
   int xsta;        /* x grid location for output (source in recip. exp.)   */
   int ysta;        /* y grid location for output (source in recip. exp.)   */
   int zsta;        /* z grid location for output (source in recip. exp.)   */
   float xmom;      /* moment strength of x-oriented force in this run      */
   float ymom;      /* moment strength of y-oriented force in this run      */
   float zmom;      /* moment strength of z-oriented force in this run      */
   };

struct sgtheaderOLD
   {
   int nt;          /* number of time points                                */
   float dt;        /* time sampling                                        */
   float h;         /* grid spacing                                         */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   float lam;       /* lambda [in dyne/(cm*cm)] at output point             */
   float mu;        /* rigidity [in dyne/(cm*cm)] at output point           */
   float rho;       /* density [in gm/(cm*cm*cm)] at output point           */
   int xsrc;        /* x grid location for source (station in recip. exp.)  */
   int ysrc;        /* y grid location for source (station in recip. exp.)  */
   int zsrc;        /* z grid location for source (station in recip. exp.)  */
   int xsta;        /* x grid location for output (source in recip. exp.)   */
   int ysta;        /* y grid location for output (source in recip. exp.)   */
   int zsta;        /* z grid location for output (source in recip. exp.)   */
   float xmom;      /* moment strength of x-oriented force in this run      */
   float ymom;      /* moment strength of y-oriented force in this run      */
   float zmom;      /* moment strength of z-oriented force in this run      */
   };

struct mtheader    /* header for moment tensor output information */
   {
   char title[128];
   int indx;
   int nt;
   float dt;
   float h;
   float modelrot;   /* rotation of y-axis from south (clockwise positive) */
   float modellat;
   float modellon;
   float mu;
   int xsrc;
   int ysrc;
   int zsrc;
   int xsta;
   int ysta;
   int zsta;
   float xmom;
   float ymom;
   float zmom;
   };
