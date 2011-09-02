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
