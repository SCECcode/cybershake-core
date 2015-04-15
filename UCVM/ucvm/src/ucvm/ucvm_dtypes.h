#ifndef UCVM_DTYPES_H
#define UCVM_DTYPES_H

/* Maximum number of supported models */
#define UCVM_MAX_MODELS 10

/* Maximum projection description length */
#define UCVM_MAX_PROJ_LEN 256

/* Maximum model label length */
#define UCVM_MAX_LABEL_LEN 64

/* Maximum model version length */
#define UCVM_MAX_VERSION_LEN 64

/* Maximum path length */
#define UCVM_MAX_PATH_LEN 256

/* NO CVM flag */
#define UCVM_MODEL_NONE -1


/* Byte order */
typedef enum { UCVM_BYTEORDER_LSB = 0, 
               UCVM_BYTEORDER_MSB } ucvm_byteorder_t;


/* Supported error codes */
typedef enum { UCVM_CODE_SUCCESS = 0, 
	       UCVM_CODE_ERROR,
	       UCVM_CODE_DATAGAP,
	       UCVM_CODE_NODATA} ucvm_code_t;


/* Supported coordinate query modes */
typedef enum { UCVM_COORD_GEO_DEPTH = 0, 
	       UCVM_COORD_UTM_DEPTH,
	       UCVM_COORD_GEO_ELEV,
	       UCVM_COORD_UTM_ELEV } ucvm_ctype_t;


/* Supported grid types */
typedef enum { UCVM_GRID_CELL_CENTER = 0, 
	       UCVM_GRID_CELL_VERTEX} ucvm_gtype_t;


/* For swapping endian */
typedef union ucvm_fdata_t {
    float f;
    unsigned char b[4];
} ucvm_fdata_t;


/* 3D point */
typedef struct ucvm_point_t 
{
  double coord[3];
} ucvm_point_t;


/* 3D dims */
typedef struct ucvm_dim_t 
{
  int dim[3];
} ucvm_dim_t;


/* Material properties */
typedef struct ucvm_prop_t 
{
  int model;
  double vp;
  double vs;
  double rho;
} ucvm_prop_t;


/* Region box */
typedef struct ucvm_region_t 
{
  ucvm_ctype_t cmode;
  double p1[3];
  double p2[3];
} ucvm_region_t;


/* Projection */
typedef struct ucvm_proj_t 
{
  char proj[UCVM_MAX_PROJ_LEN];
} ucvm_proj_t;


/* Projection transformation */
typedef struct ucvm_trans_t 
{
  double origin[3];
  double rotate;
  double translate[3];
  ucvm_gtype_t gtype;
} ucvm_trans_t;


/* Model */
typedef struct ucvm_model_t 
{
  char label[UCVM_MAX_LABEL_LEN];
  ucvm_region_t region;
  char config[UCVM_MAX_PATH_LEN];
  int (*init)(int id, ucvm_region_t *r, const char *config);
  int (*finalize)();
  const char* (*version)();
  int (*query)(ucvm_ctype_t cmode,
	       int n, ucvm_point_t *pnt, 
	       ucvm_prop_t *prop);
} ucvm_model_t;


#endif
