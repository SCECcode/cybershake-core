#ifndef DEM_DTYPES_H
#define DEM_DTYPES_H


/* Byte order */
typedef enum { DEM_BYTEORDER_LSB = 0, 
	       DEM_BYTEORDER_MSB } dem_byteorder_t;


/* Data source */
typedef enum { DEM_SRC_NED = 0, 
	       DEM_SRC_BATH } dem_src_t;


/* Geo lon/lat point */
typedef struct dem_point_t 
{
  double coord[2];
} dem_point_t;


/* Elevation value */
typedef struct dem_elev_t 
{
  char source[16];
  double elev;
  int valid;
} dem_elev_t;


/* NED file info */
typedef struct dem_info_t 
{
  int dims[2];
  double llcorner[2];
  double urcorner[2];
  double spacing;
  double no_data;
  dem_byteorder_t  byteorder;
  char file[512];
} dem_info_t;


#endif
