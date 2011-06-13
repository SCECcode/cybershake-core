#ifndef UCVM_PROJ_BILINEAR_H
#define UCVM_PROJ_BILINEAR_H


/* Bilinear parameters */
typedef struct ucvm_bilinear_t 
{
  double xi[4];
  double yi[4];
  double dims[2];
} ucvm_bilinear_t;


/* Convert lon,lat to x,y */
int bilinear_geo2xy(ucvm_bilinear_t *p,
		    double lon, double lat, double *rx, double *ry);


/* Convert x,y to lon,lat */
int bilinear_xy2geo(ucvm_bilinear_t *p,
		    double x, double y, double *lon, double *lat);



#endif
