#ifndef UCVM_UTILS_H
#define UCVM_UTILS_H

#include "ucvm_dtypes.h"

/* Determine system endian */
int system_endian();

/* Swap float endian */
float swap_endian_float(float f);

/* Parses string list into double array */
int list_parse(char *lstr, int n, double *arr);

/* Returns true if region contains point p of coord type c */
int region_contains(ucvm_region_t *r, ucvm_ctype_t c, ucvm_point_t *p);

/* Parses string region into structure */
int region_parse(char *rstr, ucvm_region_t *r);

/* Parses region structure into string */
int region_string(ucvm_region_t *r, char *rstr, int len);

/* Rotate point in 2d about origin by theta degrees */
int rot_point_2d(ucvm_point_t *p, ucvm_point_t *o, double theta);

/* Interpolate point linearly between two 1d values */
double interpolate_linear_1d(double v1, double v2, double ratio);

/* Interpolate point bilinearly between four corners */
double interpolate_bilinear_2d(double x, double y, 
			       double x1, double y1, double x2, double y2, 
			       double q11, double q21, double q12, double q22);

/* Interpolate point tri-linearly between 8 cube corners.
   Points are indexed [ll,ur][x,y,z], q is indexed[z][y][x] */
double interpolate_trilinear_1d(double x, double y, double z,
				double p[2][3], double q[2][2][2]);

#endif
