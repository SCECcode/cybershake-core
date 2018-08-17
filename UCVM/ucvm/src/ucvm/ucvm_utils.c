#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ucvm_utils.h"


/* Region coord delimiter */
#define REGION_DELIM ","

/* List delimiter */
#define LIST_DELIM ","


/* Determine system endian */
int system_endian()
{
  int num = 1;
  if(*(char *)&num == 1) {
    return UCVM_BYTEORDER_LSB;
  } else {
    return UCVM_BYTEORDER_MSB;
  }
}


/* Swap float endian-ness */
float swap_endian_float(float f)
{
  ucvm_fdata_t dat1, dat2;

  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return(dat2.f);
}


/* Parses string list into double array */
int list_parse(char *lstr, int n, double *arr)
{
  char *token;
  int i = 0;

  if ((lstr == NULL) || (n <= 0) || (arr == NULL)) {
    return(UCVM_CODE_ERROR);
  }

  token = strtok(lstr, LIST_DELIM);
  while ((token != NULL) && (i < n)) {
    arr[i] = atof(token);
    i++;
    token = strtok(NULL, LIST_DELIM);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Returns true if region contains point p of coord type c */
int region_contains(ucvm_region_t *r, ucvm_ctype_t c, ucvm_point_t *p)
{
  /* NOTE remove this */
  return(1);

  if (c != r->cmode) {
    fprintf(stderr, "Coord type conversion not supported in region_contains().\n");
    return(0);
  }

  if ((p->coord[0] < r->p1[0]) || (p->coord[0] > r->p2[0])) {
    return(0);
  }
  if ((p->coord[1] < r->p1[1]) || (p->coord[1] > r->p2[1])) {
    return(0);
  }
  return(1);
}


/* Parses string region specification into structure struct */
int region_parse(char *rstr, ucvm_region_t *r)
{
  char *token;
  int i = 0;

  r->p1[2] = 0.0;
  r->p2[2] = 0.0;

  token = strtok(rstr, REGION_DELIM);
  while ((token != NULL) && (i < 4)) {
    switch(i) {
    case 0:
      r->p1[0] = atof(token);
      break;
    case 1:
      r->p1[1] = atof(token);
      break;
    case 2:
      r->p2[0] = atof(token);
      break;
    case 3:
      r->p2[1] = atof(token);
      return(UCVM_CODE_SUCCESS);
      break;
    }
    i++;
    token = strtok(NULL, REGION_DELIM);
  }

  return(UCVM_CODE_ERROR);
}

/* Parses region structure into string*/
int region_string(ucvm_region_t *r, char *rstr, int len)
{
  snprintf(rstr, len, "lon=%lf, lat=%lf to lon=%lf, lat=%lf",
	   r->p1[0], r->p1[1], r->p2[0], r->p2[1]);
  return(UCVM_CODE_SUCCESS);
}


/* Rotate point in 2d about origin by theta radians */
int rot_point_2d(ucvm_point_t *p, ucvm_point_t *o, double theta)
{
  double x_offset, y_offset;

  /* Compute offset from origin */
  x_offset = p->coord[0] - o->coord[0];
  y_offset = p->coord[1] - o->coord[1];

  /* Rotate this offset */
  //printf("x_offset=%lf, y_offset=%lf\n", x_offset, y_offset);
  p->coord[0] = (x_offset) * cos(theta) - (y_offset) * sin(theta);
  p->coord[1] = (x_offset) * sin(theta) + (y_offset) * cos(theta);

  /* Add origin back in */
  p->coord[0] = p->coord[0] + o->coord[0];
  p->coord[1] = p->coord[1] + o->coord[1];

  return(0);
}

/* Interpolate point linearly between two 1d values */
double interpolate_linear_1d(double v1, double v2, double ratio) 
{
  return(ratio*v2 + v1*(1-ratio));
}



/* Interpolate point bilinearly between four corners */
double interpolate_bilinear_2d(double x, double y, 
			       double x1, double y1, double x2, double y2, 
			       double q11, double q21, double q12, double q22)
{
  double p = (x2 - x1) * (y2 - y1);
  double f1 = (q11 / p) * (x2 - x) * (y2 - y);
  double f2 = (q21 / p) * (x - x1) * (y2 - y);
  double f3 = (q12 / p) * (x2 - x) * (y - y1);
  double f4 = (q22 / p) * (x - x1) * (y - y1);
  return f1 + f2 + f3 + f4;
}


/* Interpolate point tri-linearly between 8 cube corners.
   Points are indexed [ll,ur][x,y,z], q is indexed[z][y][x] */
double interpolate_trilinear_1d(double x, double y, double z,
			      double p[2][3], double q[2][2][2]) 
{
  double c0, c1;
  double ratio;

  /* Top plane */
  c0 = interpolate_bilinear_2d(x, y,
			       p[0][0], p[0][1],
			       p[1][0], p[1][1],
			       q[0][0][0], q[0][0][1], 
			       q[0][1][0], q[0][1][1]);

  /* Bottom plane */
  c1 = interpolate_bilinear_2d(x, y,
			       p[0][0], p[0][1],
			       p[1][0], p[1][1],
			       q[1][0][0], q[1][0][1], 
			       q[1][1][0], q[1][1][1]);

  /* Z axis */
  ratio = (z - p[0][2])/(p[1][2] - p[0][2]); 
  return(interpolate_linear_1d(c0, c1, ratio));
}
