#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ucvm_utils.h"
#include "ucvm_config.h"
#include "ucvm_proj_bilinear.h"
#include "ucvm_model_linthurber.h"


/* Version string */
char ucvm_linthurber_version_id[UCVM_CONFIG_MAX_STR];

/* Label string */
const char *ucvm_linthurber_label_id = "lin-thurber";

/* Model ID */
int ucvm_linthurber_id = UCVM_MODEL_NONE;

/* Lin-Thurber projection parameters */
ucvm_bilinear_t linthurber_proj;

/* Model dimensions */
#define LINTHURBER_MAX_Z_DIM 100
int linthurber_z_dim = 0;
double linthurber_depths_msl[LINTHURBER_MAX_Z_DIM];
double linthurber_vp_spacing = 0.0;
double linthurber_vs_spacing = 0.0;
double linthurber_dem_spacing = 0.0;
int linthurber_vp_dims[3];
int linthurber_vs_dims[3];
int linthurber_dem_dims[3];
double linthurber_vp_origin[2];
double linthurber_vs_origin[2];


/* Model volumes */
float *vp_buf = NULL;
float *vs_buf = NULL;
float *dem_buf = NULL;


/* Property constants */
#define LINTHURBER_DEM 0
#define LINTHURBER_VP 1
#define LINTHURBER_VS 2


/* Flags */
#define LINTHURBER_USE_DEM 1


/* Init LinThurber */
int ucvm_linthurber_model_init(int m, ucvm_region_t *r, const char *mpath)
{
  int i, j, k;
  FILE *fp;
  char filename[UCVM_MAX_PATH_LEN];
  float dep, x, y, val;
  int num_read, retval;

  ucvm_config_t *cfg = NULL;
  ucvm_config_t *cfgentry = NULL;

  ucvm_linthurber_id = m;

  /* Read CONF file */
  sprintf(filename, "%s/lin-thurber.conf", mpath);
  cfg = ucvm_parse_config(filename);
  if (cfg == NULL) {
    fprintf(stderr, "Failed to read Lin-Thurber config file\n");
    return(UCVM_CODE_ERROR);
  }

  cfgentry = ucvm_find_name(cfg, "version");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find version key\n");
    return(UCVM_CODE_ERROR);
  }
  strcpy(ucvm_linthurber_version_id, cfgentry->value);

  cfgentry = ucvm_find_name(cfg, "proj_xi");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find proj_xi key\n");
    return(UCVM_CODE_ERROR);
  }
  if (list_parse(cfgentry->value, 4, 
		 linthurber_proj.xi) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to parse %s value\n", cfgentry->name);
    return(UCVM_CODE_ERROR);
  }

  cfgentry = ucvm_find_name(cfg, "proj_yi");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find proj_yi key\n");
    return(UCVM_CODE_ERROR);
  }
  if (list_parse(cfgentry->value, 4, 
		 linthurber_proj.yi) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to parse %s value\n", cfgentry->name);
    return(UCVM_CODE_ERROR);
  }

  cfgentry = ucvm_find_name(cfg, "proj_size");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find proj_size key\n");
    return(UCVM_CODE_ERROR);
  }
  if (list_parse(cfgentry->value, 2, 
		 linthurber_proj.dims) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to parse %s value\n", cfgentry->name);
    return(UCVM_CODE_ERROR);
  }

  cfgentry = ucvm_find_name(cfg, "spacing_dem");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find spacing_dem key\n");
    return(UCVM_CODE_ERROR);
  }
  linthurber_dem_spacing = atof(cfgentry->value);

  cfgentry = ucvm_find_name(cfg, "spacing_vp");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find spacing_vp key\n");
    return(UCVM_CODE_ERROR);
  }
  linthurber_vp_spacing = atof(cfgentry->value);

  cfgentry = ucvm_find_name(cfg, "spacing_vs");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find spacing_vs key\n");
    return(UCVM_CODE_ERROR);
  }
  linthurber_vs_spacing = atof(cfgentry->value);

  cfgentry = ucvm_find_name(cfg, "num_z");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find num_z key\n");
    return(UCVM_CODE_ERROR);
  }
  linthurber_z_dim = atoi(cfgentry->value);

  cfgentry = ucvm_find_name(cfg, "z_vals");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find z_vals key\n");
    return(UCVM_CODE_ERROR);
  }
  if (list_parse(cfgentry->value, linthurber_z_dim, 
		 linthurber_depths_msl) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to parse %s value\n", cfgentry->name);
    return(UCVM_CODE_ERROR);
  }

  cfgentry = ucvm_find_name(cfg, "grid_origin");
  if (cfgentry == NULL) {
    fprintf(stderr, "Failed to find grid_origin key\n");
    return(UCVM_CODE_ERROR);
  }
  if (list_parse(cfgentry->value, 2, 
		 linthurber_vp_origin) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to parse %s value\n", cfgentry->name);
    return(UCVM_CODE_ERROR);
  }
  linthurber_vs_origin[0] = linthurber_vp_origin[0];
  linthurber_vs_origin[1] = linthurber_vp_origin[1];

  ucvm_free_config(cfg);

  /* Calculate model dims */
  linthurber_dem_dims[0] = linthurber_proj.dims[0] / linthurber_dem_spacing + 1;
  linthurber_dem_dims[1] = linthurber_proj.dims[1] / linthurber_dem_spacing + 1;
  linthurber_dem_dims[2] = 1;

  linthurber_vp_dims[0] = linthurber_proj.dims[0] / linthurber_vp_spacing + 1;
  linthurber_vp_dims[1] = linthurber_proj.dims[1] / linthurber_vp_spacing + 1;
  linthurber_vp_dims[2] = linthurber_z_dim;

  linthurber_vs_dims[0] = linthurber_proj.dims[0] / linthurber_vs_spacing + 1;
  linthurber_vs_dims[1] = linthurber_proj.dims[1] / linthurber_vs_spacing + 1;
  linthurber_vs_dims[2] = linthurber_z_dim;

  //printf("Calculated:\n");
  //printf("\tDEM Dims: %d, %d, %d\n", linthurber_dem_dims[0],
  //	 linthurber_dem_dims[1],
  //	 linthurber_dem_dims[2]);
  //printf("\tVp Dims: %d, %d, %d\n", linthurber_vp_dims[0],
  //	 linthurber_vp_dims[1],
  //	 linthurber_vp_dims[2]);
  //printf("\tVs Dims: %d, %d, %d\n", linthurber_vs_dims[0],
  //	 linthurber_vs_dims[1],
  //	 linthurber_vs_dims[2]);

  /* Allocate buffers */
  vp_buf = malloc(linthurber_vp_dims[0] * linthurber_vp_dims[1] *
		  linthurber_vp_dims[2] * sizeof(float));
  vs_buf = malloc(linthurber_vs_dims[0] * linthurber_vs_dims[1] *
		  linthurber_vs_dims[2] * sizeof(float));
  dem_buf = malloc(linthurber_dem_dims[0] * linthurber_dem_dims[1] *
		   linthurber_dem_dims[2] * sizeof(float));
  if ((vp_buf == NULL) || (vs_buf == NULL) || (dem_buf == NULL)) {
    fprintf(stderr, "Failed to allocate buffers Lin-Thurber model\n");
    return(UCVM_CODE_ERROR);
  }

  for (i = 0; i < linthurber_vp_dims[0] * linthurber_vp_dims[1] *
	 linthurber_vp_dims[2]; i++) {
    vp_buf[i] = -1.0;
  }
  for (i = 0; i < linthurber_vs_dims[0] * linthurber_vs_dims[1] *
	 linthurber_vs_dims[2]; i++) {
    vs_buf[i] = -1.0;
  }
  for (i = 0; i < linthurber_dem_dims[0] * linthurber_dem_dims[1] *
	 linthurber_dem_dims[2]; i++) {
    dem_buf[i] = 0.0;
  }

  /* Load Vp velocity file*/
  sprintf(filename, "%s/lin-thurber.vp", mpath);
  num_read = 0;
  fp = fopen(filename, "r");
  while (!feof(fp)) {
    retval = fscanf(fp, "%*f %*f %f %f %f %f %*f", &dep, &y, &x, &val);
    if (retval != EOF) {
      if (retval != 4) {
	fprintf(stderr, 
		"Failed to read Lin-Thurber Vp file, line %d (parsed %d)\n",
		num_read, retval);
	return(UCVM_CODE_ERROR);
      }
      for (k = 0; k < linthurber_z_dim; k++) {
	if (linthurber_depths_msl[k] >= dep) {
	  break;
	}
      }
      /* Flip x and y axis */
      j = round((y + (-linthurber_vp_origin[0])) 
		* 1000.0 / linthurber_vp_spacing);
      i = round((linthurber_proj.dims[0]/1000.0 - 
		 (x + (-linthurber_vp_origin[1]))) 
		* 1000.0 / linthurber_vp_spacing);
      if ((i < 0) || (j < 0) || (k < 0) || 
	  (i >= linthurber_vp_dims[0]) || 
	  (j >= linthurber_vp_dims[1]) || 
	  (k >= linthurber_vp_dims[2])) {
	fprintf(stderr, "Invalid index %d,%d,%d calculated\n", i, j, k);
	return(UCVM_CODE_ERROR);
      }
      //printf("x,y,z: %d, %d, %d\n", i, j, k);
      vp_buf[k*linthurber_vp_dims[0]*linthurber_vp_dims[1] + 
	     j*linthurber_vp_dims[0] + i] = val * 1000.0;
      num_read++;
    }
  }
  fclose(fp);

  /* Load Vs velocity file */
  sprintf(filename, "%s/lin-thurber.vs", mpath);
  num_read = 0;
  fp = fopen(filename, "r");
  while (!feof(fp)) {
    retval = fscanf(fp, "%*f %*f %f %f %f %f %*f", &dep, &y, &x, &val);
    if (retval != EOF) {
      if (retval != 4) {
	fprintf(stderr, 
		"Failed to read Lin-Thurber Vs file, line %d (parsed %d)\n",
		num_read, retval);
	return(UCVM_CODE_ERROR);
      }
      for (k = 0; k < linthurber_z_dim; k++) {
	if (linthurber_depths_msl[k] >= dep) {
	  break;
	}
      }
      /* Flip x and y axis */
      j = round((y + (-linthurber_vs_origin[0])) 
		* 1000.0 / linthurber_vs_spacing);
      i = round((linthurber_proj.dims[0]/1000.0 - 
		 (x + (-linthurber_vs_origin[1]))) 
		* 1000.0 / linthurber_vs_spacing);
      if ((i < 0) || (j < 0) || (k < 0) || 
	  (i >= linthurber_vs_dims[0]) || 
	  (j >= linthurber_vs_dims[1]) || 
	  (k >= linthurber_vs_dims[2])) {
	fprintf(stderr, "Invalid index %d,%d,%d calculated\n", i, j, k);
	return(UCVM_CODE_ERROR);
      }
      //printf("x,y,z: %d, %d, %d\n", i, j, k);
      vs_buf[k*linthurber_vs_dims[0]*linthurber_vs_dims[1] + 
	     j*linthurber_vs_dims[0] + i] = val * 1000.0;
      num_read++;
    }
  }
  fclose(fp);

  /* Load DEM file */
  sprintf(filename, "%s/lin-thurber.dem", mpath);
  fp = fopen(filename, "r");
  num_read = fread(dem_buf, sizeof(float), 
		   linthurber_dem_dims[0] * linthurber_dem_dims[1] *
		   linthurber_dem_dims[2], fp);
  if (num_read != linthurber_dem_dims[0] * linthurber_dem_dims[1] *
      linthurber_dem_dims[2]) {
    fprintf(stderr, "Failed to read Lin-Thurber DEM file\n");
    return(UCVM_CODE_ERROR);
  }

  /* Swap endian from LSB to MSB */
  if (system_endian() == UCVM_BYTEORDER_MSB) {
    for (i = 0; i < linthurber_dem_dims[0] * linthurber_dem_dims[1] *
	   linthurber_dem_dims[2]; i++) {
      dem_buf[i] = swap_endian_float(dem_buf[i]);
    }
  }
  fclose(fp);

  return(UCVM_CODE_SUCCESS);
}


/* Finalize LinThurber */
int ucvm_linthurber_model_finalize()
{
  if (vp_buf != NULL) {
    free(vp_buf);
  }
  if (vs_buf != NULL) {
    free(vs_buf);
  }
  if (dem_buf != NULL) {
    free(dem_buf);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Version LinThurber */
const char *ucvm_linthurber_model_version()
{
  return(ucvm_linthurber_version_id);
}


int ucvm_linthurber_getval(double i, double j, double k, int prop,
			   double *val)
{
  int i0, j0, k0;
  int a, b, c, x, y, z;
  int *dims = NULL;
  float *buf = NULL;
  double p[2][3];
  double q[2][2][2];

  *val = -1.0;

  i0 = (int)i;
  j0 = (int)j;
  k0 = (int)k;

  a = round(i);
  b = round(j);
  c = round(k);

  //fprintf(stderr, "ijk = %lf, %lf, %lf\n", i, j, k);

  switch (prop) {
  case LINTHURBER_DEM:
    dims = linthurber_dem_dims;
    buf = dem_buf;
    break;
  case LINTHURBER_VP:
    dims = linthurber_vp_dims;
    buf = vp_buf;
    break;
  case LINTHURBER_VS:
    dims = linthurber_vs_dims;
    buf = vs_buf;
   break;
  };

  /* Check if point falls outside of model region */
  if ((a < 0) || (b < 0) || (c < 0) || 
      (a >= dims[0]) || (b >= dims[1]) || (c >= dims[2])) {
    //fprintf(stderr, "a,b,c = %d, %d, %d\n", a, b, c);
    return(UCVM_CODE_ERROR);
  }

  /* Values at corners of interpolation cube */
  for (z = 0; z < 2; z++) {
    for (y = 0; y < 2; y++) {
      for (x = 0; x < 2; x++) {

	a = i0 + x;
	b = j0 + y;
	c = k0 + z;

	if (a < 0) {
	  a = 0;
	}
	if (b < 0) {
	  b = 0;
	}
	if (c < 0) {
	  c = 0;
	}
	if (a >= dims[0]) {
	  a = dims[0] - 1;
	}
	if (b >= dims[1]) {
	  b = dims[1] - 1;
	}
	if (c >= dims[2]) {
	  c = dims[2] - 1;
	}

	q[z][y][x] = buf[c*dims[0]*dims[1] + b*dims[0] + a];
      }
    }
  }

  /* Corners of interpolation cube */
  for (b = 0; b < 2; b++) {
    for (a = 0; a < 3; a++) { 
      p[b][a] = (double)b;
    }
  }

  /* Trilinear interpolation */
  *val = interpolate_trilinear_1d(i-i0, j-j0, k-k0, p, q);

  return(UCVM_CODE_SUCCESS);
}


/* Density derived from Vp via Nafe-Drake curve, Brocher (2005) eqn 1. */
double ucvm_linthurber_rho(double f) {
  double rho;

  /* Convert m to km */
  f = f / 1000.0;
  rho = f * (1.6612 - f * (0.4721 - f * (0.0671 - f * (0.0043 - f * 0.000106))));
  if (rho < 1.0) {
    rho = 1.0;
  }
  rho = rho * 1000.0;
  return(rho);
}


/* Query LinThurber */
int ucvm_linthurber_model_query(ucvm_ctype_t cmode,
				int n, ucvm_point_t *pnt, 
				ucvm_prop_t *prop)
{
  int k, p;
  double x, y, z, depth_msl, depth_ratio;
  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Only geo depth mode supported\n");
    return(UCVM_CODE_ERROR);
  }

  /* Query model */
  for (p = 0; p < n; p++) {
    if (prop[p].model == UCVM_MODEL_NONE) {
      if (bilinear_geo2xy(&linthurber_proj,
			  pnt[p].coord[0], pnt[p].coord[1], &x, &y) == 0) {
	prop[p].vp = -1.0;
	prop[p].vs = -1.0;
	prop[p].rho = -1.0;

	if (LINTHURBER_USE_DEM) {

	  /* DEM elevation value */
	  if (ucvm_linthurber_getval(x/linthurber_dem_spacing,
				     y/linthurber_dem_spacing,
				     0.0, LINTHURBER_DEM, 
				     &depth_msl) != UCVM_CODE_SUCCESS) {
	    continue;
	  }
	  depth_msl = (pnt[p].coord[2] - depth_msl)/1000.0;
	  for (k = 0; k < linthurber_z_dim; k++) {
	    if (linthurber_depths_msl[k] >= depth_msl) {
	      break;
	    }
	  }

	  if (k == linthurber_z_dim) {
	    k = linthurber_z_dim - 1;
	    depth_ratio = 0.0;
	    z = k;
	  } else if (k == 0) {
	    depth_ratio = 0.0;
	    z = k;
	  } else {
	    depth_ratio = (depth_msl - linthurber_depths_msl[k-1]) /
	      (linthurber_depths_msl[k] - linthurber_depths_msl[k-1]); 
	    z = (k-1) + depth_ratio;
	  }
	  
	  /* Vp value */
	  ucvm_linthurber_getval(x / linthurber_vp_spacing,
				 y / linthurber_vp_spacing,
				 z, LINTHURBER_VP, &(prop[p].vp));
	  
	  /* Vs value */
	  ucvm_linthurber_getval(x / linthurber_vs_spacing,
				 y / linthurber_vs_spacing,
				 z, LINTHURBER_VS, &(prop[p].vs));
	} else {
	  depth_msl = pnt[p].coord[2]/1000.0;
	  for (k = 0; k < linthurber_z_dim; k++) {
	    if (linthurber_depths_msl[k] >= depth_msl) {
	      break;
	    }
	  }
	  if (k == linthurber_z_dim) {
	    k = linthurber_z_dim - 1;
	  }
	  z = k;

	  /* Vp value */
	  ucvm_linthurber_getval(x / linthurber_vp_spacing,
				 y / linthurber_vp_spacing,
				 z, LINTHURBER_VP, &(prop[p].vp));
	  
	  /* Vs value */
	  ucvm_linthurber_getval(x / linthurber_vs_spacing,
				 y / linthurber_vs_spacing,
				 z, LINTHURBER_VS, &(prop[p].vs));

	}

	/* Calculate density */
	if (prop[p].vp > 0.0) {
	  prop[p].model = ucvm_linthurber_id;
	  prop[p].rho = ucvm_linthurber_rho(prop[p].vp);
	}
      }
    }
  }

  return(UCVM_CODE_SUCCESS);
}


/* Fill model structure with LinThurber */
int ucvm_linthurber_get_model(ucvm_model_t *m)
{
  strcpy(m->label, ucvm_linthurber_label_id);
  memset(&(m->region), 0, sizeof(ucvm_region_t));
  strcpy(m->config, "");
  m->init = ucvm_linthurber_model_init;
  m->finalize = ucvm_linthurber_model_finalize;
  m->version = ucvm_linthurber_model_version;
  m->query = ucvm_linthurber_model_query;

  return(UCVM_CODE_SUCCESS);
}


