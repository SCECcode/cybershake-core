#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ucvm_model_1d.h"


/* SCEC 1D depth (km) -> vs (km/s) array */
#define UCVM_MAX_SCEC_1D 9
double ucvm_layer_depths[UCVM_MAX_SCEC_1D] = 
  {1.0, 5.0, 6.0, 10.0, 15.5, 16.5, 22.0, 31.0, 33.0};
double ucvm_layer_vp[UCVM_MAX_SCEC_1D] = 
  {5.0, 5.5, 6.3, 6.3, 6.4, 6.7, 6.75, 6.8, 7.8};

/* Version string */
const char *ucvm_1d_version_id = "Hadley-Kanamori 1D";

/* Label string */
const char *ucvm_1d_label_id = "1d";

/* Model ID */
int ucvm_1d_id = UCVM_MODEL_NONE;



/* Determine vp by depth */
double ucvm_scec_vp(double depth) {
  int i;
  double vp;
  double depth_ratio;
  double vp_range;

  /* Convert from m to km */
  depth = depth / 1000.0;

  /* Scale vp by depth with linear interpolation */
  for (i = 0; i < UCVM_MAX_SCEC_1D; i++) {
    if (ucvm_layer_depths[i] > depth) {
      if (i == 0) {
        vp = ucvm_layer_vp[i];
      } else {
        depth_ratio = ((depth - ucvm_layer_depths[i-1]) / 
                       (ucvm_layer_depths[i] - ucvm_layer_depths[i - 1]));
        vp_range = ucvm_layer_vp[i] - ucvm_layer_vp[i - 1];
        vp = ((vp_range * depth_ratio) + ucvm_layer_vp[i - 1]);
      }
      break;
    } 
  }
  if (i == UCVM_MAX_SCEC_1D) {
    vp = ucvm_layer_vp[UCVM_MAX_SCEC_1D - 1];
  }

  /* Convert from km/s back to m/s */
  vp = vp * 1000.0;

  return(vp);
}


/* Determine density by vp */
double ucvm_scec_rho(double vp) {
  double rho;

  /* Calculate rho */
  rho = 1865.0 + 0.1579 * vp;
  return(rho);
}


/* Determine vs by vp and density */
double ucvm_scec_vs(double vp, double rho) {
  double vs;
  double nu;

  if (rho < 2060.0) {
    nu = 0.40;
  } else if (rho > 2500.0) {
    nu = 0.25;
  } else {
    nu = 0.40 - ((rho - 2060.0) * 0.15 / 440.0);
  }

  vs = vp * sqrt( (0.5 - nu) / (1.0 - nu) );

  return(vs);
}


/* Init 1D */
int ucvm_1d_model_init(int m, ucvm_region_t *r, const char *mpath)
{
  ucvm_1d_id = m;
  return(UCVM_CODE_SUCCESS);
}


/* Finalize 1D */
int ucvm_1d_model_finalize()
{
  return(UCVM_CODE_SUCCESS);
}


/* Version 1D */
const char *ucvm_1d_model_version()
{
  return(ucvm_1d_version_id);
}


/* Query 1D */
int ucvm_1d_model_query(ucvm_ctype_t cmode,
			int n, ucvm_point_t *pnt, ucvm_prop_t *prop)
{
  int i;

  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Only geo depth mode supported\n");
    return(UCVM_CODE_ERROR);
  }

  for (i = 0; i < n; i++) {
    if (prop[i].model == UCVM_MODEL_NONE) {
      prop[i].vp = ucvm_scec_vp(pnt[i].coord[2]);
      prop[i].rho = ucvm_scec_rho(prop[i].vp);
      prop[i].vs = ucvm_scec_vs(prop[i].vp, prop[i].rho);
      prop[i].model = ucvm_1d_id;
    }
  }

  return(UCVM_CODE_SUCCESS);
}


/* Fill model structure with 1D */
int ucvm_1d_get_model(ucvm_model_t *m)
{
  strcpy(m->label, ucvm_1d_label_id);
  memset(&(m->region), 0, sizeof(ucvm_region_t));
  strcpy(m->config, "");
  m->init = ucvm_1d_model_init;
  m->finalize = ucvm_1d_model_finalize;
  m->version = ucvm_1d_model_version;
  m->query = ucvm_1d_model_query;

  return(UCVM_CODE_SUCCESS);
}

