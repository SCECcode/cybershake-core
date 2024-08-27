#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ucvm_model_1d.h"
#include "ucvm_utils.h"
#include "ucvm_config.h"

/* Init flag */
int ucvm_1d_init_flag = 0;

/* Model dimensions */
#define UCVM_1D_MAX_Z_DIM 32
int ucvm_1d_z_dim = 0;
char ucvm_1d_version_id[UCVM_MAX_VERSION_LEN];
double ucvm_1d_layer_depths[UCVM_1D_MAX_Z_DIM];
double ucvm_1d_layer_vp[UCVM_1D_MAX_Z_DIM];

/* Model ID */
int ucvm_1d_id = UCVM_SOURCE_NONE;

/* Model conf */
ucvm_modelconf_t ucvm_1d_conf;

/* Model config */
ucvm_config_t *ucvm_1d_cfg = NULL;

/* Determine vp by depth */
double ucvm_1d_scec_vp(double depth) {
  int i;
  double vp;
  double depth_ratio;
  double vp_range;

  vp = 0.0;

  /* Convert from m to km */
  depth = depth / 1000.0;

  /* Scale vp by depth with linear interpolation */
  for (i = 0; i < ucvm_1d_z_dim; i++) {
    if (ucvm_1d_layer_depths[i] > depth) {
      if (i == 0) {
        vp = ucvm_1d_layer_vp[i];
      } else {
        depth_ratio = ((depth - ucvm_1d_layer_depths[i-1]) / 
                       (ucvm_1d_layer_depths[i] - ucvm_1d_layer_depths[i - 1]));
        vp_range = ucvm_1d_layer_vp[i] - ucvm_1d_layer_vp[i - 1];
        vp = ((vp_range * depth_ratio) + ucvm_1d_layer_vp[i - 1]);
      }
      break;
    } 
  }
  if (i == ucvm_1d_z_dim) {
    vp = ucvm_1d_layer_vp[ucvm_1d_z_dim - 1];
  }

  /* Convert from km/s back to m/s */
  vp = vp * 1000.0;

  return(vp);
}


/* Determine density by vp */
double ucvm_1d_scec_rho(double vp) {
  double rho;

  /* Calculate rho */
  rho = 1865.0 + 0.1579 * vp;
  return(rho);
}


/* Determine vs by vp and density */
double ucvm_1d_scec_vs(double vp, double rho) {
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
int ucvm_1d_model_init(int m, ucvm_modelconf_t *conf)
{
  int i, len;
  ucvm_config_t *chead;
  ucvm_config_t *cptr;

  if (ucvm_1d_init_flag) {
    fprintf(stderr, "Model %s is already initialized\n", conf->label);
    return(UCVM_CODE_ERROR);
  }

  if ((conf->config == NULL) || (strlen(conf->config) == 0)) {
    fprintf(stderr, "No config path defined for model %s\n", conf->label);
    return(UCVM_CODE_ERROR);
  }

  if (!ucvm_is_file(conf->config)) {
    fprintf(stderr, "1D conf file %s is not a valid file\n", 
            conf->config);
    return(UCVM_CODE_ERROR);
  }

  /* Save model conf */
  memcpy(&ucvm_1d_conf, conf, sizeof(ucvm_modelconf_t));

  /* Read conf file */
  chead = ucvm_parse_config(conf->config);
  if (chead == NULL) {
    fprintf(stderr, "Failed to parse config file %s\n", conf->config);
    return(UCVM_CODE_ERROR);
  }

  cptr = ucvm_find_name(chead, "version");
  if (cptr == NULL) {
    fprintf(stderr, "Failed to find version in config\n");
    return(UCVM_CODE_ERROR);
  }
  len = strlen(cptr->value);
  if(len <= UCVM_MAX_VERSION_LEN) {
    len = UCVM_MAX_VERSION_LEN -1;
  }
  snprintf(ucvm_1d_version_id, len, "%s", cptr->value);

  cptr = ucvm_find_name(chead, "num_z");
  if (cptr == NULL) {
    fprintf(stderr, "Failed to find num_z key\n");
    return(1);
  }
  ucvm_1d_z_dim = atoi(cptr->value);
  if ((ucvm_1d_z_dim <= 0) || (ucvm_1d_z_dim > UCVM_1D_MAX_Z_DIM)) {
    fprintf(stderr, "Invalid 1D Z dimension size\n");
    return(UCVM_CODE_ERROR);
  }

  cptr = ucvm_find_name(chead, "vp_vals");
  if (cptr == NULL) {
    fprintf(stderr, "Failed to find vp_vals key\n");
    return(1);
  }
  if (list_parse(cptr->value, UCVM_CONFIG_MAX_STR, 
		 ucvm_1d_layer_vp, ucvm_1d_z_dim) != 0) {
    fprintf(stderr, "Failed to parse %s value\n", cptr->name);
    return(1);
  }

  cptr = ucvm_find_name(chead, "depth_vals");
  if (cptr == NULL) {
    fprintf(stderr, "Failed to find depth_vals key\n");
    return(1);
  }
  if (list_parse(cptr->value,  UCVM_CONFIG_MAX_STR, 
		 ucvm_1d_layer_depths, ucvm_1d_z_dim) != 0) {
    fprintf(stderr, "Failed to parse %s value\n", cptr->name);
    return(1);
  }

  /* Check for valid layer setup */
  for (i = 1; i < ucvm_1d_z_dim; i++) {
    //printf("%lf, %lf\n", ucvm_1d_layer_depths[i-1], ucvm_1d_layer_vp[i-1]);
    //printf("%lf, %lf\n", ucvm_1d_layer_depths[i], ucvm_1d_layer_vp[i]);
    if ((ucvm_1d_layer_vp[i] < ucvm_1d_layer_vp[i-1]) || 
	(ucvm_1d_layer_depths[i] <= ucvm_1d_layer_depths[i-1])) {
      fprintf(stderr, "1D Depth and Vp must monotonically increase.\n");
      return(UCVM_CODE_ERROR);
    }
    if ((ucvm_1d_layer_vp[i] <= 0.0) || 
	(ucvm_1d_layer_vp[i-1] <= 0.0)) {
      fprintf(stderr, "1D Vp must be positive.\n");
      return(UCVM_CODE_ERROR);
    }
      if ((ucvm_1d_layer_depths[i] < 0.0) || 
	  (ucvm_1d_layer_depths[i-1] < 0.0)) {
      fprintf(stderr, "1D Depth must be >= zero.\n");
      return(UCVM_CODE_ERROR);
    }
  }

  ucvm_1d_id = m;
  ucvm_1d_init_flag = 1;
  ucvm_1d_cfg = chead;

  return(UCVM_CODE_SUCCESS);
}


/* Finalize 1D */
int ucvm_1d_model_finalize()
{
  ucvm_1d_z_dim = 0;
  ucvm_1d_init_flag = 0;

  /* Free config file parser resources */
  if (ucvm_1d_cfg != NULL) {
    ucvm_free_config(ucvm_1d_cfg);
    ucvm_1d_cfg = NULL;
  }

  return(UCVM_CODE_SUCCESS);
}


/* Version 1D */
int ucvm_1d_model_version(int id, char *ver, int len)
{
  if (id != ucvm_1d_id) {
    fprintf(stderr, "Invalid model id\n");
    return(UCVM_CODE_ERROR);
  }

  ucvm_strcpy(ver, ucvm_1d_version_id, len);
  return(UCVM_CODE_SUCCESS);
}


/* Label 1D */
int ucvm_1d_model_label(int id, char *lab, int len)
{
  if (id != ucvm_1d_id) {
    fprintf(stderr, "Invalid model id\n");
    return(UCVM_CODE_ERROR);
  }

  ucvm_strcpy(lab, ucvm_1d_conf.label, len);
  return(UCVM_CODE_SUCCESS);
}


/* Setparam 1D */
int ucvm_1d_model_setparam(int id, int param, ...)
{
  va_list ap;

  if (id != ucvm_1d_id) {
    fprintf(stderr, "Invalid model id\n");
    return(UCVM_CODE_ERROR);
  }

  va_start(ap, param);
  switch (param) {
  default:
    break;
  }

  va_end(ap);

  return(UCVM_CODE_SUCCESS);
}


/* Query 1D */
int ucvm_1d_model_query(int id, ucvm_ctype_t cmode,
			int n, ucvm_point_t *pnt, ucvm_data_t *data)
{
  int i;
  double depth;
  int datagap = 0;

  if (id != ucvm_1d_id) {
    fprintf(stderr, "Invalid model id\n");
    return(UCVM_CODE_ERROR);
  }

  /* Check query mode */
  switch (cmode) {
  case UCVM_COORD_GEO_DEPTH:
  case UCVM_COORD_GEO_ELEV:
    break;
  default:
    fprintf(stderr, "Unsupported coord type\n");
    return(UCVM_CODE_ERROR);
    break;
  }

  for (i = 0; i < n; i++) {
    if ((data[i].crust.source == UCVM_SOURCE_NONE) && 
	((data[i].domain == UCVM_DOMAIN_INTERP) || 
	 (data[i].domain == UCVM_DOMAIN_CRUST))) {

      /* Modify pre-computed depth to account for GTL interp range */
      depth = data[i].depth + data[i].shift_cr;

      /* 1D extends from free surface on down */
      if (depth >= 0.0) {
	data[i].crust.vp = ucvm_1d_scec_vp(depth);
	data[i].crust.rho = ucvm_1d_scec_rho(data[i].crust.vp);
	data[i].crust.vs = ucvm_1d_scec_vs(data[i].crust.vp, 
					   data[i].crust.rho);
	data[i].crust.source = ucvm_1d_id;
      } else {
	datagap = 1;
      }
    } else {
      if (data[i].crust.source == UCVM_SOURCE_NONE) {
        datagap = 1;
      }
    }
  }

  if (datagap) {
    return(UCVM_CODE_DATAGAP);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Fill model structure with 1D */
int ucvm_1d_get_model(ucvm_model_t *m)
{
  m->mtype = UCVM_MODEL_CRUSTAL;
  m->init = ucvm_1d_model_init;
  m->finalize = ucvm_1d_model_finalize;
  m->setparam = ucvm_1d_model_setparam;
  m->getversion = ucvm_1d_model_version;
  m->getconfig = NULL;
  m->getlabel = ucvm_1d_model_label;
  m->query = ucvm_1d_model_query;

  return(UCVM_CODE_SUCCESS);
}

