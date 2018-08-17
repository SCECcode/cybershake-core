#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ucvm_utils.h"
#include "ucvm_model_cvms.h"
#include "cvms.h"

/* Coverage region */
ucvm_region_t cvms_region;

/* Maximum number of points to query */
#define CVMS_MAX_POINTS 1000000

/* Buffers initialized flag */
int cvms_buf_init = 0;

/* Query buffers */
int *cvms_index = NULL;
float *cvms_lon = NULL;
float *cvms_lat = NULL;
float *cvms_dep = NULL;
float *cvms_vp = NULL;
float *cvms_vs = NULL;
float *cvms_rho = NULL;


/* Version string */
const char *ucvm_cvms_version_id = "CVM-S4 + Ely/Graves/Small Mods";

/* Label string */
const char *ucvm_cvms_label_id = "cvms";

/* Model ID */
int ucvm_cvms_id = UCVM_MODEL_NONE;


/* Init CVM-S */
int ucvm_cvms_model_init(int id, ucvm_region_t *r, const char *mpath)
{
  int errcode;
  char config[UCVM_MAX_PATH_LEN];

  if (strlen(mpath) > UCVM_MAX_PATH_LEN) {
    fprintf(stderr, "Model path is too long\n");
    return(UCVM_CODE_ERROR);
  }
  strcpy(config, mpath);
  cvms_init_(config, &errcode);
  if (errcode != 0) {
    fprintf(stderr, "Failed to init CVM-S\n");
    return(UCVM_CODE_ERROR);
  }

  /* Allocate buffers */
  if (cvms_buf_init == 0) {
    cvms_index = malloc(CVMS_MAX_POINTS*sizeof(int));
    cvms_lon = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_lat = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_dep = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_vp = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_vs = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_rho = malloc(CVMS_MAX_POINTS*sizeof(float));
    cvms_buf_init = 1;
  }

  ucvm_cvms_id = id;

  /* Save coverage region */
  memcpy(&cvms_region, r, sizeof(ucvm_region_t));

  return(UCVM_CODE_SUCCESS);
}


/* Finalize CVM-S */
int ucvm_cvms_model_finalize()
{
  if (cvms_buf_init == 1) {
    free(cvms_index);
    free(cvms_lon);
    free(cvms_lat);
    free(cvms_dep);
    free(cvms_vp);
    free(cvms_vs);
    free(cvms_rho);
  }
  return(UCVM_CODE_SUCCESS);
}


/* Version CVM-S */
const char *ucvm_cvms_model_version()
{
  return(ucvm_cvms_version_id);
}


/* Query CVM-S */
int ucvm_cvms_model_query(ucvm_ctype_t cmode,
			  int n, ucvm_point_t *pnt, ucvm_prop_t *prop)
{
  int i, j;
  int nn = 0;
  int datagap = 0;

  if (cvms_buf_init == 0) {
    return(UCVM_CODE_ERROR);
  }

  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Only geo depth mode supported\n");
    return(UCVM_CODE_ERROR);
  }

  nn = 0;
  for (i = 0; i < n; i++) {

    //if (prop[i].model == UCVM_MODEL_NONE) {

    if ((prop[i].model == UCVM_MODEL_NONE) && 
      	(region_contains(&cvms_region, cmode, &(pnt[i])))) {

      /* Query point */
      cvms_index[nn] = i;
      cvms_lon[nn] = (float)(pnt[i].coord[0]);
      cvms_lat[nn] = (float)(pnt[i].coord[1]);
      cvms_dep[nn] = (float)(pnt[i].coord[2]);
      cvms_vp[nn] = 0.0;
      cvms_vs[nn] = 0.0;
      cvms_rho[nn] = 0.0;
      nn++;
      if (nn == CVMS_MAX_POINTS) {
	cvms_query_(&nn, cvms_lon, cvms_lat, cvms_dep, 
		    cvms_vp, cvms_vs, cvms_rho);
	
	for (j = 0; j < nn; j++) {
	  prop[cvms_index[j]].vp = (double)cvms_vp[j];
	  prop[cvms_index[j]].vs = (double)cvms_vs[j];
	  prop[cvms_index[j]].rho = (double)cvms_rho[j];
	  prop[cvms_index[j]].model = ucvm_cvms_id;
	}

	nn = 0;
      }
    } else {
      if (prop[i].model == UCVM_MODEL_NONE) {
	datagap = 1;
      }
    }
  }

  if (nn > 0) {
    cvms_query_(&nn, cvms_lon, cvms_lat, cvms_dep, 
		cvms_vp, cvms_vs, cvms_rho);

    for (j = 0; j < nn; j++) {
      prop[cvms_index[j]].vp = (double)cvms_vp[j];
      prop[cvms_index[j]].vs = (double)cvms_vs[j];
      prop[cvms_index[j]].rho = (double)cvms_rho[j];
      prop[cvms_index[j]].model = ucvm_cvms_id;
    }

  }

  if (datagap) {
    return(UCVM_CODE_DATAGAP);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Fill model structure with CVM-S */
int ucvm_cvms_get_model(ucvm_model_t *m)
{
  strcpy(m->label, ucvm_cvms_label_id);
  memset(&(m->region), 0, sizeof(ucvm_region_t));
  strcpy(m->config, "");
  m->init = ucvm_cvms_model_init;
  m->finalize = ucvm_cvms_model_finalize;
  m->version = ucvm_cvms_model_version;
  m->query = ucvm_cvms_model_query;

  return(UCVM_CODE_SUCCESS);
}
