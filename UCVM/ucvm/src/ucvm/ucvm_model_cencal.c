#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ucvm_utils.h"
#include "ucvm_model_cencal.h"
#include "cvmquery.h"
#include "cvmerror.h"

/* Config delimiter */
#define MPATH_DELIM ";"

/* Buffers initialized flag */
int cc_buf_init = 0;

/* Coverage region */
ucvm_region_t cencal_region;

/* Query buffers */
void *query = NULL;
void *errHandler = NULL;

/* Default cache size */
#define CC_CACHE_SIZE 128

/* Query by depth squash limit */
#define CC_SQUASH_LIMIT -100000.0

/* CenCal number of values */
#define CC_NUM_VALS 9

/* CenCal return value buffer */
double *cc_pvals = NULL;

/* Version string */
const char *ucvm_cencal_version_id = "cencalvm 0.6.6";

/* Label string */
const char *ucvm_cencal_label_id = "cencal";

/* Model ID */
int ucvm_cc_id = UCVM_MODEL_NONE;


/* Init CenCal */
int ucvm_cencal_model_init(int id, ucvm_region_t *r, const char *mpath)
{
  char buf[UCVM_MAX_PATH_LEN];
  char *mp = NULL;
  char *mp_ext = NULL;

  strcpy(buf, mpath);
  mp = strtok(buf, MPATH_DELIM);
  if (mp == NULL) {
    mp = buf;
  } else {
    mp_ext = strtok(NULL, MPATH_DELIM);
  }

  /* Create query */
  query = cencalvm_createQuery();
  if (query == NULL) {
    fprintf(stderr, "Could not create query.\n");
    return(UCVM_CODE_ERROR);
  }

  /* Get handle to error handler */
  errHandler = cencalvm_errorHandler(query);
  if (errHandler == NULL) {
    fprintf(stderr, "Could not get handle to error handler.\n");
    return(UCVM_CODE_ERROR);
  }

  /* Set database filename */
  if (cencalvm_filename(query, mp) != 0) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return(UCVM_CODE_ERROR);
  }

  /* Set extended DB here */
  if (mp_ext != NULL) {
    if (cencalvm_filenameExt(query, mp_ext) != 0) {
      fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
      return(UCVM_CODE_ERROR);
    }
  }

  /* Set cache size */
  if (cencalvm_cacheSize(query, CC_CACHE_SIZE) != 0) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return(UCVM_CODE_ERROR);
  }

  /* Turn on squashing if requested */
  if (cencalvm_squash(query, 1, CC_SQUASH_LIMIT) != 0) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return(UCVM_CODE_ERROR);
  }

  /* Open database for querying */
  if (cencalvm_open(query) != 0) {
    fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
    return(UCVM_CODE_ERROR);
  }

 /* Set query type and resolution */
  cencalvm_queryType(query, 0);

  /* Create array to hold values returned in queries */
  cc_pvals = (double*) malloc(sizeof(double)*CC_NUM_VALS);
  if (cc_pvals == NULL) {
    return(UCVM_CODE_ERROR);
  }

  ucvm_cc_id = id;

  /* Save coverage region */
  memcpy(&cencal_region, r, sizeof(ucvm_region_t));

  cc_buf_init = 1;
  return(UCVM_CODE_SUCCESS);
}


/* Finalize CenCal */
int ucvm_cencal_model_finalize()
{

  /* Close database */
  if (query != NULL) {
    cencalvm_close(query);
  }

  /* Destroy query handle */
  cencalvm_destroyQuery(query);

  /* Free data buffer */
  if (cc_pvals != NULL) {
    free(cc_pvals);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Version Cencal */
const char *ucvm_cencal_model_version()
{
  return(ucvm_cencal_version_id);
}


/* Query CenCal */
int ucvm_cencal_model_query(ucvm_ctype_t cmode,
			    int n, ucvm_point_t *pnt, ucvm_prop_t *prop)
{
  int i;
  double lon, lat, elev;
  int datagap = 0;

  if (cc_buf_init == 0) {
    return(UCVM_CODE_ERROR);
  }

  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Only geo depth mode supported\n");
    return(UCVM_CODE_ERROR);
  }

  for (i = 0; i < n; i++) {
    if ((prop[i].model == UCVM_MODEL_NONE) && 
	(region_contains(&cencal_region, cmode, &(pnt[i])))) {

      /* Query point */
      lon = pnt[i].coord[0];
      lat = pnt[i].coord[1];
      elev = -pnt[i].coord[2];

      if (cencalvm_query(query, &cc_pvals, CC_NUM_VALS, lon, lat, elev) != 0) {
	/* NOTE comment datagap to force no data return */
	datagap = 1;
	//fprintf(stderr, "%s\n", cencalvm_error_message(errHandler));
	cencalvm_error_resetStatus(errHandler);
      } else {
	prop[i].vp = cc_pvals[0];
	prop[i].vs = cc_pvals[1];
	prop[i].rho = cc_pvals[2];
	prop[i].model = ucvm_cc_id;
      }
    } else {
      if (prop[i].model == UCVM_MODEL_NONE) {
	datagap = 1;
      }
    }
  }

  if (datagap) {
    return(UCVM_CODE_DATAGAP);
  }

  return(UCVM_CODE_SUCCESS);
}


/* Fill model structure with CenCal */
int ucvm_cencal_get_model(ucvm_model_t *m)
{
  strcpy(m->label, ucvm_cencal_label_id);
  memset(&(m->region), 0, sizeof(ucvm_region_t));
  strcpy(m->config, "");
  m->init = ucvm_cencal_model_init;
  m->finalize = ucvm_cencal_model_finalize;
  m->version = ucvm_cencal_model_version;
  m->query = ucvm_cencal_model_query;

  return(UCVM_CODE_SUCCESS);
}
