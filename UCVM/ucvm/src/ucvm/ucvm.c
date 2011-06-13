#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ucvm.h"
#include "ucvm_config.h"
#include "ucvm_utils.h"
#ifdef _UCVM_ENABLE_CVMS
#include "ucvm_model_cvms.h"
#endif
#ifdef _UCVM_ENABLE_CVMH
#include "ucvm_model_cvmh.h"
#endif
#ifdef _UCVM_ENABLE_CENCAL
#include "ucvm_model_cencal.h"
#endif
#include "ucvm_model_1d.h"
#include "ucvm_model_linthurber.h"


/* Init flag */
int ucvm_init_flag = 0;


/* Current query mode */
ucvm_ctype_t ucvm_cur_qmode = UCVM_COORD_GEO_DEPTH;


/* Model list */
int ucvm_num_models = 0;
ucvm_model_t ucvm_model_list[UCVM_MAX_MODELS];


/* UCVM config */
ucvm_config_t *ucvm_cfg = NULL;


/* Initializer */
int ucvm_init(const char *config)
{
  ucvm_init_flag = 0;

  ucvm_num_models = 0;
  memset(ucvm_model_list, 0, sizeof(ucvm_model_t)*UCVM_MAX_MODELS);

  /* Read in config file */
  ucvm_cfg = ucvm_parse_config(config);
  if (ucvm_cfg == NULL) {
    fprintf(stderr, "Failed to read config file\n");
    return(UCVM_CODE_ERROR);
  }

  ucvm_init_flag = 1;
  return(UCVM_CODE_SUCCESS);
}


/* Finalizer */
int ucvm_finalize()
{
  int i;

  /* Call all model finalizers */
  for (i = 0; i < ucvm_num_models; i++) {
    (ucvm_model_list[i].finalize)();
  }

  /* Free config file parser resources */
  if (ucvm_cfg != NULL) {
    ucvm_free_config(ucvm_cfg);
  }

  ucvm_num_models = 0;
  ucvm_init_flag = 0;
  return(UCVM_CODE_SUCCESS);
}


/* Set query mode */
int ucvm_query_mode(ucvm_ctype_t c)
{
  if (ucvm_init_flag == 0) {
    fprintf(stderr, "UCVM not initialized\n");
    return(UCVM_CODE_ERROR);
  }

  ucvm_cur_qmode = c;
  return(UCVM_CODE_SUCCESS);
}


/* Enable specific model, by label */
int ucvm_add_model(const char *label) {
  ucvm_model_t m;
  int retval = UCVM_CODE_ERROR;

#ifdef _UCVM_ENABLE_CVMS
  if (strcmp(label, UCVM_MODEL_CVMS) == 0) {
    retval = ucvm_cvms_get_model(&m);
  }
#endif

#ifdef _UCVM_ENABLE_CVMH
  if (strcmp(label, UCVM_MODEL_CVMH) == 0) {
    retval = ucvm_cvmh_get_model(&m);
  }
#endif

#ifdef _UCVM_ENABLE_CENCAL
  if (strcmp(label, UCVM_MODEL_CENCAL) == 0) {
    retval = ucvm_cencal_get_model(&m);
  }
#endif

  if (strcmp(label, UCVM_MODEL_1D) == 0) {
    retval = ucvm_1d_get_model(&m);
  }
  if (strcmp(label, UCVM_MODEL_LINTHURBER) == 0) {
    retval = ucvm_linthurber_get_model(&m);
  }

  if (retval != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to find model %s\n", label);
    fprintf(stderr, "Model may not be enabled in UCVM\n");
    return(UCVM_CODE_ERROR);
  }

  /* Register the model */
  return(ucvm_add_user_model(&m));
}


/* Enable specific model, by ucvm_model_t */
int ucvm_add_user_model(ucvm_model_t *m)
{
  int i;
  ucvm_model_t *mptr;
  char key[UCVM_CONFIG_MAX_STR];
  ucvm_config_t *cfgentry = NULL;

  if (ucvm_init_flag == 0) {
    fprintf(stderr, "UCVM not initialized\n");
    return(UCVM_CODE_ERROR);
  }

  if (ucvm_num_models >= UCVM_MAX_MODELS) {
    fprintf(stderr, "Maximum number of models reached\n");
    return(UCVM_CODE_ERROR);
  }

  /* Check if model already enabled */
  for (i = 0; i < ucvm_num_models; i++) {
    mptr = &(ucvm_model_list[i]);
    if (strcmp(mptr->label, m->label) == 0) {
      fprintf(stderr, "Model %s already enabled\n", mptr->label);
      return(UCVM_CODE_ERROR);
    }
  }

  mptr = &(ucvm_model_list[ucvm_num_models]);
  
  /* Enable this model */
  memcpy(mptr, m, sizeof(ucvm_model_t));

  /* Over-ride default config with config file values */
  sprintf(key, "%s_region", mptr->label);
  cfgentry = ucvm_find_name(ucvm_cfg, key);
  if (cfgentry != NULL) {
    region_parse(cfgentry->value, &(mptr->region));
  }

  sprintf(key, "%s_modelpath", mptr->label);
  cfgentry = ucvm_find_name(ucvm_cfg, key);
  if (cfgentry != NULL) {
    strcpy(mptr->config, cfgentry->value);
    sprintf(key, "%s_extmodelpath", mptr->label);
    cfgentry = ucvm_find_name(ucvm_cfg, key);
    if (cfgentry != NULL) {
      strcat(mptr->config, ";");
      strcat(mptr->config, cfgentry->value);
    }
  }

  /* Perform init */
  if ((mptr->init)(ucvm_num_models, 
		   &(mptr->region), 
		   mptr->config) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Failed to init model %s with cfg %s\n", 
	    mptr->label, mptr->config);
    return(UCVM_CODE_ERROR);
  }

  ucvm_num_models++;

  return(UCVM_CODE_SUCCESS);
}


/* Get label for a model */
int ucvm_model_label(int m, char *label, int len)
{
  int labellen;
  ucvm_model_t *mptr;

  if (ucvm_init_flag == 0) {
    fprintf(stderr, "UCVM not initialized\n");
    return(UCVM_CODE_ERROR);
  }

  if (m == UCVM_MODEL_NONE) {
    labellen = strlen("none");
    if (labellen > len) {
      labellen = len;
    }
    strncpy(label, "none", labellen);
    return(UCVM_CODE_SUCCESS);
  }

  if (m >= ucvm_num_models) {
    fprintf(stderr, "Invalid model ID %d\n", m);
    return(UCVM_CODE_ERROR);
  }

  mptr = &(ucvm_model_list[m]);
  labellen = strlen(mptr->label);
  if (labellen > len) {
    labellen = len;
  }
  strncpy(label, mptr->label, len);
 
  return(UCVM_CODE_SUCCESS);
}


/* Get version for a model */
int ucvm_model_version(int m, char *ver, int len)
{
  ucvm_model_t *mptr;
  const char *vstr = NULL;
  int verlen;

  if (ucvm_init_flag == 0) {
    fprintf(stderr, "UCVM not initialized\n");
    return(UCVM_CODE_ERROR);
  }

  if (m == UCVM_MODEL_NONE) {
    verlen = strlen("none");
    if (verlen > len) {
      verlen = len;
    }
    strncpy(ver, "none", verlen);
    return(UCVM_CODE_SUCCESS);
  }

  if (m >= ucvm_num_models) {
    fprintf(stderr, "Invalid model ID %d\n", m);
    return(UCVM_CODE_ERROR);
  }

  mptr = &(ucvm_model_list[m]);
  vstr = (mptr->version)();
  verlen = strlen(vstr);
  if (verlen > len) {
    verlen = len;
  }

  strncpy(ver, vstr, verlen);
  return(UCVM_CODE_SUCCESS);
}


/* Query underlying models */
int ucvm_query(int n, ucvm_point_t *pnt, ucvm_prop_t *prop)
{
  int i;
  ucvm_model_t *mptr;

  if (ucvm_init_flag == 0) {
    fprintf(stderr, "UCVM not initialized\n");
    return(UCVM_CODE_ERROR);
  }

  if (ucvm_num_models == 0) {
    fprintf(stderr, "No models enabled\n");
    return(UCVM_CODE_ERROR);
  }

  /* Initialize properties array */
  for (i = 0; i < n; i++) {
    prop[i].model = UCVM_MODEL_NONE;
    prop[i].vp = 0.0;
    prop[i].vs = 0.0;
    prop[i].rho = 0.0;
  }

  for (i = 0; i < ucvm_num_models; i++) {
    mptr = &(ucvm_model_list[i]);
    if ((mptr->query)(ucvm_cur_qmode, n, pnt, prop) == UCVM_CODE_SUCCESS) {
      return(UCVM_CODE_SUCCESS);
    }
  }

  return(UCVM_CODE_SUCCESS);
}
