#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ucvm_utils.h"
#include "ucvm_model_cvmh.h"
#include "vx_sub.h"


/* Coverage region */
ucvm_region_t cvmh_region;

/* Version string */
const char *ucvm_cvmh_version_id = "CVM-H (RC)";

/* Label string */
const char *ucvm_cvmh_label_id = "cvmh";

/* Model ID */
int ucvm_cvmh_id = UCVM_MODEL_NONE;


/* Init CVM-H */
int ucvm_cvmh_model_init(int id, ucvm_region_t *r, const char *mpath)
{
  if (vx_setup(mpath) != 0) {
    return(UCVM_CODE_ERROR);
  }

  /* Register SCEC 1D background model */
  vx_register_scec();

  /* Query by depth */
  vx_setzmode(VX_ZMODE_DEPTH);
  
  ucvm_cvmh_id = id;

  /* Save coverage region */
  memcpy(&cvmh_region, r, sizeof(ucvm_region_t));

  return(UCVM_CODE_SUCCESS);
}


/* Finalize CVM-H */
int ucvm_cvmh_model_finalize()
{
  return(UCVM_CODE_SUCCESS);
}


/* Version CVM-H */
const char *ucvm_cvmh_model_version()
{
  return(ucvm_cvmh_version_id);
}


/* Query CVM-H */
int ucvm_cvmh_model_query(ucvm_ctype_t cmode,
			  int n, ucvm_point_t *pnt, ucvm_prop_t *prop)
{
  int i;
  vx_entry_t entry;
  int datagap = 0;

  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Only geo depth mode supported\n");
    return(UCVM_CODE_ERROR);
  }

  /* Query by lat/lon */
  entry.coor_type = VX_COORD_GEO;

  for (i = 0; i < n; i++) {
    if ((prop[i].model == UCVM_MODEL_NONE) && 
	(region_contains(&cvmh_region, cmode, &(pnt[i])))) {

      /* Setup point to query */
      entry.coor[0] = pnt[i].coord[0];
      entry.coor[1] = pnt[i].coord[1];
      entry.coor[2] = pnt[i].coord[2];
      
      /* Query the point */
      vx_getcoord(&entry);
      
      /* NOTE comment if and datagap to force no data return */
      if (entry.data_src != VX_SRC_NR) {
	prop[i].vp = entry.vp;
	prop[i].vs = entry.vs;
	prop[i].rho = entry.rho;
	prop[i].model = ucvm_cvmh_id;
      } else {
	datagap = 1;
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


/* Fill model structure with CVM-H */
int ucvm_cvmh_get_model(ucvm_model_t *m)
{
  strcpy(m->label, ucvm_cvmh_label_id);
  memset(&(m->region), 0, sizeof(ucvm_region_t));
  strcpy(m->config, "");
  m->init = ucvm_cvmh_model_init;
  m->finalize = ucvm_cvmh_model_finalize;
  m->version = ucvm_cvmh_model_version;
  m->query = ucvm_cvmh_model_query;

  return(UCVM_CODE_SUCCESS);
}
