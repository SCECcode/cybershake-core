#ifndef UCVM_CVMS_H
#define UCVM_CVMS_H

#include "ucvm_dtypes.h"


/* Init CVM-S */
int ucvm_cvms_model_init(int id, ucvm_region_t *r, const char *mpath);


/* Finalize CVM-S */
int ucvm_cvms_model_finalize();


/* Version CVM-S */
const char *ucvm_cvms_model_version();


/* Query CVM-S */
int ucvm_cvms_model_query(ucvm_ctype_t cmode,
			  int n, ucvm_point_t *pnt, 
			  ucvm_prop_t *prop);


/* Fill model structure with CVM-S */
int ucvm_cvms_get_model(ucvm_model_t *m);

#endif
