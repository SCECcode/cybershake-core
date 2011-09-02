#ifndef UCVM_CENCAL_H
#define UCVM_CENCAL_H

#include "ucvm_dtypes.h"


/* Init CenCal */
int ucvm_cencal_model_init(int id, ucvm_region_t *r, const char *mpath);


/* Finalize CenCal */
int ucvm_cencal_model_finalize();


/* Version CenCal */
const char *ucvm_cencal_model_version();


/* Query CenCal */
int ucvm_cencal_model_query(ucvm_ctype_t cmode,
			    int n, ucvm_point_t *pnt, 
			    ucvm_prop_t *prop);


/* Fill model structure with CenCal */
int ucvm_cencal_get_model(ucvm_model_t *m);


#endif
