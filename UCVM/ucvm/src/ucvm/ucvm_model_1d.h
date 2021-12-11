#ifndef UCVM_1D_H
#define UCVM_1D_H

#include "ucvm_dtypes.h"


/* Init 1D */
int ucvm_1d_model_init(int id, ucvm_region_t *r, const char *mpath);


/* Finalize 1D */
int ucvm_1d_model_finalize();


/* Version 1D */
const char *ucvm_1d_model_version();


/* Query 1D */
int ucvm_1d_model_query(ucvm_ctype_t cmode,
			int n, ucvm_point_t *pnt, ucvm_prop_t *prop);


/* Fill model structure with 1D */
int ucvm_1d_get_model(ucvm_model_t *m);


#endif
