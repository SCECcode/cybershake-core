#ifndef UCVM_CVMH_H
#define UCVM_CVMH_H

#include "ucvm_dtypes.h"


/* Init CVM-H */
int ucvm_cvmh_model_init(int id, ucvm_region_t *r, const char *mpath);


/* Finalize CVM-H */
int ucvm_cvmh_model_finalize();


/* Version CVM-H */
const char *ucvm_cvmh_model_version();


/* Query CVM-H */
int ucvm_cvmh_model_query(ucvm_ctype_t cmode,
			  int n, ucvm_point_t *pnt, 
			  ucvm_prop_t *prop);


/* Fill model structure with 1D */
int ucvm_cvmh_get_model(ucvm_model_t *m);


#endif
