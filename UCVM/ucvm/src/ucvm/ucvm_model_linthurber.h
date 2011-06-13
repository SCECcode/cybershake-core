#ifndef UCVM_LINTHURBER_H
#define UCVM_LINTHURBER_H

#include "ucvm_dtypes.h"


/* Init LinThurber */
int ucvm_linthurber_model_init(int id, ucvm_region_t *r, const char *mpath);


/* Finalize LinThurber */
int ucvm_linthurber_model_finalize();


/* Version LinThurber */
const char *ucvm_linthurber_model_version();


/* Query LinThurber */
int ucvm_linthurber_model_query(ucvm_ctype_t cmode,
			int n, ucvm_point_t *pnt, ucvm_prop_t *prop);


/* Fill model structure with LinThurber */
int ucvm_linthurber_get_model(ucvm_model_t *m);


#endif
