#ifndef UCVM_H
#define UCVM_H

#include "ucvm_dtypes.h"

/* Predefined models */
#define UCVM_MODEL_CVMH "cvmh"
#define UCVM_MODEL_CVMS "cvms"
#define UCVM_MODEL_CENCAL "cencalvm"
#define UCVM_MODEL_1D "1d"
#define UCVM_MODEL_LINTHURBER "lin-thurber"


/* Initializer */
int ucvm_init(const char *config);


/* Finalizer */
int ucvm_finalize();


/* Set query mode */
int ucvm_query_mode(ucvm_ctype_t c);


/* Enable specific model, by label or by ucvm_model_t */
int ucvm_add_model(const char *label);
int ucvm_add_user_model(ucvm_model_t *m);


/* Get label for a model */
int ucvm_model_label(int m, char *label, int len);


/* Get version for a model */
int ucvm_model_version(int m, char *ver, int len);


/* Query underlying models */
int ucvm_query(int n, ucvm_point_t *pnt, ucvm_prop_t *prop);


#endif
