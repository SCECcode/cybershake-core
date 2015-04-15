#ifndef UCVM_GRID_H
#define UCVM_GRID_H

#include "ucvm_dtypes.h"

/* Generate grid from projection and dimensions */
int ucvm_grid_gen(ucvm_proj_t *iproj, ucvm_trans_t *trans,
		  ucvm_proj_t *oproj,
		  ucvm_dim_t *dims, double spacing, 
		  ucvm_point_t *pnts);


/* Generate grid from projection and dimensions */
int ucvm_grid_gen_file(ucvm_proj_t *iproj, ucvm_trans_t *trans,
		       ucvm_proj_t *oproj,
		       ucvm_dim_t *dims, double spacing, 
		       const char *filename);


/* Convert point list from one projection to another */
int ucvm_grid_convert(ucvm_proj_t *iproj, 
		      ucvm_proj_t *oproj, 
		      size_t n, ucvm_point_t *pnts);


/* Convert point list from one projection to another */
int ucvm_grid_convert_file(ucvm_proj_t *iproj, 
			   ucvm_proj_t *oproj, 
			   size_t n, const char *filename);


#endif
