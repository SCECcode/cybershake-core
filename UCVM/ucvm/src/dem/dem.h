#ifndef DEM_H
#define DEM_H

#include "dem_dtypes.h"


/* Initializer */
int dem_init(const char *neddir, const char *bathdir);


/* Query underlying models */
int dem_query(int n, dem_point_t *pnt, dem_elev_t *elev);


#endif
