#ifndef DEM_CONFIG_H
#define DEM_CONFIG_H

/* Maximum string length */
#define DEM_CONFIG_MAX_STR 512

/* Config entry */
typedef struct dem_config_t {
  char name[DEM_CONFIG_MAX_STR];
  char value[DEM_CONFIG_MAX_STR];
  struct dem_config_t *next;
} dem_config_t;

/* Parse config file */
dem_config_t *dem_parse_config(const char *file);

/* Search for a name in the list */
dem_config_t *dem_find_name(dem_config_t *chead, const char *name);

/* Free config list */
int dem_free_config(dem_config_t *chead);

#endif
