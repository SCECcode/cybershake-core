#ifndef RUPGEN_CONFIG_H
#define RUPGEN_CONFIG_H

/* Maximum string length */
#define RG_CONFIG_MAX_STR 512

/* Config entry */
typedef struct rg_config_t {
  char name[RG_CONFIG_MAX_STR];
  char value[RG_CONFIG_MAX_STR];
  struct rg_config_t *next;
} rg_config_t;

/* Parse config file */
rg_config_t *rg_parse_config(const char *file);

/* Return next entry containing name as a key */
rg_config_t *rg_find_name(rg_config_t *chead, const char *name);

/* Dump config to screen */
int rg_dump_config(rg_config_t *chead);

/* Free config list */
int rg_free_config(rg_config_t *chead);

#endif
