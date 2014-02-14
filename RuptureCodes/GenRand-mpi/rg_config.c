#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "rg_config.h"
#include "function.h"


/* Whitespace characters */
const char *RG_WHITESPACE = " \t\n";


/* Strip whitespace from string */
void rg_strip_whitespace(char *str)
{
  int i1, i2;
  int len;

  i1 = 0;
  i2 = 0;
  len = strlen(str);

  for (i2 = 0; i2 < len; i2++) {
    if (strchr(RG_WHITESPACE, str[i2]) == NULL) {
      str[i1++] = str[i2];
    }
  }
  str[i1] = '\0';
  return;
}


/* Strip trailing whitespace from string */
void rg_strip_trailing_whitespace(char *str)
{
  int i;

  i = strlen(str);
  while (strchr(RG_WHITESPACE, str[i-1]) != NULL) {
    str[i-1] = '\0';
    i = i - 1;
  }
  return;
}


/* Parse config file */
rg_config_t *rg_parse_config(const char *file)
{
  FILE *fp;
  char line[RG_CONFIG_MAX_STR];
  char *name, *value;
  rg_config_t celem;
  rg_config_t *chead = NULL;
  rg_config_t *cnew;

  if (!rg_is_file(file)) {
    fprintf(stderr, "Config file %s is not a valid file\n", file);
    return(NULL);
  }

  fp = fopen(file, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed to open config %s\n", file);
    return(NULL);
  }

  while (!feof(fp)) {
    if ((fgets(line, RG_CONFIG_MAX_STR, fp) != NULL) && 
	(strlen(line) > 0)) {
      memset(&celem, 0, sizeof(rg_config_t));
      name = strtok(line, "=");
      if (name == NULL) {
	continue;
      }
      rg_strcpy(celem.name, name, RG_CONFIG_MAX_STR);
      rg_strip_whitespace(celem.name);
      if ((celem.name[0] == '#') || (strlen(celem.name) == 0)) {
        continue;
      }
      value = name + strlen(name) + 1;
      if (strlen(value) == 0) {
	continue;
      }
      rg_strcpy(celem.value, value, RG_CONFIG_MAX_STR);
      rg_strip_trailing_whitespace(celem.value);
      cnew = (rg_config_t *)malloc(sizeof(rg_config_t));
      memcpy(cnew, &celem, sizeof(rg_config_t));
      cnew->next = chead;
      chead = cnew;
    }
  }

  fclose(fp);
  return(chead);
}


/* Search for a name in the list */
rg_config_t *rg_find_name(rg_config_t *chead, const char *name)
{
  rg_config_t *cptr;

  cptr = chead;
  while (cptr != NULL) {
    if (strcmp(name, cptr->name) == 0) {
      break;
    }
    cptr = cptr->next;
  }

  return(cptr);
}


/* Dump config to screen */
int rg_dump_config(rg_config_t *chead)
{
  rg_config_t *cptr;

  cptr = chead;
  fprintf(stderr, "Parsed Config:\n");
  while (cptr != NULL) {
    fprintf(stderr, "\t%s : %s\n", cptr->name, cptr->value);
    cptr = cptr->next;
  }
  return(0);
}


/* Free config list */
int rg_free_config(rg_config_t *chead)
{
  rg_config_t *cptr;
  rg_config_t *ctmp;

  cptr = chead;
  while (cptr != NULL) {
    ctmp = cptr;
    cptr = cptr->next;
    free(ctmp);
  }

  return(0);
}

