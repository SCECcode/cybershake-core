#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include "dem_config.h"

/* Whitespace characters */
const char *DEM_WHITESPACE = " \t\r\n";
const char *DEM_SEPARATOR = " ";


/* Strip whitespace from string */
void dem_strip_whitespace(char *str)
{
  int i1, i2;
  int len;

  i1 = 0;
  i2 = 0;
  len = strlen(str);

  for (i2 = 0; i2 < len; i2++) {
    if (strchr(DEM_WHITESPACE, str[i2]) == NULL) {
      str[i1++] = str[i2];
    }
  }
  str[i1] = '\0';
  return;
}


/* Strip trailing whitespace from string */
void dem_strip_trailing_whitespace(char *str)
{
  int i;

  i = strlen(str);
  while (strchr(DEM_WHITESPACE, str[i-1]) != NULL) {
    str[i-1] = '\0';
    i = i - 1;
  }
  return;
}


/* Convert to lowercase */
void dem_lowercase(char *str)
{
  int i;

  for(i = strlen(str); i >= 0; i--) {
    str[i] = tolower(str[i]);
  }
  
  return;
}


/* Parse config file */
dem_config_t *dem_parse_config(const char *file)
{
  FILE *fp;
  char line[DEM_CONFIG_MAX_STR];
  char *token;
  dem_config_t celem;
  dem_config_t *chead = NULL;
  dem_config_t *cnew;

  fp = fopen(file, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed to open config %s\n", file);
    return(NULL);
  }

  while (!feof(fp)) {
    fgets(line, DEM_CONFIG_MAX_STR, fp);
    if (strlen(line) > 0) {
      token = strtok(line, DEM_SEPARATOR);
      if (token == NULL) {
	continue;
      }
      dem_lowercase(token);
      strcpy(celem.name, token);
      dem_strip_whitespace(celem.name);
      if (celem.name[0] == '#') {
	continue;
      }
      strcpy(celem.value, &(line[strlen(token) + 1]));
      dem_strip_whitespace(celem.value);
      cnew = (dem_config_t *)malloc(sizeof(dem_config_t));
      memcpy(cnew, &celem, sizeof(dem_config_t));
      cnew->next = chead;
      chead = cnew;
    }

  }

  fclose(fp);
  return(chead);
}


/* Search for a name in the list */
dem_config_t *dem_find_name(dem_config_t *chead, const char *name)
{
  dem_config_t *cptr;

  cptr = chead;
  while (cptr != NULL) {
    if (strcmp(name, cptr->name) == 0) {
      break;
    }
    cptr = cptr->next;
  }

  return(cptr);
}


/* Free config list */
int dem_free_config(dem_config_t *chead)
{
  dem_config_t *cptr;
  dem_config_t *ctmp;

  cptr = chead;
  while (cptr != NULL) {
    ctmp = cptr;
    cptr = cptr->next;
    free(ctmp);
  }

  return(0);
}

