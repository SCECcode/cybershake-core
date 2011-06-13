#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ucvm_config.h"

/* Whitespace characters */
const char *UCVM_WHITESPACE = " \t\n";


/* Strip whitespace from string */
void ucvm_strip_whitespace(char *str)
{
  int i1, i2;
  int len;

  i1 = 0;
  i2 = 0;
  len = strlen(str);

  for (i2 = 0; i2 < len; i2++) {
    if (strchr(UCVM_WHITESPACE, str[i2]) == NULL) {
      str[i1++] = str[i2];
    }
  }
  str[i1] = '\0';
  return;
}


/* Strip trailing whitespace from string */
void ucvm_strip_trailing_whitespace(char *str)
{
  int i;

  i = strlen(str);
  while (strchr(UCVM_WHITESPACE, str[i-1]) != NULL) {
    str[i-1] = '\0';
    i = i - 1;
  }
  return;
}


/* Parse config file */
ucvm_config_t *ucvm_parse_config(const char *file)
{
  FILE *fp;
  char line[UCVM_CONFIG_MAX_STR];
  char *token;
  ucvm_config_t celem;
  ucvm_config_t *chead = NULL;
  ucvm_config_t *cnew;

  fp = fopen(file, "r");
  if (fp == NULL) {
    fprintf(stderr, "Failed to open config %s\n", file);
    return(NULL);
  }

  while (!feof(fp)) {
    fgets(line, UCVM_CONFIG_MAX_STR, fp);
    if (strlen(line) > 0) {
      token = strtok(line, "=");
      if (token == NULL) {
	continue;
      }
      strcpy(celem.name, token);
      ucvm_strip_whitespace(celem.name);
      if (celem.name[0] == '#') {
	continue;
      }
      strcpy(celem.value, &(line[strlen(token) + 1]));
      ucvm_strip_trailing_whitespace(celem.value);
      cnew = (ucvm_config_t *)malloc(sizeof(ucvm_config_t));
      memcpy(cnew, &celem, sizeof(ucvm_config_t));
      cnew->next = chead;
      chead = cnew;
    }

  }

  fclose(fp);
  return(chead);
}


/* Search for a name in the list */
ucvm_config_t *ucvm_find_name(ucvm_config_t *chead, const char *name)
{
  ucvm_config_t *cptr;

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
int ucvm_free_config(ucvm_config_t *chead)
{
  ucvm_config_t *cptr;
  ucvm_config_t *ctmp;

  cptr = chead;
  while (cptr != NULL) {
    ctmp = cptr;
    cptr = cptr->next;
    free(ctmp);
  }

  return(0);
}

