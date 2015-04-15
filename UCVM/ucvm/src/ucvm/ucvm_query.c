#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "ucvm.h"

/* Constants */
#define STRTOK_DELIM ","
#define NUM_POINTS 20000


/* Getopt flags */
extern char *optarg;
extern int optind, opterr, optopt;


/* Usage function */
void usage() {
  printf("\n     ucvm_query - (c) SCEC\n");
  printf("Extract velocities from any Community Velocity Model. Accepts\n");
  printf("geographic coordinates and UTM Zone 11, NAD27 coordinates in\n");
  printf("X Y Z columns.\n\n");
  printf("\tusage: ucvm_query [-m models] [-c coordtype] [-f config] < file.in\n\n");
  printf("Flags:\n");
  printf("\t-m Models to query in order of preference: cvms, cvmh, cencal, 1d (default).\n");
  printf("\t-c Coordinate mode: geo-depth (gd, default), geo-elev (ge), utm-depth (ud), utm-elev (ue).\n");
  printf("\t-f Configuration file. Default is ./ucvm.conf.\n");
  printf("\t-h This help message.\n\n");
  printf("Output format is:\n");
  printf("\tX Y Z model vp vs rho\n\n");
  printf("Notes:\n");
  printf("\t- If running interactively, type Cntl-D to end input coord list.\n");
  exit (0);
}


int main(int argc, char **argv)
{
  int i;
  int opt;
  char modelstr[UCVM_MAX_LABEL_LEN];
  char configfile[UCVM_MAX_PATH_LEN];
  char *token = NULL;

  ucvm_ctype_t cmode;
  int num_models = 0;
  char models[UCVM_MAX_MODELS][UCVM_MAX_LABEL_LEN];

  ucvm_point_t *pnt1;
  ucvm_prop_t *prop1;
  int numread = 0;
  char label[UCVM_MAX_LABEL_LEN];


  cmode = UCVM_COORD_GEO_DEPTH;
  strcpy(configfile, "./ucvm.conf");

  /* Parse options */
  while ((opt = getopt(argc, argv, "m:c:f:h")) != -1) {
    switch (opt) {
    case 'm':
      if (strlen(optarg) < UCVM_MAX_LABEL_LEN) {
	strcpy(modelstr, optarg);
	token = strtok(modelstr, STRTOK_DELIM);
	while (token != NULL) {
	  if (num_models == UCVM_MAX_MODELS) {
	      fprintf(stderr, "Maximum number of models exceeded\n");
	      return(1);
	  }
	  strcpy(models[num_models], token);
	  for (i = 0; i < num_models; i++) {
	    if (strcmp(models[i], token) == 0) {
	      fprintf(stderr, "Repeated model: %s.\n", token);
	      usage();
	      exit(1);
	    }
	  }
	  num_models++;
	  token = strtok(NULL, STRTOK_DELIM);
	}
      } else {
	fprintf(stderr, "Invalid model list specified: %s.\n", optarg);
	usage();
	exit(1);
      }
      break;
    case 'c':
      if (strcmp(optarg, "gd") == 0) {
	cmode = UCVM_COORD_GEO_DEPTH;
      } else if (strcmp(optarg, "ge") == 0) {
	cmode = UCVM_COORD_GEO_ELEV;
      } else if (strcmp(optarg, "ud") == 0) {
	cmode = UCVM_COORD_UTM_DEPTH;
      } else if (strcmp(optarg, "ue") == 0) {
	cmode = UCVM_COORD_UTM_ELEV;
      } else {
	fprintf(stderr, "Invalid z mode %s.\n", optarg);
	usage();
	exit(1);
      }
      break;
    case 'f':
      if (strlen(optarg) < UCVM_MAX_PATH_LEN) {
	strcpy(configfile, optarg);
      } else {
	fprintf(stderr, "Invalid config file specified: %s.\n", optarg);
	usage();
	exit(1);
      }
      break;
    case 'h':
      usage();
      exit(0);
      break;
    default: /* '?' */
      usage();
      exit(1);
    }
  }

  if (num_models == 0) {
    fprintf(stderr, "Using 1D as default model.\n");
    if (ucvm_add_model(UCVM_MODEL_1D) != UCVM_CODE_SUCCESS) {
      fprintf(stderr, "Enable model failed for model 1d\n");
      return(1);
    }
  }

  if (cmode != UCVM_COORD_GEO_DEPTH) {
    fprintf(stderr, "Unsupported coordinate mode.\n");
    usage();
    return(1);
  }
  fprintf(stderr, "Using Geo Depth coordinates as default mode.\n");

  if (ucvm_init(configfile) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Init failed\n");
    return(1);
  }

  if (ucvm_query_mode(cmode) != UCVM_CODE_SUCCESS) {
    fprintf(stderr, "Set query mode failed\n");
    return(1);
  }

  for (i = 0; i < num_models; i++) {
    if (ucvm_add_model(models[i]) != UCVM_CODE_SUCCESS) {
      fprintf(stderr, "Enable model failed for model %s\n", models[i]);
      return(1);
    }
  }

  pnt1 = malloc(NUM_POINTS * sizeof(ucvm_point_t));
  prop1 = malloc(NUM_POINTS * sizeof(ucvm_prop_t));

  /* Read in coords */
  while (!feof(stdin)) {
    if (fscanf(stdin,"%lf %lf %lf",
               &(pnt1[numread].coord[0]),
	       &(pnt1[numread].coord[1]),
	       &(pnt1[numread].coord[2])) == 3) {

      numread++;
      if (numread == NUM_POINTS) {
	/* Query the UCVM */
	if (ucvm_query(numread, pnt1, prop1) != UCVM_CODE_SUCCESS) {
	  fprintf(stderr, "Query CVM failed\n");
	  return(1);
	}

	/* Display results */
	for (i = 0; i < numread; i++) {
	  ucvm_model_label(prop1[i].model, label, UCVM_MAX_LABEL_LEN);
	  printf("%lf %lf %lf %s %lf %lf %lf\n", 
		 pnt1[i].coord[0], pnt1[i].coord[1], pnt1[i].coord[2],
		 label, prop1[i].vp, prop1[i].vs, prop1[i].rho);
	}

	numread = 0;
      }
    }
  }

  if (numread > 0) {
    /* Query the UCVM */
    if (ucvm_query(numread, pnt1, prop1) != UCVM_CODE_SUCCESS) {
      fprintf(stderr, "Query CVM failed\n");
      return(1);
    }
    
    /* Display results */
    for (i = 0; i < numread; i++) {
      ucvm_model_label(prop1[i].model, label, UCVM_MAX_LABEL_LEN);
      printf("%lf %lf %lf %s %lf %lf %lf\n", 
	     pnt1[i].coord[0], pnt1[i].coord[1], pnt1[i].coord[2],
	     label, prop1[i].vp, prop1[i].vs, prop1[i].rho);
    }
    
    numread = 0;
  }

  ucvm_finalize();
  free(pnt1);
  free(prop1);

  return(0);
}
