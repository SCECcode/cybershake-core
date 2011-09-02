#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "dem.h"


/* Getopt flags */
extern char *optarg;
extern int optind, opterr, optopt;


/* Usage function */
void usage() {
  printf("\n     dem_query - (c) SCEC\n");
  printf("Extract elevation from USGS 1 arcs NED given lon/lat input coord.\n\n");
   printf("\tusage: dem_query -e neddir -b bathdir < file.in\n\n");
  printf("Flags:\n");
  printf("\t-e NED root directory. Default is ./ned.\n");
  printf("\t-b Bathymetry root directory. Default is ./bath.\n");
  printf("\t-h This help message.\n\n");
  printf("Output format is:\n");
  printf("\tlon lat elev valid\n\n");
  printf("Notes:\n");
  printf("\t- If running interactively, type Cntl-D to end input coord list.\n");
  exit (0);
}


int main(int argc, char **argv)
{
  int opt;
  char neddir[512], bathdir[512];
  dem_point_t pnt;
  dem_elev_t elev;

  strcpy(neddir, "./ned");
  strcpy(bathdir, "./bath");

  /* Parse options */
  while ((opt = getopt(argc, argv, "b:e:h")) != -1) {
    switch (opt) {
    case 'b':
      strcpy(bathdir, optarg);
      break;
    case 'e':
      strcpy(neddir, optarg);
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

  if (dem_init(neddir, bathdir) != 0) {
    fprintf(stderr, "Failed to init DEM interface\n");
    return(1);
  }

  /* Read in coords */
  while (!feof(stdin)) {
    if (fscanf(stdin,"%lf %lf", &(pnt.coord[0]), &(pnt.coord[1])) == 2) {
      if (dem_query(1, &pnt, &elev) != 0) {
	fprintf(stderr, "Failed to query point\n");
	return(1);
      }
    }

    printf("%9.5lf %9.5lf %s %9.5lf %d\n", 
	   pnt.coord[0], pnt.coord[1], elev.source, elev.elev, elev.valid);
  }

  return(0);
}
