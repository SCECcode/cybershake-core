#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>
//#include "proj_api.h"
#include "dem_config.h"
#include "dem.h"

/* Maximum number of dems to allow */
#define DEM_MAX_DEMS 100


/* Default NED/BATH projections */
//#define PROJ_GEO "+proj=latlong +datum=WGS84"
//#define PROJ_NED "+proj=latlong +datum=NAD83 +ellps=GRS80"


/* NED header list */
int num_neds = 0;
dem_info_t nedlist[DEM_MAX_DEMS];

/* BATH header list */
int num_baths = 0;
dem_info_t bathlist[DEM_MAX_DEMS];

/* Currently open DEMs */
dem_info_t *ned_hdr = NULL;
FILE *ned_fp = NULL;
dem_info_t *bath_hdr = NULL;
FILE *bath_fp = NULL;


/* Determine system endianness */
int systemEndian()
{
  int num = 1;
  if(*(char *)&num == 1) {
    return DEM_BYTEORDER_LSB;
  } else {
    return DEM_BYTEORDER_MSB;
  }
}


/* Check if file exists */
int fileExists(const char *file)
{
  struct stat st;

  if (stat(file, &st) == 0) {
    return(1);
  } else {
    return(0);
  }
}


/* Read the header file */
int dem_read_hdr(const char *hfile, dem_info_t *hdr)
{
  dem_config_t *dem_cfg = NULL;
  dem_config_t *cfgentry = NULL;

  /* Read in the header */
  dem_cfg = dem_parse_config(hfile);
  if (dem_cfg == NULL) {
    fprintf(stderr, "Failed to read config file\n");
    return(1);
  }

  cfgentry = dem_find_name(dem_cfg, "ncols");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: ncols not found\n", hfile);
    return(1);
  }
  hdr->dims[0] = atoi(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "nrows");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: nrows not found\n", hfile);
    return(1);
  }
  hdr->dims[1] = atoi(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "xllcorner");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: xllcorner not found\n", hfile);
    return(1);
  }
  hdr->llcorner[0] = atof(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "yllcorner");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: yllcorner not found\n", hfile);
    return(1);
  }
  hdr->llcorner[1] = atof(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "cellsize");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: cellsize not found\n", hfile);
    return(1);
  }
  hdr->spacing = atof(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "nodata_value");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: NODATA_value not found\n", hfile);
    return(1);
  }
  hdr->no_data = atof(cfgentry->value);

  cfgentry = dem_find_name(dem_cfg, "byteorder");
  if (cfgentry == NULL) {
    fprintf(stderr, "%s: byteorder not found\n", hfile);
    return(1);
  }
  if (strcmp(cfgentry->value, "LSBFIRST") == 0) {
    hdr->byteorder = DEM_BYTEORDER_LSB;
  } else if (strcmp(cfgentry->value, "MSBFIRST") == 0) {
    hdr->byteorder = DEM_BYTEORDER_MSB;
  } else {
    fprintf(stderr, "Invalid byte order '%s'\n", cfgentry->value);
    return(1);
  }

  dem_free_config(dem_cfg);

  /* Compute and save upper corner */
  hdr->urcorner[0] = hdr->llcorner[0] + hdr->dims[0]*hdr->spacing;
  hdr->urcorner[0] = hdr->llcorner[1] + hdr->dims[1]*hdr->spacing;

  /* Save name of NED GridFloat formatted file */
  strcpy(hdr->file, hfile);
  strcpy(strstr(hdr->file, "hdr"), "flt");
  if (!fileExists(hdr->file)) {
    fprintf(stderr, "GridFloat file %s does not exist\n", hdr->file);
    return(1);
  }

  return(0);
}


/* Initializer */
int dem_init(const char *neddir, const char *bathdir)
{
  DIR *dir;
  struct dirent *ent;
  char hfile[512];

  /* Get list of NED files from directory */
  dir = opendir (neddir);
  if (dir == NULL) {
    fprintf(stderr, "Failed to open neddir %s\n", neddir);
    return(1);
  }

  fprintf(stderr, "Reading NED DEMs:\n");
  while ((ent = readdir(dir)) != NULL) {
    if (strstr(ent->d_name, ".hdr") != NULL) {
      sprintf(hfile, "%s/%s", neddir, ent->d_name);
      if (dem_read_hdr(hfile, &(nedlist[num_neds])) != 0) {
	fprintf(stderr, "Failed to read in header %s\n", hfile);
	return(1);
      }

      fprintf(stderr, "NED %s: %d, %d, %lf, %lf, %lf, %lf %d %s\n",
	     hfile, 
	     nedlist[num_neds].dims[0], nedlist[num_neds].dims[1], 
	     nedlist[num_neds].llcorner[0], nedlist[num_neds].llcorner[1], 
	     nedlist[num_neds].spacing, nedlist[num_neds].no_data,
	     nedlist[num_neds].byteorder, nedlist[num_neds].file);

      num_neds++;
    }

  }   
  closedir(dir);

  if (num_neds == 0) {
    fprintf(stderr, "No NEDs found in %s\n", neddir);
    return(0);
  }

  /* Get list of Bath files from directory */
  dir = opendir (bathdir);
  if (dir == NULL) {
    fprintf(stderr, "Failed to open bathdir %s\n", bathdir);
    return(1);
  }

  fprintf(stderr, "Reading Bathymetry:\n");
  while ((ent = readdir(dir)) != NULL) {
    if (strstr(ent->d_name, ".hdr") != NULL) {
      sprintf(hfile, "%s/%s", bathdir, ent->d_name);
      if (dem_read_hdr(hfile, &(bathlist[num_baths])) != 0) {
	fprintf(stderr, "Failed to read in header %s\n", hfile);
	return(1);
      }
      fprintf(stderr, "BAT %s: %d, %d, %lf, %lf, %lf, %lf %d %s\n",
	     hfile, 
	     bathlist[num_baths].dims[0], bathlist[num_baths].dims[1], 
	     bathlist[num_baths].llcorner[0], bathlist[num_baths].llcorner[1], 
	     bathlist[num_baths].spacing, bathlist[num_baths].no_data,
	     bathlist[num_baths].byteorder, bathlist[num_baths].file);
      num_baths++;
    }
  }   
  closedir(dir);

  if (num_baths == 0) {
    fprintf(stderr, "Warning: No BATs found in %s\n", bathdir);
    return(0);
  }

  return(0);
}


/* Find header that contains point of interest */
int dem_find_hdr(dem_src_t src, dem_point_t *p,
		 dem_info_t **hdr, FILE **fp, int *x, int *y)
{
  int i;
  int num_hdr;
  dem_info_t *thdr = NULL;
  dem_info_t *hlist = NULL;
  dem_info_t *oldhdr = NULL;
  FILE *oldfp;
  dem_point_t tp;

  *x = -1;
  *y = -1;
  oldhdr = *hdr;
  oldfp = *fp;
  *hdr = NULL;
  *fp = NULL;

  switch (src) {
  case DEM_SRC_NED:
    num_hdr = num_neds;
    hlist = nedlist;
    break;
  case DEM_SRC_BATH:
    num_hdr = num_baths;
    hlist = bathlist;
    break;
  default:
    return(1);
  }

  /* Convert lon,lat datums */
  tp.coord[0] = p->coord[0];
  tp.coord[1] = p->coord[1];

  /* Find header that contains this point */
  for (i = 0; i < num_hdr; i++) {
    thdr = &(hlist[i]);
    *x = (int)((tp.coord[0] - thdr->llcorner[0]) / thdr->spacing);
    *y = (int)((tp.coord[1] - thdr->llcorner[1]) / thdr->spacing);
    if ((*x >= 0) && (*y >= 0) && 
	(*x < thdr->dims[0]) && (*y < thdr->dims[1])) {
      *hdr = thdr;
      break;
    }
  }

  /* Found a header */
  if (*hdr != NULL) {
    if (oldhdr != *hdr) {
      /* Close old data file */
      if (oldfp != NULL) {
	fclose(oldfp);
      }

      /* Open new data file */
      *fp = fopen((*hdr)->file, "rb");
      if (*fp == NULL) {
	return(1);
      }
    } else {
      *fp = oldfp;
    }
  }

  return(0);
}


/* Get elevation at offset x,y from current file */
int dem_get_elev(dem_info_t *hdr, FILE *fp, int x, int y, float *elev)
{
  size_t offset;

  *elev = 0.0;

  /* Handle endian-ness here */
  if (hdr->byteorder != systemEndian()) {
    fprintf(stderr, "Unsupported byte order, add byte swap here\n");
    return(1);
  }

  /* Find elevation at this point */
  offset = ((hdr->dims[1]-y-1) * hdr->dims[0] + x) * 
    sizeof(float);
  if (fseek(fp, offset, SEEK_SET) != 0) {
    fprintf(stderr, "Failed to seek within DEM %s\n", hdr->file);
    return(1);
  }
  if (fread(elev, sizeof(float), 1, fp) != 1) {
    fprintf(stderr, "Failed to read from DEM %s\n", hdr->file);
    return(1);
  }

  return(0);
}


/* Query underlying models */
int dem_query(int n, dem_point_t *pnt, dem_elev_t *elev)
{
  int i;
  int ned_x, ned_y, bath_x, bath_y;
  float ned_val, bath_val;

  for (i = 0; i < n; i++) {
    /* Find NED/Bath files containing point */
    strcpy(elev[i].source, "none");
    elev[i].elev = 0.0;
    elev[i].valid = 0;
    if (dem_find_hdr(DEM_SRC_NED, &(pnt[i]), 
		     &ned_hdr, &ned_fp, &ned_x, &ned_y) != 0) {
      fprintf(stderr, "Failed to find NED header\n");
      return(1);
    }
    if (dem_find_hdr(DEM_SRC_BATH, &(pnt[i]), 
		     &bath_hdr, &bath_fp, &bath_x, &bath_y) != 0) {
      fprintf(stderr, "Failed to find BATH header\n");
      return(1);
    }

    /* Save NED data if valid */
    if (ned_hdr != NULL) {
      dem_get_elev(ned_hdr, ned_fp, ned_x, ned_y, &ned_val);
      if (fabs(ned_val - ned_hdr->no_data) > 0.01) {
	strcpy(elev[i].source, "ned");
	elev[i].elev = (double)ned_val;
	elev[i].valid = 1;
      }
    }

    /* If NED is 0.0 or invalid, check BATH data and save it if 
       it is valid and negative */
    if ((elev[i].valid == 0) || (fabs(elev[i].elev - 0.0) < 0.01)) {
      if (bath_hdr != NULL) {
	dem_get_elev(bath_hdr, bath_fp, bath_x, bath_y, &bath_val);
	if ((bath_val < 0.0) && (fabs(bath_val - bath_hdr->no_data) > 0.01)) {
	  strcpy(elev[i].source, "bath");
	  elev[i].elev = (double)bath_val;
	  elev[i].valid = 1;
	}
      }
    }
  }

  return(0);
}

