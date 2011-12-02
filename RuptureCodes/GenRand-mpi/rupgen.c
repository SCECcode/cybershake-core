#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <mpi.h>
#include <unistd.h>
#include <dirent.h>
#include "structure.h"
#include "func_mpi.h"
#include "genslip.h"
#include "function.h"

/* Constants */
#define MAX_STR_ARGV 512
#define MAX_NUM_FILES 50000
#define RUP_FILE_EXT ".txt"
#define RUP_FILE_DELIM "_"


/* Parse rupture file names */
int parsename(const char *fname, int *src, int *rup)
{
  char tmpstr[MAX_FILENAME];
  char *tok;

  *src = -1;
  *rup = -1;

  memset(tmpstr, 0, MAX_FILENAME);
  strncpy(tmpstr, fname, strstr(fname, RUP_FILE_EXT) - fname);

  tok = strtok(tmpstr, RUP_FILE_DELIM);
  if (tok != NULL) {
    *src = atoi(tok);
    tok = strtok(NULL, RUP_FILE_DELIM);
    if (tok != NULL) {
      *rup = atoi(tok);
    } else {
      return(1);
    }
  } else {
    return(1);
  }

  return(0);
}


/* Get number of rupture files */
int getnumrups(const char *path, int *numrup)
{
  FILE *fp;
  char buf[MAX_FILENAME];
  char index[MAX_FILENAME];

  *numrup = 0;
  sprintf(index, "%s/index.list", path);
  fp = fopen(index, "r");
  if (fp == NULL) {
    return(1);
  }

  while (!feof(fp)) {
    if (fgets(buf, MAX_FILENAME, fp) != NULL) {
      if (strstr(buf, RUP_FILE_EXT) != NULL) {
	(*numrup)++;
      }
    }
  }
  fclose(fp);
  return(0);
}


/* Get number of rupture files */
int getnumrups2(const char *path, int *numrup)
{
  DIR *dir;
  struct dirent *ent;

  *numrup = 0;
  dir = opendir (path);
  if (dir != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      if (strstr(ent->d_name, RUP_FILE_EXT) != NULL) {
	(*numrup)++;
      }
    }
    closedir(dir);
  } else {
    /* could not open directory */
    perror ("");
    return(1);
  }
  return(0);
}


/* Get list of rupture files */
int getrups(const char *path, rg_rfile_t *rups, int maxrup, int *numrup)
{
  FILE *fp;
  char buf[MAX_FILENAME];
  char index[MAX_FILENAME];

  *numrup = 0;
  sprintf(index, "%s/index.list", path);
  fp = fopen(index, "r");
  if (fp == NULL) {
    return(1);
  }

  while (!feof(fp)) {
    if (fgets(buf, MAX_FILENAME, fp) != NULL) {
      if (strstr(buf, RUP_FILE_EXT) != NULL) {
	if (buf[strlen(buf) - 1] == '\n') {
	  buf[strlen(buf) - 1] = '\0';
	}
	strcpy(rups[*numrup].filename, buf);
	if (parsename(buf, &(rups[*numrup].src), 
		      &(rups[*numrup].rup))!= 0) {
	  fprintf(stderr, "Failed to parse rupture filename\n");
	  return(1);
	}
	rups[(*numrup)].index = *numrup;
	rups[(*numrup)].stats.numslip = 0;
	rups[(*numrup)].stats.numhypo = 0;
	(*numrup)++;
      }
    }
  }
  fclose(fp); 

  return(0);
}


/* Get list of rupture files */
int getrups2(const char *path, rg_rfile_t *rups, int maxrup, int *numrup)
{
  DIR *dir;
  struct dirent *ent;
  
  *numrup = 0;

  dir = opendir (path);
  if (dir != NULL) {
    /* print all the files and directories within directory */
    while (((ent = readdir (dir)) != NULL) && (*numrup < maxrup)) {
      if (strstr(ent->d_name, RUP_FILE_EXT) != NULL) {
	strcpy(rups[*numrup].filename, ent->d_name);
	if (parsename(ent->d_name, &(rups[*numrup].src), 
		      &(rups[*numrup].rup))!= 0) {
	  fprintf(stderr, "Failed to parse rupture filename\n");
	  return(1);
	}
	rups[(*numrup)].stats.numslip = 0;
	rups[(*numrup)].stats.numhypo = 0;
	(*numrup)++;
      }
    }
    if ((*numrup) == maxrup) {
      fprintf(stderr, "Too many files in rupture directory\n");
      return(1);
    }
    closedir(dir);
  } else {
    /* could not open directory */
    perror ("");
    return(1);
  }

  return(0);
}


/* Run Rob Graves GenSlip rupture generator */
int run_genslip(char *infile, char *outfile, char *logfile,
		rg_stats_t *stats)
{
  int j, rgargc;
  char **rgargv = NULL;

  /* Create pseudo-program args */
  rgargc = 4;
  rgargv = malloc(rgargc * sizeof(char*));
  for (j = 0; j < rgargc; j++) {
    rgargv[j] = malloc(MAX_STR_ARGV);
  }
  sprintf(rgargv[0], "%s", "rupgen");
  sprintf(rgargv[1], "infile=%s", infile);
  sprintf(rgargv[2], "outfile=%s", outfile);
  sprintf(rgargv[3], "logfile=%s", logfile);

  /* Run rupture generator */
  genslip(rgargc, rgargv, stats);
  
  /* Free pseudo-program args */
  rgargc = 4;
  for (j = 0; j < rgargc; j++) {
    free(rgargv[j]);
  }

  free(rgargv);
  return(0);
}


int main(int argc, char **argv)
{
  /* MPI stuff and distributed computation variables */
  int myid, nproc, pnlen;
  char procname[128];
  MPI_Status status;

  /* MPI Rupture info vars */
  MPI_Datatype MPI_RUP_T;
  int num_fields_rup;

  char rupdir[MAX_FILENAME];
  char logdir[MAX_FILENAME];
  char infile[MAX_FILENAME];
  char outfile[MAX_FILENAME];
  char logfile[MAX_FILENAME];
  rg_rfile_t *rups = NULL;
  int i, n, numrup, numvar;
  rg_rfile_t rup;
  int currup = 0;
  int done = 0;

  if (argc != 3) {
    printf("usage: %s rupdir logdir", argv[0]);
    return(0);
  }

  /* Init MPI */
  mpi_init(&argc, &argv, &nproc, &myid, procname, &pnlen);

  /* Register new mesh data types */
  if (mpi_register_rupinfo(&MPI_RUP_T, &num_fields_rup) != 0) {
      fprintf(stderr, "[%d] Failed to register data type\n", myid);
      return(1);
  }

  if (nproc < 2) {
      fprintf(stderr, "[%d] nproc must be at least 2\n", myid);
      return(1);
  }

  strcpy(rupdir, argv[1]);
  strcpy(logdir, argv[2]);
  numrup = 0;

  if (myid == 0) {

    printf("[%d] rupdir=%s\n", myid, rupdir);
    printf("[%d] logdir=%s\n", myid, logdir);
    fflush(stdout);

    /* Get number of rupture files */
    if (getnumrups(rupdir, &numrup) != 0) {
      fprintf(stderr, "[%d] Failed to count rupture files\n", myid);
      return(1);
    }

    if (numrup == 0) {
      fprintf(stderr, "[%d] No rupture files to process\n", myid);
      return(1);
    } else {
      printf("[%d] Retrieved %d rupture files\n", myid, numrup);
    }

    /* Round num ruptures up to nearest multiple of nproc */
    rups = malloc(numrup * sizeof(rg_rfile_t));
    if (rups == NULL) {
      fprintf(stderr, "[%d] Failed to allocate rupture file buffer\n",
	      myid);
      return(1);
    }
    memset(rups, 0, sizeof(rg_rfile_t)*numrup);

    /* Get rupture files */
    if (getrups(rupdir, rups, numrup, &n) != 0) {
      fprintf(stderr, "[%d] Failed to find rupture files\n", myid);
      return(1);
    }
  }

  /* Main loop */
  memset(&rup, 0, sizeof(rg_rfile_t));
  if (myid == 0) {
    /* Dispatch rupture files to worker pool */
    while (currup < numrup) {
      //printf("[%d] currup=%d, numrup=%d\n", myid, currup, numrup);
      /* Blocking receive on worker pool */
      MPI_Recv(&rup, 1, MPI_RUP_T, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD, &status);

      /* Update rupture record if necessary */
      if (strlen(rup.filename) > 0) {
	memcpy(&(rups[rup.index]), &rup, sizeof(rg_rfile_t));
      }

      /* Send next available rupture to worker */
      MPI_Send(&(rups[currup++]), 1, MPI_RUP_T, status.MPI_SOURCE, 
	       0, MPI_COMM_WORLD);
    }

    /* Send stop message to all workers */
    printf("[%d] All ruptures dispatched, stopping workers\n", myid);
    for (i = 1; i < nproc; i++) {
      /* Blocking receive on worker pool */
      MPI_Recv(&rup, 1, MPI_RUP_T, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	       MPI_COMM_WORLD, &status);

      /* Update rupture record if necessary */
      if (strlen(rup.filename) > 0) {
	memcpy(&(rups[rup.index]), &rup, sizeof(rg_rfile_t));
      }

      /* Send stop work to worker */
      memset(&rup, 0, sizeof(rg_rfile_t));
      MPI_Send(&(rup), 1, MPI_RUP_T, status.MPI_SOURCE, 
	       0, MPI_COMM_WORLD);
    }
  } else {
    while (!done) {
      /* Send work request to rank 0 */
      MPI_Send(&(rup), 1, MPI_RUP_T, 0, 0, MPI_COMM_WORLD);

      /* Blocking receive on next rupture */
      memset(&rup, 0, sizeof(rg_rfile_t));
      MPI_Recv(&rup, 1, MPI_RUP_T, 0, MPI_ANY_TAG, 
	       MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      /* Check to see if we are done */
      if (strlen(rup.filename) == 0) {
	done = 1;
      } else {
	sprintf(infile, "%s/%d/%d/%s", rupdir, rup.src,
		rup.rup, rup.filename);

	printf("[%d] Processing rupture %s\n", myid, infile);

	/* Create log directory */
	sprintf(logfile, "%s/%d", logdir, rup.src);
	mkdir(logfile, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	sprintf(logfile, "%s/%d/rupgen_%d_%d.log", logdir, rup.src,
		rup.src, rup.rup);

	/* Generate the rupture variations */
	sprintf(outfile, "%s/%d/%d/%s.variation", rupdir, 
		rup.src, rup.rup, rup.filename);
	//printf("[%d] Processing file %s\n", myid, infile);
	if (run_genslip(infile, outfile, logfile, &(rup.stats)) != 0) {
	  fprintf(stderr, "Failed to run rupture generator\n");
	  return(1);
	}
      }
    }
  }

  /* Print summary statistics */
  if (myid == 0) {
    numvar = 0;
    for (i = 0; i < numrup; i++) {
      numvar += rups[i].stats.numslip * rups[i].stats.numhypo;
    }
    printf("[%d] Ruptures: %d\n", myid, numrup);
    printf("[%d] Rupture Variations: %d\n", myid, numvar);

    free(rups);
  }

  /* Final sync */
  mpi_barrier();
  mpi_final("MPI Done");

  return(0);
}
