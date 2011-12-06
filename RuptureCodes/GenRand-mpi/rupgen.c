#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <mpi.h>
#include <unistd.h>
#include <dirent.h>
#include <libgen.h>
#include "structure.h"
#include "func_mpi.h"
#include "genslip.h"
#include "function.h"
#include "rg_config.h"


/* Constants */
#define MAX_STR_ARGV 512
#define RUP_FILE_EXT ".txt"
#define NUM_GENSLIP_ARGS 7
#define GENSLIP_WILDCARD -1
#define GENSLIP_NOWRITE 0
#define GENSLIP_WRITE 1


/* MPI Rupture info vars */
MPI_Datatype MPI_RUP_T;
int num_fields_rup;
MPI_Datatype MPI_CONF_T;
int num_fields_conf;


/* Rupture comparator, sort by slip descending order */
int rupcomp(const void *p1, const void *p2)
{
  rg_rfile_t *r1;
  rg_rfile_t *r2;

  r1 = (rg_rfile_t *)p1;
  r2 = (rg_rfile_t *)p2;

  /* Reverse sort by slip number */
  if (r1->stats.numslip < r2->stats.numslip) {
    return(1);
  } else if (r1-> mag == r2->mag) {
    return(0);
  } else {
    return(-1);
  }
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


/* Get list of rupture files */
int getrups(const char *path, rg_rfile_t *rups, int maxrup, 
	    int *numrup)
{
  FILE *fp;
  char buf[MAX_FILENAME];
  char index[MAX_FILENAME];
  int i;

  *numrup = 0;
  sprintf(index, "%s/index.list", path);
  fp = fopen(index, "r");
  if (fp == NULL) {
    return(1);
  }

  while (!feof(fp)) {
    if (fgets(buf, MAX_FILENAME, fp) != NULL) {
      if (sscanf(buf, "%d %d %f %s", &(rups[*numrup].src),
		 &(rups[*numrup].rup),
		 &(rups[*numrup].mag),
		 rups[*numrup].filename) == 4) {
	if (rups[*numrup].filename[strlen(rups[*numrup].filename) - 1] 
	    == '\n') {
	  rups[*numrup].filename[strlen(rups[*numrup].filename) - 1] = '\0';
	}
	if ((strlen(rups[*numrup].filename) == 0) || 
	    (rups[*numrup].src < 0) || 
	    (rups[*numrup].rup < 0) || 
	    (rups[*numrup].mag <= 0.0)) {
	  fprintf(stderr, "Invalid rupture file rec: %s\n", buf);
	  return(1);
	}
	rups[(*numrup)].stats.numslip = 0;
	rups[(*numrup)].stats.numhypo = 0;
	rups[(*numrup)].request = RUPGEN_STAGE_CALC;
	(*numrup)++;
      }
    }
  }
  fclose(fp); 

  /* Reassign indices */
  for (i = 0; i < *numrup; i++) {
    rups[i].index = i;
  }
  return(0);
}


/* Run Rob Graves GenSlip rupture generator */
int run_genslip(char *infile, char *outfile, char *logfile,
		rg_stats_t *stats, int doslip, int dohypo, 
		int writeout)
{
  int j, rgargc;
  char **rgargv = NULL;

  /* Create pseudo-program args */
  rgargc = NUM_GENSLIP_ARGS;
  rgargv = malloc(rgargc * sizeof(char*));
  for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
    rgargv[j] = malloc(MAX_STR_ARGV);
  }
  sprintf(rgargv[0], "%s", "rupgen");
  sprintf(rgargv[1], "infile=%s", infile);
  sprintf(rgargv[2], "outfile=%s", outfile);
  sprintf(rgargv[3], "logfile=%s", logfile);
  sprintf(rgargv[4], "doslip=%d", doslip);
  sprintf(rgargv[5], "dohypo=%d", dohypo);
  sprintf(rgargv[6], "writeout=%d", writeout);

  /* Run rupture generator */
  genslip(rgargc, rgargv, stats);
  
  /* Free pseudo-program args */
  for (j = 0; j < NUM_GENSLIP_ARGS; j++) {
    free(rgargv[j]);
  }

  free(rgargv);
  return(0);
}


/* Server loop */
int server(int myid, int nproc, rg_rfile_t *rups, int numrup, int s)
{
  MPI_Status status;
  int n;
  rg_rfile_t rup;
  int currup = 0;

  /* Dispatch rupture files to worker pool */
  memset(&rup, 0, sizeof(rg_rfile_t));
  currup = 0;
  while (currup < numrup) {
    /* Blocking receive on worker pool */
    MPI_Recv(&rup, 1, MPI_RUP_T, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &status);
    
    /* Update rupture record if necessary */
    if (strlen(rup.filename) > 0) {
      if ((rup.index < 0) || (rup.src < 0) || (rup.rup < 0) || 
	  (rup.stats.numslip <= 0) || (rup.stats.numhypo <= 0)) {
	fprintf(stderr, "[%d] Worker %d returned bad data for rupture %d\n", 
		myid, status.MPI_SOURCE, rup.index);
	return(1);
      }      
      memcpy(&(rups[rup.index]), &rup, sizeof(rg_rfile_t));
    }
    
    /* Send next available rupture to worker */
    memcpy(&rup, &(rups[currup]), sizeof(rg_rfile_t));
    rup.request = s;
    MPI_Send(&rup, 1, MPI_RUP_T, status.MPI_SOURCE, 
	     0, MPI_COMM_WORLD);
    
    if ((currup % 100 == 0) && (currup > 0)) {
      if (s == RUPGEN_STAGE_CALC) {
	printf("[%d] Dispatched %d rupture(s)\n", myid, currup);
      } else {
	printf("[%d] Dispatched %d slip(s)\n", myid, currup);
      }
    }
    
    currup++;  
  }
  
  /* Send stop message to all workers */
  printf("[%d] All ruptures dispatched, stopping workers\n", myid);
  for (n = 1; n < nproc; n++) {
    /* Blocking receive on worker pool */
    MPI_Recv(&rup, 1, MPI_RUP_T, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	     MPI_COMM_WORLD, &status);
    
    /* Update rupture record if necessary */
    if (strlen(rup.filename) > 0) {
      if ((rup.index < 0) || (rup.src < 0) || (rup.rup < 0) || 
	  (rup.stats.numslip <= 0) || (rup.stats.numhypo <= 0)) {
	fprintf(stderr, "[%d] Worker %d returned bad data for rupture %d\n", 
		myid, status.MPI_SOURCE, rup.index);
	return(1);
      }
      memcpy(&(rups[rup.index]), &rup, sizeof(rg_rfile_t));
    }
    
    /* Send stop work to worker */
    memset(&rup, 0, sizeof(rg_rfile_t));
    MPI_Send(&(rup), 1, MPI_RUP_T, status.MPI_SOURCE, 
	     0, MPI_COMM_WORLD);
    if ((n % 100 == 0) || ((nproc - n) < 20))  {
      printf("[%d] Stopped %d worker(s)\n", myid, n);
    }
  }
    
  return(0);
}


/* Worker loop */
int worker(int myid, int nproc, rg_conf_t *conf, int s)
{
  char infile[MAX_FILENAME];
  char outfile[MAX_FILENAME];
  char logfile[MAX_FILENAME];
  rg_rfile_t rup;
  int done = 0;

  memset(&rup, 0, sizeof(rg_rfile_t));
  done = 0;
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
      sprintf(infile, "%s%s", conf->rupdir, rup.filename);
      if (!file_exists(infile)) {
	fprintf(stderr, "[%d] Rupture file %s not found\n",
		myid, infile);
	return(1);
      }
      
      /* Create log directory */
      sprintf(logfile, "%s/%d", conf->logdir, rup.src);
      mkdir(logfile, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (s == RUPGEN_STAGE_WRITE) {
	sprintf(logfile, "%s/%d/rupgen_%d_%d_%d_%d.log", conf->logdir, 
		rup.src, rup.src, rup.rup, 
		rup.stats.numslip, rup.stats.numhypo);
      } else {
	sprintf(logfile, "%s/%d/rupgen_%d_%d.log", conf->logdir, 
		rup.src, rup.src, rup.rup);
      }

      /* Create slip directory */
      if (s == RUPGEN_STAGE_WRITE) {
	sprintf(outfile, "%s/%d/%d/%d", conf->rupdir, rup.src, rup.rup,
		rup.stats.numslip);
	mkdir(outfile, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	//sprintf(outfile, "%s/%d_%d.txt.variation", outfile,
	//	rup.src, rup.rup);
	sprintf(outfile, "%s/%s.variation", outfile,
		basename(rup.filename)); 
      } else {
	sprintf(outfile, "%s%s.variation", conf->rupdir, rup.filename);
      }

      /* Generate the rupture variations */
      if (s == RUPGEN_STAGE_CALC) {
	if (run_genslip(infile, outfile, logfile, &(rup.stats),
			GENSLIP_WILDCARD, GENSLIP_WILDCARD, 
			GENSLIP_NOWRITE) != 0) {
	  fprintf(stderr, "[%d] Failed to run rupture generator\n", myid);
	  return(1);
	}
      } else {
	if (run_genslip(infile, outfile, logfile, &(rup.stats),
			rup.stats.numslip, rup.stats.numhypo, 
			GENSLIP_WRITE) != 0) {
	  fprintf(stderr, "[%d] Failed to run rupture generator\n", myid);
	  return(1);
	}
      }
    }
  }
  
  return(0);
}


/* Parse config file */
int parseconfig(int myid, const char *path, rg_conf_t *conf)
{
  rg_config_t *rg_cfg = NULL;
  rg_config_t *cfgentry = NULL;

  if (myid == 0) {
    /* Parse conf */

    rg_cfg = rg_parse_config(path);
    if (rg_cfg == NULL) {
      fprintf(stderr, "[%d] Failed to read conf file %s\n", 
	      myid, path);
      return(1);
    }

    cfgentry = rg_find_name(rg_cfg, "rupdir");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key rupdir not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->rupdir, cfgentry->value, MAX_FILENAME);

    cfgentry = rg_find_name(rg_cfg, "logdir");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key logdir not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->logdir, cfgentry->value, MAX_FILENAME);

    cfgentry = rg_find_name(rg_cfg, "erfid");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key erfid not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->erfid, cfgentry->value, MAX_ID_LEN);

    cfgentry = rg_find_name(rg_cfg, "rvid");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key rvid not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->rvid, cfgentry->value, MAX_ID_LEN);

    cfgentry = rg_find_name(rg_cfg, "gridftp");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key gridftp not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->gridftp, cfgentry->value, MAX_FILENAME);

    cfgentry = rg_find_name(rg_cfg, "pool");
    if (cfgentry == NULL) {
      fprintf(stderr, "[%d] Key pool not found in %s\n", 
	      myid, path);
      return(1);
    }
    rg_strcpy(conf->pool, cfgentry->value, MAX_FILENAME);

    rg_free_config(rg_cfg);

    printf("[%d] rupdir=%s\n", myid, conf->rupdir);
    printf("[%d] logdir=%s\n", myid, conf->logdir);
    printf("[%d] erf=%s\n", myid, conf->erfid);
    printf("[%d] rv=%s\n", myid, conf->rvid);
    printf("[%d] gridftp server=%s\n", myid, conf->gridftp);
    printf("[%d] pool=%s\n", myid, conf->pool);
    fflush(stdout);
  }

  /* Broadcast config */
  if (MPI_Bcast(conf, 1, MPI_CONF_T, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
    fprintf(stderr, "[%d] Failed to broadcast config\n", myid);
    return(1);
  }

  return(0);
}


int main(int argc, char **argv)
{
  /* MPI stuff and distributed computation variables */
  int myid, nproc, pnlen;
  char procname[128];

  rg_conf_t conf;
  rg_rfile_t *rups = NULL;
  rg_rfile_t *rupvars = NULL;
  int i, j, n, numrup, numvar, numslip;
  FILE *sfp = NULL;
  char outfile[MAX_FILENAME];

  /* Init MPI */
  mpi_init(&argc, &argv, &nproc, &myid, procname, &pnlen);

  /* Register new MPI data types */
  if (mpi_register_rupinfo(&MPI_RUP_T, &num_fields_rup) != 0) {
      fprintf(stderr, "[%d] Failed to register rup data type\n", myid);
      return(1);
  }

  /* Register new MPI data types */
  if (mpi_register_conf(&MPI_CONF_T, &num_fields_conf) != 0) {
      fprintf(stderr, "[%d] Failed to register conf data type\n", myid);
      return(1);
  }

  memset(&conf, 0 , sizeof(rg_conf_t));
  numrup = 0;

  if (myid == 0) {
    if (argc != 2) {
      printf("usage: %s rupgen.conf\n", argv[0]);
      return(0);
    }

    if (parseconfig(myid, argv[1], &conf) != 0) {
      fprintf(stderr, "[%d] Failed to parse config %s\n", myid, argv[1]);
      return(1);
    }

    if (nproc < 2) {
      fprintf(stderr, "[%d] nproc must be at least 2\n", myid);
      return(1);
    }

    /* Create the log directory */
    mkdir(conf.logdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Get number of rupture files */
    if (getnumrups(conf.rupdir, &numrup) != 0) {
      fprintf(stderr, "[%d] Failed to count rupture files\n", myid);
      return(1);
    }
    if (numrup == 0) {
      fprintf(stderr, "[%d] No rupture files to process\n", myid);
      return(1);
    } else {
      printf("[%d] Retrieved %d rupture files\n", myid, numrup);
    }

    /* Allocate rupture buffer */
    rups = malloc(numrup * sizeof(rg_rfile_t));
    if (rups == NULL) {
      fprintf(stderr, "[%d] Failed to allocate rupture file buffer\n",
	      myid);
      return(1);
    }
    memset(rups, 0, sizeof(rg_rfile_t)*numrup);

    /* Get rupture files */
    if (getrups(conf.rupdir, rups, numrup, &n) != 0) {
      fprintf(stderr, "[%d] Failed to find rupture files\n", myid);
      return(1);
    }

    /* Calculate number of variations */
    if (server(myid, nproc, rups, numrup, RUPGEN_STAGE_CALC) != 0) {
      return(1);
    }

    /* Print statistics */
    sprintf(outfile, "%s/variations.list", conf.rupdir);
    sfp = fopen(outfile, "w");
    if (sfp == NULL) {
      fprintf(stderr, "[%d] Failed to open statistics file %s\n", 
	      myid, outfile);
      return(1);
    }
    numvar = 0;
    numslip = 0;
    for (n = 0; n < numrup; n++) {
      numslip = numslip + rups[n].stats.numslip;
      numvar += (rups[n].stats.numslip * rups[n].stats.numhypo);
      fprintf(sfp, "%d %d %d %d\n", rups[n].src, rups[n].rup,
	      rups[n].stats.numslip, rups[n].stats.numhypo);
    }
    
    fclose(sfp);  
    printf("[%d] Ruptures: %d\n", myid, numrup);
    printf("[%d] Rupture Slips: %d\n", myid, numslip);
    printf("[%d] Rupture Variations: %d\n", myid, numvar);
    fflush(stdout);
  
    /* Save RLS logical finame mappings */
    sprintf(outfile, "%s/rls.list", conf.rupdir);
    printf("[%d] Saving RLS filename mappings to %s\n", myid, outfile);
    sfp = fopen(outfile, "w");
    if (sfp == NULL) {
      fprintf(stderr, "[%d] Failed to open RLS mapping file %s\n", 
	      myid, outfile);
      return(1);
    }
    for (n = 0; n < numrup; n++) {
      for (i = 0; i < rups[n].stats.numslip; i++) {
	for (j = 0; j < rups[n].stats.numhypo; j++) {
	  fprintf(sfp, "e%s_rv%s_%s.variation-s%04d-h%04d, %s%s/%d/%d/%d/%s.variation-s%04d-h%04d, %s\n", 
		  conf.erfid, conf.rvid, basename(rups[n].filename), 
		  i, j, conf.gridftp, conf.rupdir, rups[n].src, rups[n].rup,
		  i, basename(rups[n].filename), i, j, conf.pool);
	}
      }
    }
    fclose(sfp);

    mpi_barrier();

    /* Generate all rupture variations by slip */
    rupvars = malloc(numslip * sizeof(rg_rfile_t));
    if (rupvars == NULL) {
      fprintf(stderr, "[%d] Failed to allocate rupture var file buffer\n",
	      myid);
      return(1);
    }
    memset(rupvars, 0, sizeof(rg_rfile_t)*numslip);

    numvar = 0;
    for (n = 0; n < numrup; n++) {
      for (i = 0; i < rups[n].stats.numslip; i++) {
	memcpy(&(rupvars[numvar]), &(rups[n]), sizeof(rg_rfile_t));
	rupvars[numvar].request = RUPGEN_STAGE_WRITE;
	rupvars[numvar].stats.numslip = i;
	rupvars[numvar].stats.numhypo = GENSLIP_WILDCARD; 
	numvar++;
      }
    }

    /* Sort rupture slips */
    qsort(rupvars, numvar, sizeof(rg_rfile_t), rupcomp);

    /* Reassign indices after sort */
    for (n = 0; n < numvar; n++) {
      rupvars[n].index = n;
    }

    /* Dispatch slips to workers */
    if (server(myid, nproc, rupvars, numvar, RUPGEN_STAGE_WRITE) != 0) {
      return(1);
    }

    printf("[%d] Complete\n", myid);
    fflush(stdout);

    free(rups);
    free(rupvars);
  } else {

    if (parseconfig(myid, argv[1], &conf) != 0) {
      fprintf(stderr, "[%d] Failed to parse config %s\n", myid, argv[1]);
      return(1);
    }

    if (worker(myid, nproc, &conf, RUPGEN_STAGE_CALC) != 0) {
      return(1);
    }
    mpi_barrier();
    if (worker(myid, nproc, &conf, RUPGEN_STAGE_WRITE) != 0) {
      return(1);
    }
  }

  /* Final sync */
  mpi_barrier();
  mpi_final("MPI Done");

  return(0);
}
