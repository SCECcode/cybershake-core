/*
 * defs.h
 *
 *  Created on: Dec 12, 2014
 *      Author: scott
 */
#ifndef DEFS_H
#define DEFS_H

//Debug flag
extern int debug;
extern int my_global_id;

//component flags
extern const int X_COMP_FLAG;
extern const int Y_COMP_FLAG;
extern const int Z_COMP_FLAG;

//MPI message tags
#define NUM_RUPTURES_TAG 3
#define VARIATION_INFO_TAG 5
#define SGT_HEADER_TAG 10
#define DATA_FILENAME_TAG 15
#define DATA_TAG 20
#define SGT_REQUEST_TAG 30
#define SGT_HEADER_DATA_TAG 40
#define SGT_RAW_DATA_TAG 50
#define WORK_REQUEST_TAG 60
#define WORK_RESPONSE_TAG 70
#define WORK_COMPLETE_TAG 80

//Message types to master
#define WORKERS_COMPLETE 1
#define DATA_FILE 2

//Message types to handler
#define SGT_REQUEST 3

//Message types to task manager
#define WORK_REQUEST 4

//Message types to workers
#define WORK_RESPONSE 5
#define WORK_COMPLETE 6

//Worker status
#define WORKER_IDLE 0
#define WORKER_WORKING 1
#define WORKER_COMPLETE 2

//For SGT MPI I/O
#define N_SGTvars (6)

//Number of periods in SCEC array in surfseis_rspectra code
#define NUM_SCEC_PERIODS 44
//MAXPERIODS as defined in surfseis_rspectra
#define MAXPERIODS 305

//num periods as defined in rotd.c
#define NUM_ROTD_PERIODS 22

//Checkpoint filename
#define CHECKPOINT_FILE "checkpoint.out"

//Dimensions of lookup table with number of files per rupture
#define MAX_SOURCE_ID 11000
#define MAX_RUPTURE_ID 1300

//Defs from jbsim3d code
#define         NTMAX              10000
#define         SOMERVILLE_FLAG    1
#define         MAI_FLAG           2
#define         MINSLIP            1.0e-02

#define   DHYPO_FRAC       0.75     /* hypo at 0.75 down-dip width */
#define   SHYPO_STEP       20.0     /* hypo spacing at 20 km along strike */
#define   SHYPO_MIN_OFF    1.0      /* hypos start at 1.0 km along strike */
#define   SLIPS_TO_HYPOS   2    /* no. slip models = 2 times no. of hypos */

#define   U1FLAG   1
#define   U2FLAG   2
#define   U3FLAG   3

#if _FILE_OFFSET_BITS == 64

#define RDONLY_FLAGS    O_RDONLY | O_LARGEFILE
#define RDWR_FLAGS      O_RDWR | O_LARGEFILE
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR | O_LARGEFILE

#else

#define RDONLY_FLAGS    O_RDONLY
#define RDWR_FLAGS      O_RDWR
#define CROPTR_FLAGS    O_CREAT | O_TRUNC | O_RDWR

#endif

//integ_diff constants
#define FD_STATCHAR 8
#define STATCHAR 12
#define COMPCHAR 4
#define TITLCHAR 64

#define mmv 30000

#define SWAP_FLAG -12345
#define MAX_VAR_LIST 100
#define TAP_PERC 0.05

#endif
