/*
 ============================================================================
 Name        : DirectSynth.c
 Author      : Scott Callaghan
 Version     :
 Description : Performs CyberShake post-processing without writing SGT extraction data
 ============================================================================
 */

#include "include.h"
#include "defs.h"
#include "structure.h"
#include "duration.h"
#include "functions.h"

int debug = 0;
int my_global_id = -1;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	int my_id, num_procs;

	//Number of processes responsible for reading the SGTs, including rank 0, but it just coordinates
	int num_sgt_handlers = 0;

	getMPIInfo(&my_id, &num_procs);
	my_global_id = my_id;

	setpar(argc, argv);
	mstpar("sgt_handlers", "d", &num_sgt_handlers);

	//station parameters
	float slat, slon;
	float sdep = 0.0;
	char stat[64];
	mstpar("slat", "f", &slat);
	mstpar("slon", "f", &slon);
	mstpar("stat","s",stat);

	//Input files
	struct sgtfileparams sgtfilepar;
	init_sgtfileparams(&sgtfilepar);
	int num_comps = 0;
	num_comps += getpar("sgt_xfile","s",sgtfilepar.xfile);
	num_comps += getpar("sgt_yfile","s",sgtfilepar.yfile);
	num_comps += getpar("sgt_zfile","s",sgtfilepar.zfile);

	if(num_comps==0) {
		if (my_id==0) {
			fprintf(stderr,"*** need to specify at least one of sgt_xfile, sgt_yfile, or sgt_zfile; exiting ...\n");
		}
		MPI_Finalize();
		exit(-1);
	}

	if (sgtfilepar.xfile[0] != '\0' && sgtfilepar.zfile[0] != '\0') {
		if (my_id==0) {
			fprintf(stderr, "We don't support only X and Z components yet. Exiting.\n");
		}
		MPI_Finalize();
		exit(-2);
	}

	getpar("x_header","s",sgtfilepar.xfile_header);
	getpar("y_header","s",sgtfilepar.yfile_header);
	getpar("z_header","s",sgtfilepar.zfile_header);

	//Needed for output filenames
	int run_id;
	mstpar("run_id","d",&run_id);

	int run_PSA, run_rotd, run_duration;
	//Default is to run everything
	run_PSA = 1;
	run_rotd = 1;
	run_duration = 1;
	getpar("run_psa","d",&run_PSA);
	getpar("run_rotd","d",&run_rotd);
	getpar("run_duration","d",&run_duration);

	//track timing
	int timing = 0;
	struct timeval tv_start, tv_end;
	getpar("timing","d",&timing);


	getpar("debug","d",&debug);
	if (debug) {
		open_log(my_id);
	}

	//Directory hierarchy at the source level
        int src_dir_hierarchy = 0;
        getpar("dir_hierarchy","d",&src_dir_hierarchy);

	//Determine core, hostname
	int core = sched_getcpu();
	char hostname[256];
	gethostname(hostname, 256);
	if (debug) {
		char buf[512];
		sprintf(buf, "Running on core %d of host %s", core, hostname);
		write_log(buf);
	}

	MPI_Comm sgt_handler_comm, sgt_readers_comm;
	//Includes ranks 0 through num_sgt_handlers - 1
        constructSGTHandlerComm(num_sgt_handlers, &sgt_handler_comm);
	//Includes ranks 1 through num_sgt_handlers - 1, for reading in SGT files
	constructSGTReaderComm(num_sgt_handlers, &sgt_readers_comm);

	if (my_id<num_sgt_handlers) {
		if (my_id==0) {
			if (debug) write_log("Entering master.");
			master(&sgtfilepar, &sgt_handler_comm, num_sgt_handlers, stat, run_id, run_PSA, run_rotd, run_duration, src_dir_hierarchy);
		} else if (my_id<num_sgt_handlers) {
			if (debug) write_log("Entering sgt handler.");
			sgt_handler(&sgtfilepar, num_comps, &sgt_handler_comm, num_sgt_handlers, &sgt_readers_comm, my_id);
		}
	} else {
		//Frequency information
		float det_max_freq = 0.5;
		float stoch_max_freq = -1.0;
		getpar("det_max_freq", "f", &det_max_freq);
		getpar("stoch_max_freq", "f", &stoch_max_freq);

		//Runtime info
		int max_buf_mb = 1024;
		getpar("max_buf_mb","d",&max_buf_mb);
		long long MAX_BUFFER_SIZE = 1024*(long long)1024*max_buf_mb;

		int num_workers = num_procs - num_sgt_handlers - 1;

		if (my_id==num_sgt_handlers) {
			if (debug) write_log("Entering task manager.");
			task_manager(num_sgt_handlers, num_workers, num_procs, MAX_BUFFER_SIZE, my_id);
		} else if (my_id<num_procs){
			//Output options
			int ntout;

			mstpar("ntout","d",&ntout);

			if (debug) write_log("Entering worker.");
			worker(argc, argv, num_sgt_handlers, &sgtfilepar, stat, slat, slon, run_id, det_max_freq, stoch_max_freq, run_PSA, run_rotd, run_duration, my_id);
		}
	}

	if (debug) {
		close_log();
	}
	MPI_Finalize();
	return 0;
}
