/*
 * task_manager.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "defs.h"
#include "structure.h"
#include "functions.h"

/*
 * if my_id==N:
 *  Obtain sgtmast, sgtindx
	Obtain process/point mapping
	Read in list of ruptures
	Construct lists of tasks from ruptures (might not send all rupture variations to a core due to memory limitations)
	Service requests for next task
	When no more tasks, send messages to processors that request
	When sent no more task messages to all processors, send message to processors 0 to N-1
	exit
 */

void broadcast_completion(int num_sgt_handlers, int num_workers, int num_procs, int my_id);
void manager_listen(int num_workers, worker_task* task_list, int num_tasks, int my_id);
int handle_work_request(manager_msg msg, worker_task* task_list, int* current_task, int num_tasks, MPI_Datatype worker_msg_type, int my_id);
int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int rup_var_spacing, int my_id);
void get_point_mapping(int num_sgt_readers);
int parse_rv_info(char* rv_info_file, rv_info** rvinfo_array);
int cmp_rvinfo(const void* a, const void* b);

int task_manager(int num_sgt_handlers, int num_workers, int num_procs, long long MAX_BUFFER_SIZE, int rup_var_spacing, MPI_Comm* manager_plus_workers_comm, int my_id) {
	//Construct MPI datatype
	MPI_Datatype sgtmast_type, sgtindx_type;
	construct_sgtmast_datatype(&sgtmast_type);
	construct_sgtindx_datatype(&sgtindx_type);

	//Need dtout to calculate rupture variation size requirements
	float dtout = 0.1;
	getpar("dtout", "f", &dtout);

	//Get sgtmast, sgtindx info
	struct sgtmaster sgtmast;
	struct sgtindex* sgtindx;
	if (debug) write_log("Receiving sgtmast from master.");
	check_bcast(&sgtmast, 1, sgtmast_type, 0, MPI_COMM_WORLD, "Error receiving sgtmast, aborting.", my_id);
	sgtindx = check_malloc(sizeof(struct sgtindex)*sgtmast.globnp);
	if (debug) write_log("Receiving sgtindx from master.");
	check_bcast(sgtindx, sgtmast.globnp, sgtindx_type, 0, MPI_COMM_WORLD, "Error receiving sgtindx, aborting.", my_id);

	get_point_mapping(num_sgt_handlers);

    //See if rvfrac and the rv seed are provided
    char rv_info_file[256];
    rv_info_file[0] = '\0';
    getpar("rv_info_file","s",rv_info_file);
    rv_info* rvinfo_array = NULL;
    if (rv_info_file[0]!='\0') {
        if (debug) write_log("Parsing rvinfo file.");
        int num_rv_infos = parse_rv_info(rv_info_file, &rvinfo_array);
        //Broadcast to workers
        MPI_Datatype rvinfo_type;
        construct_rvinfo_datatype(&rvinfo_type);
        check_bcast(&num_rv_infos, 1, MPI_INT, 0, *manager_plus_workers_comm, "Error sending length of rvinfo data, aborting.", my_id);
        check_bcast(rvinfo_array, num_rv_infos, rvinfo_type, 0, *manager_plus_workers_comm, "Error sending rvinfo data, aborting.", my_id);
    }

	//Get rupture list
	char rup_list_file[256];
	worker_task* task_list;
	mstpar("rup_list_file","s",rup_list_file);

	int num_tasks = parse_rupture_list(rup_list_file, &task_list, MAX_BUFFER_SIZE, dtout, sgtmast.nt, sgtmast.globnp, rup_var_spacing, my_id);

	manager_listen(num_workers, task_list, num_tasks, my_id);

	//Sleep for 5 seconds to allow any final files going to the master to get accepted and written
	sleep(5);

	broadcast_completion(num_sgt_handlers, num_workers, num_procs, my_id);

	if (debug) write_log("Shutting down.");

	free(sgtindx);
	free(task_list);
	return 0;
}

void broadcast_completion(int num_sgt_handlers, int num_workers, int num_procs, int my_id) {
	//Send out completed message to master + SGT handlers
	//Need to send point-to-point messages since handlers are just listening
	if (debug) write_log("Notifying master + SGT handlers that all tasks are complete.");
	//We are constructing a handler message, but the master message format is identical
	handler_msg w_msg;
	w_msg.msg_type = WORKERS_COMPLETE;
	w_msg.msg_src = my_id;
	int i;
	for (i=0; i<num_sgt_handlers; i++) {
		if (debug) {
			char buf[256];
			sprintf(buf, "Sending complete message to process %d.", i);
			write_log(buf);
		}
		check_send(&w_msg, 3, MPI_INT, i, WORK_COMPLETE_TAG, MPI_COMM_WORLD, "Error sending work complete message, aborting.", my_id);
	}
}

void manager_listen(int num_workers, worker_task* task_list, int num_tasks, int my_id) {
	MPI_Datatype worker_msg_type;
	construct_worker_message_datatype(&worker_msg_type);

	manager_msg msg;
	int current_task = 0;
	//Keep track of the workers and who is working and who is complete
	int workers_working = 0;
	int* worker_status = check_malloc(sizeof(int)*num_workers);
	int i;
	for (i=0; i<num_workers; i++) {
		worker_status[i] = WORKER_IDLE;
	}
	if (num_tasks==0) {
		//Need to tell everyone that there's no more work
		int num_to_notify = num_workers;
		while(num_to_notify>0) {
			check_recv(&msg, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, "Error listening for message, aborting.", my_id);
			if (msg.msg_type==WORK_REQUEST) {
	                        if (debug) {
	                                char buf[256];
	                                sprintf(buf, "Received request for work from process %d.", msg.msg_src);
	                                write_log(buf);
	                        }
	                        handle_work_request(msg, task_list, &current_task, num_tasks, worker_msg_type, my_id);
				num_to_notify--;
			}
		}
		return;
	}

	while (1) {
		if (debug) write_log("Waiting for messages.");
		check_recv(&msg, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, "Error listening for message, aborting.", my_id);
		if (msg.msg_type==WORK_REQUEST) {
			if (debug) {
				char buf[256];
				sprintf(buf, "Received request for work from process %d.", msg.msg_src);
				write_log(buf);
			}
			int completed_worker = handle_work_request(msg, task_list, &current_task, num_tasks, worker_msg_type, my_id);
			//Convert the message source into a worker offset for worker_status
			int worker_offset = msg.msg_src - my_id - 1;

			if (completed_worker==-1) {
				if (worker_status[worker_offset]!=WORKER_WORKING) {
					worker_status[worker_offset] = WORKER_WORKING;
					workers_working++;
				}
			} else {
				//This means work is completing
				if (worker_status[worker_offset]==WORKER_WORKING) {
					workers_working--;
				}
				worker_status[worker_offset] = WORKER_COMPLETE;
			}
		} else {
			fprintf(stderr, "%d) Task manager received message of unknown type %d from %d, aborting.", my_id, msg.msg_type, msg.msg_src);
			if (debug) close_log();
			MPI_Finalize();
			exit(4);
		}
		if (debug) {
			char buf[256];
			sprintf(buf, "%d workers are working on tasks.", workers_working);
			write_log(buf);
		}
		if (workers_working==0) {
			if (debug) write_log("All workers complete.");
			break;
		}
	}
	free(worker_status);
	MPI_Type_free(&worker_msg_type);
}

int handle_work_request(manager_msg msg, worker_task* task_list, int* current_task, int num_tasks, MPI_Datatype worker_msg_type, int my_id) {
	//Send task back if there are any
	worker_msg out_msg;
	if (*current_task<num_tasks) {
		if (debug) {
			char buf[256];
			sprintf(buf, "Assigning task %d to process %d.", *current_task, msg.msg_src);
			write_log(buf);
		}
		worker_msg out_msg;
		out_msg.msg_type = WORK_RESPONSE;
		out_msg.task = task_list[(*current_task)];
		check_send(&out_msg, 1, worker_msg_type, msg.msg_src, WORK_RESPONSE_TAG, MPI_COMM_WORLD, "Error sending work message to worker, aborting.", my_id);
		(*current_task)++;
		return -1;
	} else {
		//No more work; notify worker
		if (debug) {
			char buf[256];
			sprintf(buf, "Notifying process %d there is no more work.", msg.msg_src);
			write_log(buf);
		}
		out_msg.msg_type = WORK_COMPLETE;
		check_send(&out_msg, 1, worker_msg_type, msg.msg_src, WORK_COMPLETE_TAG, MPI_COMM_WORLD, "Error sending work complete message to worker, aborting.", my_id);
		return msg.msg_src;
	}
}

int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int rup_var_spacing, int my_id) {
	int i, j;
	/*File has format
	 * <rupture file> <# of slips> <# of hypos> <num_points>
	 */
	if (debug) write_log("Parsing rupture list to construct task list.");
	FILE* rup_list_in;
	int num_ruptures;
	fopfile_ro(rup_list_file, &rup_list_in);
	//first line contains number of ruptures
	fscanf(rup_list_in, "%d", &num_ruptures);
	if (num_ruptures==0) {
		fprintf(stderr, "First line of rupture list file must contain number of ruptures, aborting.\n");
		if (debug) close_log();
		MPI_Finalize();
		exit(5);
	}
	//Build (source id, rupture id, num rup vars) tuples to send to master, for checkpointing
	int* tuples_for_master = check_malloc(sizeof(int)*num_ruptures*3);

	//Check for checkpoint file
	char** completed_list = NULL;
	int num_completed = 0;
	if (access(CHECKPOINT_FILE, R_OK)==0) {
		completed_list = check_malloc(sizeof(char*)*num_ruptures);
		FILE* checkpoint_fp;
		fopfile_ro(CHECKPOINT_FILE, &checkpoint_fp);
		char line[16];
		char* ret = fgets(line, 16, checkpoint_fp);
		while (ret!=NULL) {
			char* tok;
			tok = strtok(line, "\n");
			completed_list[num_completed] = check_malloc(sizeof(char)*16);
			strcpy(completed_list[num_completed], tok);
			num_completed++;
			ret = fgets(line, 16, checkpoint_fp);
		}
		fclose(checkpoint_fp);
	}
	if (debug) {
		char buf[256];
		sprintf(buf, "%d ruptures already completed in checkpoint file.", num_completed);
		write_log(buf);
	}

	//Start here, but will need to grow
	int task_list_length = num_ruptures;
	(*task_list) = check_malloc(sizeof(worker_task)*task_list_length);
	char rupture_file[256], string[256];
	int num_slips, num_hypos, num_rows, num_cols, num_points, num_tasks;
	float mag;
	//num ruptures to process is different, because some might have already completed
	int num_ruptures_to_process = 0;
	num_tasks = 0;

	#ifdef _V3_3_1
	//Define constants for estimating memory
	//a = 32.63, b = -81.26, c = 104.3, d = 54.6      
        //Apply scale factor of .6*dt^-.17
	float CONST_A = 32.63;
	float CONST_B = -81.26;
	float CONST_C = 104.3;
	float CONST_D = 54.6;
	float DT_SCALE_FAC = 0.6*pow(dtout, -0.17);
	#endif

	for (i=0; i<num_ruptures; i++) {
		fscanf(rup_list_in, "%s %d %d %d %d %f", rupture_file, &num_slips, &num_hypos, &num_rows, &num_cols, &mag);
		num_points = num_rows * num_cols;

		//Determine source and rupture ID
		//Files in format e<ERF_ID>_rv<RV_ID>_<source>_<rupture>.txt
		//Could also be just <source>_<rupture>.txt
		//Need 2nd to last, 3rd to last
		char* tok, *last, *ntl, *ttl;
		tok = last = ntl = ttl = 0;
		strcpy(string, rupture_file);
		tok = strtok(string, "/_.");
		while (tok!=NULL) {
			ttl = ntl;
			ntl = last;
			last = tok;
			tok = strtok(NULL, "/_.");
		}
		int source_id = atoi(ttl);
		int rupture_id = atoi(ntl);

		//If we have a checkpoint file, see if we've done this one already
		if (num_completed>0) {
			char search_string[16];
			sprintf(search_string, "%d %d", source_id, rupture_id);
			int flag = 0;
			for (j=0; j<num_completed; j++) {
				if (strcmp(search_string, completed_list[j])==0) {
					//We found it
					flag = 1;
					break;
				}
			}
			if (flag==1) {
				if (debug) {
					char buf[256];
					sprintf(buf, "Skipping source %d, rupture %d.", source_id, rupture_id);
					write_log(buf);
				}
				//Move to next entry in rupture list
				continue;
			}
		}

		//Want to make sure we're not exceeding 1.8 GB per task
		long long full_sgt_data = (long long)num_points*(3*sizeof(float)*N_SGTvars*nt + sizeof(struct sgtheader));
		long long sgt_size = full_sgt_data;
		if (sgt_size>MAX_BUFFER_SIZE) {
			sgt_size = MAX_BUFFER_SIZE;
		}
		//printf("sgt_size = %ld\n", sgt_size);
		//printf("num_points = %d\n", num_points);
		//printf("dtout = %f\n", dtout);
		#ifdef _V3_3_1
                //Each rupture variation adds roughly
                //14.8 * log10(rupture_points) * rupture_points^1.14 MB worth of storage * 1.1(tolerance) * 0.1/dtout
		/*long long single_rv_size = (long long)(14.8 * (long long)log10(num_points) * pow(num_points, 1.14) * 1.1);
                //Cap rupture size at 80 MB
                if (single_rv_size > 80*1024*1024) {
                        single_rv_size = 80*1024*1024;
                }
                //Take dtout into consideration
                single_rv_size = (long long)((float)single_rv_size * 0.1/dtout);
		*/
		
		/*New formula, taking magnitude into account:
		For dt = 0.1
		shifted mag = mag - 5.35
		mem per point = ax^3 + bx^2 + cx + d
		x = shifted mag, a = 32.63, b = -81.26, c = 104.3, d = 54.6
		
		Apply scale factor of .6*dt^-.17
		Multiply by # of points
		*/
		float shifted_mag = mag - 5.35;
		float mem_per_pt = CONST_A*shifted_mag*shifted_mag*shifted_mag + CONST_B*shifted_mag*shifted_mag + CONST_C*shifted_mag + CONST_D;
		long long single_rv_size = (long long)(mem_per_pt * num_points * DT_SCALE_FAC * 1.1);
		//Need to include for compatibility
		long long generation_size = 0;
		printf("shifted_mag=%f, mem_per_pt=%f, single RV size for src %d rup %d is %ld.\n", shifted_mag, mem_per_pt, source_id, rupture_id, single_rv_size);
		#else
		//For v5.2.3; assumes dtout=0.05.  size = 6.491 * num_points ^ 1.314 * 1.1 (tolerance)
		long long single_rv_size = (long long)(6.491*pow(num_points, 1.314)*1.1);
		//Cap at 120 MB
		if (single_rv_size > 120*1024*1024) {
			single_rv_size = 120*1024*1024;
		}
		//Also need to include space for rough_c, tsfac2_c, rtime2_c arrays, used when generating ruptures
        //Each of these has size sizeof(struct complex)*(3.3*nstk)^2, plus 10% slop
		int size_struct_complex = 8;
		int nstk = num_cols;
		if (num_rows>num_cols) {
			nstk = num_rows;
		}
        long long generation_size = 3*size_struct_complex*10.89*nstk*nstk*1.1;
		#endif
		//Do not permit more than this amount to be used
		//Usage is rupture variations, SGTs, memcached (32 MB/core), and gfmech; everything else is tiny
		//Summit has 3.04 GB per hardware thread, or 12.16 GB per core
		long long MAX_ALLOWED = (long long)(1.5 * 1024.0 * 1024.0 * 1024.0);
		int memcached = 32*1024*1024;
		//gfmech uses 1st power of 2 larger than 4*nt
		int power = ceil(log2(nt))+2;
		int gfmech_size = (int)(3.0*12.0*sizeof(float)*pow(2, power));
		//printf("MAX_ALLOWED = %ld\n", MAX_ALLOWED);
		int num_vars_per_task = (MAX_ALLOWED - sgt_size - memcached - generation_size)/(single_rv_size + gfmech_size);
		if (num_vars_per_task<1) {
			if (debug) {
				char buf[256];
				sprintf(buf, "num vars per task is %d, adjusting to minimum of 1.", num_vars_per_task);
				write_log(buf);
			}
			num_vars_per_task = 1;
		}
		//int num_vars_per_task = (MAX_ALLOWED - sgt_size)/single_rv_size;
		//Change for debugging
		//num_vars_per_task = 2;
		printf("MAX_ALLOWED=%ld, sgt_size=%ld, single_rv_size=%ld\n", MAX_ALLOWED, sgt_size, single_rv_size);
		printf("num_vars_per_task = %d\n", num_vars_per_task);
		int num_vars = num_slips * num_hypos;
		int tasks_for_rupture = ceil(((float)num_vars)/((float)num_vars_per_task));
		if (debug) {
			char buf[256];
			sprintf(buf, "rupture file %s (source_id %d, rupture_id %d) has %d vars, which yield %d tasks.", rupture_file, source_id, rupture_id, num_vars, tasks_for_rupture);
			write_log(buf);
		}

		tuples_for_master[3*num_ruptures_to_process] = source_id;
		tuples_for_master[3*num_ruptures_to_process+1] = rupture_id;
		tuples_for_master[3*num_ruptures_to_process+2] = num_vars;

		num_ruptures_to_process++;

		for (j=0; j<tasks_for_rupture; j++) {
			strcpy((*task_list)[num_tasks].rupture_filename, rupture_file);
			//printf("rupture_filename = %s\n", (*task_list)[num_tasks].rupture_filename);
			(*task_list)[num_tasks].source_id = source_id;
			(*task_list)[num_tasks].rupture_id = rupture_id;
			(*task_list)[num_tasks].num_slips = num_slips;
			(*task_list)[num_tasks].num_hypos = num_hypos;
			(*task_list)[num_tasks].rup_var_spacing = rup_var_spacing;
			(*task_list)[num_tasks].starting_rv_id = j*num_vars_per_task;
			int ending_rv_id = (j+1)*num_vars_per_task;
			if (ending_rv_id > num_vars) {
				ending_rv_id = num_vars;
			}
			(*task_list)[num_tasks].ending_rv_id = ending_rv_id;
			num_tasks++;
			if (num_tasks>=task_list_length) {
				//Needs to be at least enough longer to handle the rest of the tasks
				//At least the rest of the tasks for this rupture + the # of rupture left
				task_list_length += tasks_for_rupture-j + num_ruptures-i - 1;
				if (debug) {
					char buf[256];
					sprintf(buf, "Increasing task list to length %d", task_list_length);
					write_log(buf);
				}
				(*task_list) = check_realloc(*task_list, sizeof(worker_task)*task_list_length);
			}
		}
	}
	if (debug) {
		char buf[256];
		sprintf(buf, "%d tasks constructed from %d ruptures.", num_tasks, num_ruptures);
		write_log(buf);
	}

	for (i=0; i<num_completed; i++) {
		free(completed_list[i]);
	}
	free(completed_list);

	//Send info to master, for checkpointing
	if (debug) write_log("Sending task info to master.");
	check_send(&num_ruptures_to_process, 1, MPI_INT, 0, NUM_RUPTURES_TAG, MPI_COMM_WORLD, "Error sending num ruptures to master, aborting.", my_id);
	check_send(tuples_for_master, 3*num_ruptures_to_process, MPI_INT, 0, VARIATION_INFO_TAG, MPI_COMM_WORLD, "Error sending (source, rupture, #rvs) tuples to master, aborting.", my_id);

	free(tuples_for_master);
	return num_tasks;
}

void get_point_mapping(int num_sgt_readers) {
	//Even though we don't need it, the point-to-process mapping is broadcast to all
	int* proc_points = check_malloc(sizeof(int)*(num_sgt_readers+1));
	if (debug) write_log("Receiving SGT point-to-process mapping from master.");
	check_bcast(proc_points, num_sgt_readers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error broadcasting SGT point-to-process mapping to all, aborting.", 0);

	free(proc_points);
}

int parse_rv_info(char* rv_info_file, rv_info** rvinfo_array) {
    int num_rv_infos = 0;
    int src_id, rup_id, rv_id, seed;
    float rvfrac;
    FILE* rv_info_fp_in;
    fopfile_ro(rv_info_file, &rv_info_fp_in);
    //File has format
    //Number of entries
    //<src id> <rup id> <rv_id> <rvfrac> <seed>
    //If rvfrac or seed shouldn't be used, they're set to -1 in the file
    fscanf(rv_info_fp_in, "%d", &num_rv_infos);
	(*rvinfo_array) = check_malloc(sizeof(rv_info)*num_rv_infos);
	int i;
	for (i=0; i<num_rv_infos; i++) {
    	fscanf(rv_info_fp_in, "%d %d %d %f %d", &src_id, &rup_id, &rv_id, &rvfrac, &seed);
        (*rvinfo_array)[i].source_id = src_id;
        (*rvinfo_array)[i].rupture_id = rup_id;
        (*rvinfo_array)[i].rup_var_id = rv_id;
        (*rvinfo_array)[i].rvfrac = rvfrac;
        (*rvinfo_array)[i].seed = seed;
    }
    fclose(rv_info_fp_in);
    //Sort for easy searching
    qsort(&((*rvinfo_array)[0]), num_rv_infos, sizeof(rv_info), cmp_rvinfo);
    return num_rv_infos;
}

int cmp_rvinfo(const void* a, const void* b) {
	rv_info* rvi_a = (rv_info*)a;
	rv_info* rvi_b = (rv_info*)b;
	if (rvi_a->source_id != rvi_b->source_id) {
		return (rvi_a->source_id - rvi_b->source_id);
	} else {
		if (rvi_a->rupture_id != rvi_b->rupture_id) {
			return (rvi_a->rupture_id - rvi_b->rupture_id);
		} else {
			return (rvi_a->rup_var_id - rvi_b->rup_var_id);
		}
	}
}
