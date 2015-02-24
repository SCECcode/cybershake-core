/*
 * task_manager.c
 *
 *  Created on: Dec 11, 2014
 *      Author: scott
 */

#include "include.h"
#include "structure.h"
#include "functions.h"
#include "defs.h"

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
int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int rup_var_spacing);
void get_point_mapping(int num_sgt_readers);

int task_manager(int num_sgt_handlers, int num_workers, int num_procs, long long MAX_BUFFER_SIZE, int rup_var_spacing, int my_id) {
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

	//Get rupture list
	char rup_list_file[256];
	worker_task* task_list;
	mstpar("rup_list_file","s",rup_list_file);

	int num_tasks = parse_rupture_list(rup_list_file, &task_list, MAX_BUFFER_SIZE, dtout, sgtmast.nt, sgtmast.globnp, rup_var_spacing);

	manager_listen(num_workers, task_list, num_tasks, my_id);

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
			int worker_offset = msg.msg_src - my_id;

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

int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int rup_var_spacing) {
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
	//Start here, but will need to grow
	int task_list_length = num_ruptures;
	(*task_list) = check_malloc(sizeof(worker_task)*task_list_length);
	int i, j;
	char rupture_file[256], string[256];
	int num_slips, num_hypos, num_points, num_tasks;
	num_tasks = 0;
	for (i=0; i<num_ruptures; i++) {
		fscanf(rup_list_in, "%s %d %d %d", rupture_file, &num_slips, &num_hypos, &num_points);

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

		//Want to make sure we're not exceeding 1.8 GB per task
		long long full_sgt_data = (long long)num_points*(3*sizeof(float)*N_SGTvars*nt + sizeof(struct sgtheader));
		long long sgt_size = full_sgt_data;
		if (sgt_size>MAX_BUFFER_SIZE) {
			sgt_size = MAX_BUFFER_SIZE;
		}
		printf("sgt_size = %ld\n", sgt_size);
		printf("num_points = %d\n", num_points);
		printf("dtout = %f\n", dtout);
		//Each rupture variation adds roughly
		//14.8 * log10(rupture_points) * rupture_points^1.14 MB worth of storage * 1.1(tolerance) * 0.1/dtout
		long long single_rv_size = (long long)(14.8 * (long long)log10(num_points) * pow(num_points, 1.14) * 1.1);
		//Cap rupture size at 80 MB
		if (single_rv_size > 80*1024*1024) {
			single_rv_size = 80*1024*1024;
		}
		//Take dtout into consideration
		single_rv_size = (long long)((float)single_rv_size * 0.1/dtout);
		printf("single_rv_size = %ld\n", single_rv_size);
		//Do not permit more than this amount to be used
		long long MAX_ALLOWED = (long long)(1.8 * 1024.0 * 1024.0 * 1024.0);
		printf("MAX_ALLOWED = %ld\n", MAX_ALLOWED);
		int num_vars_per_task = (MAX_ALLOWED - sgt_size)/single_rv_size;
		//Change for debugging
		//num_vars_per_task = 2;
		printf("num_vars_per_task = %d\n", num_vars_per_task);
		int num_vars = num_slips * num_hypos;
		int tasks_for_rupture = ceil(((float)num_vars)/((float)num_vars_per_task));
		if (debug) {
			char buf[256];
			sprintf(buf, "rupture file %s (source_id %d, rupture_id %d) has %d vars, which yield %d tasks.", rupture_file, source_id, rupture_id, num_vars, tasks_for_rupture);
			write_log(buf);
		}

		for (j=0; j<tasks_for_rupture; j++) {
			strcpy((*task_list)[num_tasks].rupture_filename, rupture_file);
			printf("rupture_filename = %s\n", (*task_list)[num_tasks].rupture_filename);
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
	return num_tasks;
}

void get_point_mapping(int num_sgt_readers) {
	//Even though we don't need it, the point-to-process mapping is broadcast to all
	int* proc_points = check_malloc(sizeof(int)*(num_sgt_readers+1));
	if (debug) write_log("Receiving SGT point-to-process mapping from master.");
	check_bcast(proc_points, num_sgt_readers+1, MPI_INT, 0, MPI_COMM_WORLD, "Error broadcasting SGT point-to-process mapping to all, aborting.", 0);

	free(proc_points);
}
