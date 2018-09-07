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
	Read in list of RSQSim ruptures - 1 task per rupture
	Service requests for next task
	When no more tasks, send messages to processors that request
	When sent no more task messages to all processors, send message to processors 0 to N-1
	exit
 */

void broadcast_completion(int num_sgt_handlers, int num_workers, int num_procs, int my_id);
void manager_listen(int num_workers, worker_task* task_list, int num_tasks, int my_id);
int handle_work_request(manager_msg msg, worker_task* task_list, int* current_task, int num_tasks, MPI_Datatype worker_msg_type, int my_id);
int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int my_id);
void get_point_mapping(int num_sgt_readers);

int task_manager(int num_sgt_handlers, int num_workers, int num_procs, long long MAX_BUFFER_SIZE, int my_id) {
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

	int num_tasks = parse_rupture_list(rup_list_file, &task_list, MAX_BUFFER_SIZE, dtout, sgtmast.nt, sgtmast.globnp, my_id);

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

int parse_rupture_list(char rup_list_file[256], worker_task** task_list, long long MAX_BUFFER_SIZE, float dtout, int nt, int globnp, int my_id) {
	int i, j;
	/*File has format
	* <RSQSim SRF>
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
	int num_slips, num_hypos, num_points, num_tasks;
	float mag;
	//num ruptures to process is different, because some might have already completed
	int num_ruptures_to_process = 0;
	num_tasks = 0;
	for (i=0; i<num_ruptures; i++) {
		fscanf(rup_list_in, "%s", rupture_file);

		//Determine source ID - for RSQSim ruptures, assume rupture ID = rv ID = 0
		//e42_rv8_1000_0_event1090675.srf
		//Need third to last, next to last
		char* tok, *last, *ntl, *ttl;
		tok = last = ntl = ttl = 0;
		strcpy(string, rupture_file);
		tok = strtok(string, "/_");
		while (tok!=NULL) {
			ttl = ntl;
			ntl = last;
			last = tok;
			tok = strtok(NULL, "/_");
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
		//Removed task size check, since each task is just 1 SRF
		int tasks_for_rupture = 1;
		//Do not permit more than this amount to be used
		long long MAX_ALLOWED = (long long)(1.7 * 1024.0 * 1024.0 * 1024.0);
		//For RSQSim events, just 1 slip and 1 hypo
		num_slips = 1;
		num_hypos = 1;
		int num_vars = num_slips * num_hypos;

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
			//rup_var_spacing is not used for RSQSim ruptures, set to 0
			(*task_list)[num_tasks].rup_var_spacing = 0;
			//Only 1 rupture variation each
			(*task_list)[num_tasks].starting_rv_id = 0;
			(*task_list)[num_tasks].ending_rv_id = 1;
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
